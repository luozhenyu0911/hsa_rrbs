#!/usr/bin/env python3
"""
Prepare metilene input matrix from two groups of sorted bedGraph files.
Supports parallel processing by chromosome.
Automatically skips chromosomes: chrUn_*, *_random, *_alt, *_fix.
Filtering rule: keep a site only if BOTH groups have non-missing ratio >= --ratio.
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip

def pretty_time():
    from datetime import datetime
    return datetime.now().strftime("%a %b %d, %H:%M:%S, %Y")


def log(msg, level="INFO"):
    print(f"[{level}]\t{pretty_time()}\t{msg}", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate metilene input from two groups of sorted bedGraph files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Notes:
  - Input bedGraph files MUST be sorted (bedtools sort -i file.bg > file.sorted.bg)
  - Prefer --in1_list / --in2_list when you have many samples
  - Chromosomes matching chrUn_*, *_random, *_alt, *_fix are automatically skipped
  - A site is kept only if BOTH group1 and group2 have non-missing ratio >= --ratio
        """
    )

    # Input
    parser.add_argument("--in1", help="Comma-separated group1 bedGraph files")
    parser.add_argument("--in2", help="Comma-separated group2 bedGraph files")
    parser.add_argument("--in1_list", help="File with group1 paths (one per line)")
    parser.add_argument("--in2_list", help="File with group2 paths (one per line)")

    # Output & groups
    parser.add_argument("--out", help="Output file (default: metilene_<h1>_<h2>.input)")
    parser.add_argument("--h1", default="g1", help="Group1 name (default: g1)")
    parser.add_argument("--h2", default="g2", help="Group2 name (default: g2)")

    # Filtering & performance
    parser.add_argument("--ratio", type=float, default=0.0,
                        help="Min non-missing ratio required in EACH group (0-1, default: 0)")
    parser.add_argument("--NA", default=".", help="Missing value symbol (default: .)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Parallel jobs (default: 4; set 1 for single-pass)")
    parser.add_argument("-b", "--bedtools", default="bedtools",
                        help="Path to bedtools")
    parser.add_argument("--tmpdir", default="./tmp/",
                        help="Directory for temporary files (default: system temp)")

    args = parser.parse_args()

    if not ((args.in1 or args.in1_list) and (args.in2 or args.in2_list)):
        parser.error("Must provide --in1/--in2 or --in1_list/--in2_list")
    if not (0 <= args.ratio <= 1):
        parser.error("--ratio must be between 0 and 1")
    if args.threads < 1:
        parser.error("--threads must be >= 1")

    return args


def read_file_list(path):
    files = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                files.append(line)
    if not files:
        raise ValueError(f"Empty file list: {path}")
    return files


def check_files(files):
    resolved = []
    for f in files:
        p = Path(f).resolve()
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
        if not os.access(p, os.R_OK):
            raise PermissionError(f"File not readable: {p}")
        resolved.append(str(p))
    return resolved


def is_valid_chrom(chrom):
    """Skip chrUn_*, *_random, *_alt, *_fix."""
    if chrom.startswith("chrUn_"):
        return False
    if chrom.endswith(("_random", "_alt", "_fix")):
        return False
    return True


def get_chromosomes(bedgraph):
    chroms = set()
    with gzip.open(bedgraph, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom and is_valid_chrom(chrom):
                chroms.add(chrom)

    def key(c):
        m = re.search(r"(\d+)", c)
        return (0, int(m.group(1)), c) if m else (1, 0, c)

    return sorted(chroms, key=key)


def run_unionbedg(bedtools, names, na, file_list, awk_filter, outfile):
    """Run bedtools unionbedg + cut/sed/awk pipeline."""
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".txt") as lf:
        lf.write("\n".join(file_list) + "\n")
        listfile = lf.name

    cmd = (
        f"{bedtools} unionbedg -header -names {names} -filler '{na}' "
        f"-i $(cat {listfile}) "
        f"| cut -f1,3- "
        f"| sed 's/end/pos/' "
        f"| awk '{awk_filter}' "
        f"> {outfile}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    finally:
        os.unlink(listfile)


def filter_to_chrom(src, dest, chrom):
    """Extract lines of one chromosome."""
    with gzip.open(src, "rt") as fin, gzip.open(dest, "wt") as fout:
        for line in fin:
            if line.startswith(("#", "track", "browser")):
                continue
            if line.startswith(chrom + "\t"):
                fout.write(line)


def filter_valid_chroms(src, dest, valid_chroms):
    """Keep only lines from valid chromosomes."""
    with gzip.open(src, "rt") as fin, gzip.open(dest, "wt") as fout:
        for line in fin:
            if line.startswith(("#", "track", "browser")):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom in valid_chroms:
                fout.write(line)


def process_one_chrom(chrom, short_paths, bedtools, names, na, awk_filter, out_dir):
    """Worker: process a single chromosome."""
    work_dir = Path(tempfile.mkdtemp(prefix=f"chr_{chrom}_", dir=out_dir))
    try:
        chr_files = []
        for i, src in enumerate(short_paths):
            dest = work_dir / f"{i}.bg"
            filter_to_chrom(src, dest, chrom)
            chr_files.append(str(dest))

        chr_out = Path(out_dir) / f"{chrom}.out"
        run_unionbedg(bedtools, names, na, chr_files, awk_filter, str(chr_out))
        return chrom, str(chr_out)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def main():
    args = parse_args()

    log("Input files must be sorted (bedtools sort)", "WARNING")
    log("Checking input...")

    # Load groups
    g1 = read_file_list(args.in1_list) if args.in1_list else \
         [x.strip() for x in args.in1.split(",") if x.strip()]
    g2 = read_file_list(args.in2_list) if args.in2_list else \
         [x.strip() for x in args.in2.split(",") if x.strip()]

    g1 = check_files(g1)
    g2 = check_files(g2)

    n1 = len(g1)
    n2 = len(g2)

    names = " ".join([args.h1] * n1 + [args.h2] * n2)
    out_file = str(Path(args.out or f"metilene_{args.h1}_{args.h2}.input").resolve())

    if shutil.which(args.bedtools) is None and not Path(args.bedtools).exists():
        sys.exit(f"ERROR: bedtools not found: {args.bedtools}")

    log(f"group1 = {n1}, group2 = {n2}")
    log(f"ratio (per group) = {args.ratio}, threads = {args.threads}")
    log(f"Output: {out_file}")

    # Awk filter: keep site only if BOTH groups meet the non-missing ratio
    # After cut -f1,3- : $1=chrom, $2=pos, $3..(2+n1)=group1, then group2
    awk_filter = (
        f'NR==1{{print; next}} '
        f'{{ '
        f'  n1={n1}; n2={n2}; '
        f'  miss1=0; miss2=0; '
        f'  for(i=3; i<=2+n1; i++) if($i=="{args.NA}") miss1++; '
        f'  for(i=3+n1; i<=2+n1+n2; i++) if($i=="{args.NA}") miss2++; '
        f'  if( (n1-miss1)/n1 >= {args.ratio} && (n2-miss2)/n2 >= {args.ratio} ) print; '
        f'}}'
    )

    # Temporary directory
    if args.tmpdir:
        root = Path(args.tmpdir).resolve()
        root.mkdir(parents=True, exist_ok=True)
        tmpdir = Path(tempfile.mkdtemp(prefix="metilene_", dir=root))
        need_cleanup = True
    else:
        tmp_ctx = tempfile.TemporaryDirectory(prefix="metilene_")
        tmpdir = Path(tmp_ctx.name)
        need_cleanup = False

    try:
        log(f"Temp dir: {tmpdir}")

        # Short symlinks
        short_paths = []
        for i, f in enumerate(g1 + g2):
            link = tmpdir / f"{i}.bg"
            link.symlink_to(f)
            short_paths.append(str(link))

        # Always get valid chromosomes
        log("Collecting valid chromosomes (skip chrUn_*/_random/_alt/_fix)...")
        chroms = get_chromosomes(short_paths[0])
        valid_set = set(chroms)
        log(f"Valid chromosomes: {len(chroms)}")

        if args.threads == 1:
            # Single-pass with chromosome filter
            log("Single-pass mode (with chromosome filter)...")
            filtered = []
            for i, src in enumerate(short_paths):
                dest = tmpdir / f"filt_{i}.bg"
                filter_valid_chroms(src, dest, valid_set)
                filtered.append(str(dest))
            run_unionbedg(args.bedtools, names, args.NA, filtered, awk_filter, out_file)

        else:
            # Parallel by chromosome
            log(f"Parallel mode ({args.threads} threads)...")
            chr_dir = tmpdir / "chr_out"
            chr_dir.mkdir()

            results = {}
            with ProcessPoolExecutor(max_workers=args.threads) as pool:
                futures = {
                    pool.submit(process_one_chrom, chrom, short_paths,
                                args.bedtools, names, args.NA, awk_filter, str(chr_dir)): chrom
                    for chrom in chroms
                }
                for fut in as_completed(futures):
                    chrom, path = fut.result()
                    results[chrom] = path
                    log(f"Done: {chrom}")

            # Merge in order, keep one header
            log("Merging results...")
            header_done = False
            with open(out_file, "w") as out:
                for chrom in chroms:
                    path = results.get(chrom)
                    if not path or not Path(path).exists():
                        continue
                    with open(path) as fh:
                        for line in fh:
                            if line.startswith(("chrom\t", "chr\t")):
                                if not header_done:
                                    out.write(line)
                                    header_done = True
                                continue
                            out.write(line)

    finally:
        if need_cleanup and tmpdir.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)

    print("\n*****\n", file=sys.stderr)
    print(f"[BASIC CALL] metilene -t 4 -a {args.h1} -b {args.h2} {out_file} > out.file",
          file=sys.stderr)
    print("Adjust -t / -a / -b as needed.\n", file=sys.stderr)


if __name__ == "__main__":
    main()