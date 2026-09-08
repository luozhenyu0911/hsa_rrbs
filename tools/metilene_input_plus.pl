#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use File::Basename;
use File::Temp qw/tempdir tempfile/;

# -----------------------------------------------------------------------------
# Variables

my ($USAGE, $call, $ret, $help);
my ($in1, $in2, $in1_list, $in2_list, $out, $h1, $h2, $bedtools);
my ($NA, $out_metilene, $g1_h, $g2_h);
my ($ratio);
my (@g1, @g2);
my $SCRIPTNAME = basename($0);

# -----------------------------------------------------------------------------
# OPTIONS

$USAGE = << "USE";

    usage:  perl $SCRIPTNAME [options]

    === 输入文件（二选一） ===
        --in1 <list>        逗号分隔的 group1 bedGraph 文件（样本少时用）
        --in2 <list>        逗号分隔的 group2 bedGraph 文件（样本少时用）

        --in1_list <file>   group1 文件列表（每行一个路径，推荐！样本多时用这个）
        --in2_list <file>   group2 文件列表（每行一个路径，推荐！样本多时用这个）

    === 其他参数 ===
        --out <file>        输出文件（默认: metilene_<h1>_<h2>.input）
        --h1 <string>       group1 名称（默认: g1）
        --h2 <string>       group2 名称（默认: g2）
        --ratio <float>     保留位点的最小非缺失比例（0~1，默认: 0，即不过滤）
                            例: 0.75 表示至少 75% 样本有值才保留该位点
        --NA <string>       缺失值符号（默认: .）
        -b <path>           bedtools 路径（默认从 PATH 查找）
        -h|--help           显示帮助

    注意: 输入 bedGraph 必须已经排序！
          bedtools sort -i file.bg > file.sorted.bg

USE

if (!@ARGV) {
    printf STDERR $USAGE;
    exit -1;
}

unless (GetOptions(
    "in1=s"       => \$in1,
    "in2=s"       => \$in2,
    "in1_list=s"  => \$in1_list,
    "in2_list=s"  => \$in2_list,
    "out=s"       => \$out,
    "h1=s"        => \$h1,
    "h2=s"        => \$h2,
    "ratio=f"     => \$ratio,
    "NA=s"        => \$NA,
    "b=s"         => \$bedtools,
    "h|help"      => \$help
)){
    printf STDERR $USAGE;
    exit -1;
}
if (defined $help){
    printf STDERR $USAGE;
    exit -1;
}

# -----------------------------------------------------------------------------
# MAIN

############
## checks ##
############
print STDERR ("[WARNING]" . prettyTime() . "Input files need to be SORTED, i.e., use \"bedtools sort -i file >file.sorted\"\n");
print STDERR ("[INFO]" . prettyTime() . "Checking flags\n");

# ----- group 1 -----
if (defined $in1_list) {
    # 从文件列表读取（推荐，避免 ARG_MAX）
    $in1_list = File::Spec->rel2abs($in1_list);
    die "##### AN ERROR: --in1_list file not found: $in1_list\n" unless -e $in1_list;
    die "##### AN ERROR: --in1_list file not readable: $in1_list\n" unless -r $in1_list;

    open(my $fh, "<", $in1_list) or die "##### Cannot open $in1_list: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g;          # trim
        next if $line eq "" || $line =~ /^#/;  # skip empty / comment
        push @g1, $line;
    }
    close $fh;
    die "##### AN ERROR: --in1_list is empty\n" if scalar(@g1) == 0;
}
elsif (defined $in1) {
    @g1 = split(/,/, $in1);
}
else {
    die "##### AN ERROR: must provide either --in1 or --in1_list\n";
}

for (my $i = 0; $i < scalar(@g1); $i++){
    $g1[$i] = File::Spec->rel2abs($g1[$i]);
    if (-e $g1[$i]) {
        die "##### AN ERROR: $g1[$i] not readable\n" unless -r $g1[$i];
    } else {
        die "##### AN ERROR: file not found: $g1[$i]\n";
    }
}

# ----- group 2 -----
if (defined $in2_list) {
    $in2_list = File::Spec->rel2abs($in2_list);
    die "##### AN ERROR: --in2_list file not found: $in2_list\n" unless -e $in2_list;
    die "##### AN ERROR: --in2_list file not readable: $in2_list\n" unless -r $in2_list;

    open(my $fh, "<", $in2_list) or die "##### Cannot open $in2_list: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq "" || $line =~ /^#/;
        push @g2, $line;
    }
    close $fh;
    die "##### AN ERROR: --in2_list is empty\n" if scalar(@g2) == 0;
}
elsif (defined $in2) {
    @g2 = split(/,/, $in2);
}
else {
    die "##### AN ERROR: must provide either --in2 or --in2_list\n";
}

for (my $i = 0; $i < scalar(@g2); $i++){
    $g2[$i] = File::Spec->rel2abs($g2[$i]);
    if (-e $g2[$i]) {
        die "##### AN ERROR: $g2[$i] not readable\n" unless -r $g2[$i];
    } else {
        die "##### AN ERROR: file not found: $g2[$i]\n";
    }
}

# ----- defaults -----
if (!defined $h1) { $h1 = "g1"; }
if (!defined $h2) { $h2 = "g2"; }
if (!defined $NA) { $NA = "."; }
if (!defined $ratio) { $ratio = 0; }

if ($ratio < 0 || $ratio > 1) {
    die "##### AN ERROR: --ratio must be between 0 and 1 (got $ratio)\n";
}

# build header name strings
$g1_h = $h1;
$g2_h = $h2;
for (my $i = 1; $i < scalar(@g1); $i++){ $g1_h .= " $h1"; }
for (my $i = 1; $i < scalar(@g2); $i++){ $g2_h .= " $h2"; }

# output file
if (defined $out){
    $out_metilene = File::Spec->rel2abs($out);
} else {
    $out_metilene = "metilene_${h1}_${h2}.input";
}

## bedtools ##
if (defined $bedtools){
    $bedtools = File::Spec->rel2abs($bedtools);
    die "##### AN ERROR: bedtools not found: $bedtools\n" unless -e $bedtools;
    die "##### AN ERROR: bedtools is a directory\n" if -d $bedtools;
    die "##### AN ERROR: bedtools not executable\n" unless -x $bedtools;
} else {
    $bedtools = "bedtools";
}
$ret = system("command -v $bedtools &> /dev/null");
if ($ret != 0){
    die "##### AN ERROR: No bedtools found. Use -b to specify path.\n";
}

#####################
## union + filter  ##
#####################
my $total_samples = scalar(@g1) + scalar(@g2);
print STDERR ("[INFO]" . prettyTime() . "group1 = " . scalar(@g1) . " samples, group2 = " . scalar(@g2) . " samples\n");
print STDERR ("[INFO]" . prettyTime() . "Total samples = $total_samples  |  ratio threshold = $ratio\n");
print STDERR ("[INFO]" . prettyTime() . "Write metilene input to $out_metilene\n");

# Create short symlinks + file list to avoid ARG_MAX inside bedtools
my $tmpdir = tempdir(CLEANUP => 1);
print STDERR ("[INFO]" . prettyTime() . "Temporary directory: $tmpdir\n");

my @all_files = (@g1, @g2);
my @short_paths;
for (my $i = 0; $i < scalar(@all_files); $i++) {
    my $link = File::Spec->catfile($tmpdir, "$i.bg");
    symlink($all_files[$i], $link)
        or die "##### Cannot create symlink $link -> $all_files[$i]\n";
    push @short_paths, $link;
}

my $listfile = File::Spec->catfile($tmpdir, "filelist.txt");
open(my $lfh, ">", $listfile) or die "##### Cannot write $listfile: $!\n";
print $lfh join("\n", @short_paths), "\n";
close $lfh;

# awk filter
my $awk_code = 'NR==1{print;next}{non_miss=0;for(i=3;i<=NF;i++)if($i!="'.$NA.'")non_miss++;if(non_miss/(NF-2)>='.$ratio.')print}';

# Write pipeline to a temporary shell script (command line stays short)
my ($script_fh, $script_path) = tempfile(
    "metilene_unionXXXXXX",
    DIR    => $tmpdir,
    SUFFIX => ".sh",
    UNLINK => 1
);

print $script_fh "#!/bin/bash\n";
print $script_fh "set -euo pipefail\n\n";
print $script_fh "mapfile -t FILES < '$listfile'\n";
print $script_fh "$bedtools unionbedg -header -names $g1_h $g2_h -filler '$NA' -i \"\${FILES[@]}\" \\\n";
print $script_fh "  | cut -f1,3- \\\n";
print $script_fh "  | sed 's/end/pos/' \\\n";
print $script_fh "  | awk '$awk_code' \\\n";
print $script_fh "  > '$out_metilene'\n";

close $script_fh;
chmod 0755, $script_path;

print STDERR ("[INFO]" . prettyTime() . "Running unionbedg pipeline...\n");
$ret = system($script_path);
if ($ret != 0) {
    die "##### AN ERROR while running unionbedg pipeline (exit $ret)\n";
}

print STDERR ("\n*****\n\n"
    . "[BASIC CALL:] metilene -t 4 -a $h1 -b $h2 $out_metilene > out.file\n"
    . "Please adjust threads (-t) and group names (-a/-b) as needed.\n");


# -----------------------------------------------------------------------------
# FUNCTIONS

sub prettyTime {
    my @months   = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    return "\t$weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $hour:$minute:$second, $year\t";
}