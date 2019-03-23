#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-11-5
#Modefied by Duan Zhongqu, 2018-11-27
#Extract both the fully unaligned contigs and partially unaligned contigs.
#Extract the coordinate of partially unaligned contigs locate in reference genome.
package unalnCtg;

sub getUnaln{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_p);
getopts("hp:");

my $usage="\nUsage: hupan getUnalnCtg [options]  <assembly_directory> <QUAST_assess_directory> <output_directory>

hupan getUnalnCtg is used to collect unaligned contigs.

Necessary input description:

  assembly_directory      <string>    This directory should contain many sub-directories
                                      named by sample names, such as sample1, sample2,etc.
                                      In each sub-directory, assembly results, including 
                                      file *.contig, should exist.

  QUAST_assess_directory  <string>    This directory should contain many sub-directories 
                                      named by sample names, such as sample1, sample2,etc.
                                      In each sub-directory, quast assessment, including 
                                      directory file contigs_reports, should exist.

  output_directory        <string>    Results will be output to this directory.To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

Options:
     -h                               Print this usage page.

     -p                               The suffix of contigs file in assembly directory.
               

";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($contig_dir,$quast_dir,$out_dir)=@ARGV;

#detect perl script to collect contigs for single sample
my $execp="getUnalnCtg.pl";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
    $p.="/".$execp;
    if(-e $p && -x $p){
	$fpflag=1;
	last;
    }
}
die("Executable getUnalnCtg.pl cannot be found in your PATH!\n
") unless($fpflag);

#Adjust input directories
$contig_dir.="/" unless($contig_dir=~/\/$/);
$quast_dir.="/" unless($quast_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);
#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

my $thread_num=1;
my $suffix=".contig.gz";
$suffix=$opt_p if defined $opt_p;
#Create output directory and sub-directories
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#Read samples
opendir(DATA,$contig_dir) || die("Error03: cannot open data directory: $contig_dir\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    print STDERR "Process sample $s ...\n";
    my $cond=$contig_dir.$s."/";
    my $quad=$quast_dir.$s."/";
    unless(-d $cond && -d $quad){
	print STDERR "Warning: cannot find $s in $cond or $quad. Not processing.\n";
	next;
    }

#obtain *.contig file within the assembly directory
    my $contig_file="";
    opendir(ASS,$cond) || die("Error04: cannot open data directory: $cond\n");
    my @ass_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@ass_files){
	if($f=~/$suffix$/){
	    $contig_file=$f;
	    last;
	}
    }
    if($contig_file eq ""){
	print STDERR "Warnings: cannot find assembly file(*$suffix) in $cond\n";
	next; 
    }
    $contig_file=$cond.$contig_file;

#obtain the unaligned info file within quast dir
    $quad.="contigs_reports/";
    my $unalign_info_file="";
    opendir(ASS,$quad) || die("Error04: cannot open data directory: $quad\n");
    my @q_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@q_files){
	if($f=~/\.unaligned\.info$/){
	    $unalign_info_file=$f;
	    last;
	}
    }
    if($unalign_info_file eq ""){
	print STDERR "Warnings: cannot find unaligned info file(*.unaligned.info) in $quad\n";
	next; 
    }
    $unalign_info_file=$quad.$unalign_info_file;

#obtain the filtered coords file within quast dir
    $quad.="nucmer_output/";
    my $coords_file="";
    opendir(ASS,$quad) || die("Error04: cannot open data directory: $quad\n");
    my @nucmer_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@nucmer_files){
        if($f=~/\.coords\.filtered$/){
            $coords_file=$f;
            last;
        }
    }
    if($coords_file eq ""){
        print STDERR "Warnings: cannot find the filtered coords file(*.coords.filtered) in $quad\n";
        next;
    }
    $coords_file=$quad.$coords_file;

#generate command
    my $out_data_dir=$out_data.$s."/";
    mkdir($out_data_dir);
    my $out_prefix=$out_data_dir.$s;
    my $com="$execp $contig_file $unalign_info_file $coords_file $out_prefix $s";
    system($com);
}

my $out_total_dir=$out_dir."total/";
mkdir($out_total_dir);
mergeUnaln($out_data,$out_total_dir);
1;
}

sub mergeUnaln{
use strict;
use warnings;

my ($dir,$out)=@_;

die("$dir doesn't exist!
") unless -e $dir;
die("$dir is not a directory!
") unless -d $dir;

my $total_full=$out."total.fully.fa";
my $total_part=$out."total.partilly.fa";
my $total_coords=$out."total.partilly.coords";

$dir.="/" unless $dir=~/\/$/;
$out.="/" unless $out=~/\/$/;
opendir(DIR,$dir);
my @samples=readdir(DIR);
closedir DIR;
my $com1="cat";
my $com2="cat";
my $com3="cat";
foreach my $sample (@samples){
    next if $sample=~/^\.+$/;
    my $data_dir=$dir.$sample."/";
    opendir(DATA,$data_dir) || die("Error03: cannot open data directory: $data_dir\n");
    my @files=readdir(DATA);
    closedir DATA;
    my $full_unalign_file="";
    my $part_unalign_file="";
    my $part_coords_file="";   
    foreach my $file (@files){
        next if $sample=~/^\.+$/;
        $full_unalign_file=$data_dir.$file if $file=~/\.fully\.contig/;
        $part_unalign_file=$data_dir.$file if $file=~/\.partially\.contig/;
        $part_coords_file=$data_dir.$file if $file=~/\.partially\.coords/;
    }
    if($full_unalign_file eq ""){
        print STDERR "Warnings: cannot find fully unaligned file(*.fully.contig) from sample: $sample\n";
        next;
    }
    if($part_unalign_file eq ""){
        print STDERR "Warnings: cannot find partilly unaligned file(*.partilly.contig) from sample: $sample\n";
        next;
    }
    if($part_coords_file eq ""){
        print STDERR "Warnings: cannot find partilly unaligned coords(*..partially.coords) from sample: $sample\n";
        next;
    }
    $com1.=" $full_unalign_file";
    $com2.=" $part_unalign_file";
    $com3.=" $part_coords_file";
}

$com1.=" >$total_full";
$com2.=" >$total_part";
$com3.=" >$total_coords";
system($com1);
system($com2);
system($com3);
}
1;
