#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-9-5
#Modefied by Duan Zhongqu, 2018-7-2
#Added the function of de novo assembly by sga with lower memory.
#Modefied by Duan Zhongqu, 2020-4-25
#Fix some bugs, for example, the suffix of sequencing data files.

package assembly;

sub assemble{
    my $usage="
hupanSLURM assemble [commands] ...

Commands:
\tsoapdenovo    Aseembly with SOAPdenovo 2.
\tlinearK       Assembly with an iterative use of SOAPdenovo 2 (Recommended).
\tsga           Assembly with SGA (Recommended.)
";
    die $usage if @ARGV<1;
    my $com=shift @ARGV;
    if($com eq "soapdenovo"){
	soap(@ARGV);
    }
    elsif($com eq "linearK"){
	linearK(@ARGV);
    }
    elsif($com eq "sga"){
        sga(@ARGV);
    }
    else{
	print STDERR "Unknown command: $com\n";
	die($usage);
    }
}

sub soap{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_s $opt_n $opt_m $opt_k $opt_s $opt_c $opt_g $opt_q);
getopts("hs:n:m:t:k:s:c:gq:");

my $usage="\nUsage: hupanSLURM assemble soapdenovo [options] <fastq_data_directory> <output_directory> <soapdenovo_directory>

hupanSLURM assemble soapdenovo is used to assemble high-quality reads on large scale.

Necessary input description:

  fastq_data_directory    <string>    This directory should contain many sub-directories
                                      named by sample names, such as Sample1, Sample2,etc.
                                      In each sub-directory, there should be several 
                                      sequencing files ended by \".fq.gz\" or \".fastq.gz\".

  output_directory        <string>    Alignment results will be output to this directory.
                                      To avoid overwriting of existing files, we kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  sopadenovo_directory    <string>    Directory where soapdenovo2 executable files exists   

Options:
     -h                               Print this usage page.

     -t                   <int>       Threads used.
                                      Default: 1

     -s                   <string>    Suffix of files within data_directory.
                                      Default: \".fastq.gz\" 

     -k                   <int>       Kmer.
                                      Default: 35

     -c                   <string>    Parameters of soapdenovo2 config file. 8 parameters ligated by comma
                                        1)maximal read length
                                        2)average insert size
                                        3)if sequence needs to be reversed
                                        4)in which part(s) the reads are used
                                        5)use only first N bps of each read
                                        6)in which order the reads are used while scaffolding
                                        7)cutoff of pair number for a reliable connection (at least 3 for 
                                          short insert size)
                                        8)minimum aligned length to contigs for a reliable read location 
                                          (at least 32 for short insert size)
                                      Default: 80,460,0,3,80,1,3,32

     -g                               Enable gapcloser 

     -q                  <string>     The queue name for job submiting. 
                                      Default: default queue
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$tool_dir)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n");
}

$tool_dir.="/" unless $tool_dir=~/\/$/;
my $exec63=$tool_dir."SOAPdenovo-63mer";
die("Cannot find $exec63 in directory $tool_dir\n") unless -e $exec63;
my $exec127=$tool_dir."SOAPdenovo-127mer";
die("Cannot find $exec127 in directory $tool_dir\n") unless -e $exec127;
my $execgap=$tool_dir."GapCloser";
if(defined $opt_g){
die("Cannot find $execgap in directory $tool_dir\n") unless -e $execgap;
}
#read threads
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}
#define kmer
my $kmer=35;
$kmer=$opt_k if defined $opt_k;

#define file suffix
my $suffix=".fastq.gz";
$suffix=$opt_s if defined($opt_s);

my ($max_rd_len,$avg_ins,$reverse_seq,$asm_flags,$rd_len_cutoff,$rank,$pair_num_cutoff,$map_len)
    =(80,460,0,3,80,1,3,32);
($max_rd_len,$avg_ins,$reverse_seq,$asm_flags,$rd_len_cutoff,$rank,$pair_num_cutoff,$map_len)
    =split /,/,$opt_c if defined($opt_c);

#define the min length of contigs 
my $min_contig_len=500;
$min_contig_len=$opt_n if defined($opt_n);


#Adjust directory names and create output directory

$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#read samples
opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
my @sample=readdir(DATA);
closedir DATA;

#************** Might be modified for different task submission system *******************
my $job_out=$out_dir."job";       
mkdir($job_out);
my $script_out=$job_out."/scripts"; #job script directory
mkdir($script_out);
my $stderr_out=$job_out."/err";     #stdout directory
mkdir($stderr_out);
my $stdout_out=$job_out."/out";     #sdterr directory
mkdir($stdout_out);
#*****************************************************************************************

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $s_dir=$data_dir.$s."/";
    my $o_dir=$out_data.$s."/";
    mkdir($o_dir);
    opendir(DATA,$s_dir) || die("Error: can not read input data directory: $s_dir\n");
    my @files=readdir(DATA);
    closedir DATA;

#parse fastq pairs
    my %fq_base;
    foreach my $f (@files){
	next if $f=~/^\.+$/;
	next if $f=~/^single/;
	print STDERR "Warnig: $f without suffix: $suffix\n" unless $f=~/$suffix$/;
	next unless $f=~/$suffix$/;
	my $fb=substr($f,0,length($f)-length($suffix)-1);
	$fq_base{$fb}=1 unless defined($fq_base{$fb});
    }

#generate SOAPdenovo config
    my $config_file=$o_dir."soap.config";
    open(CONFIG,">$config_file") || die("Error02: can not write soap config file:$config_file\n");
    print CONFIG "\#maximal read length\nmax_rd_len=$max_rd_len\n[LIB]\n\#average insert size\navg_ins=$avg_ins\n\#if sequence needs to be reversed\nreverse_seq=$reverse_seq\n\#in which part(s) the reads are used\nasm_flags=$asm_flags\n\#use only first 100 bps of each read\nrd_len_cutoff=$rd_len_cutoff\n\#in which order the reads are used while scaffolding\nrank=$rank\n\# cutoff of pair number for a reliable connection (at least 3 for short insert size)\npair_num_cutoff=$pair_num_cutoff\n\#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=$map_len\n\#a pair of fastq file, read 1 file should always be followed by read 2 file\n";
    foreach my $b (keys(%fq_base)){
	my $forward=$s_dir.$b."1".$suffix;
	my $reverse=$s_dir.$b."2".$suffix;
	print STDERR "Warning: missed file: $forward\n" unless -e $forward;
	print STDERR "Warning: missed file: $reverse\n" unless -e $reverse;
	print STDERR "Warning: missed file: $reverse\n" unless(-e $forward && -e $reverse);
	print CONFIG "q1=$forward\nq2=$reverse\n" if(-e $forward && -e $reverse);
    }
    close CONFIG;
    my $com;
    if($kmer<=63){
        $com=$exec63;
    }
    else{
	$com=$exec127;
    }
    $com.=" all -s $config_file -o $o_dir/K$kmer -K $kmer -R -F -p $thread_num >$o_dir/soap.log 2>&1\n";
    
#run gapcloser
    if(defined $opt_g){
	$com.=$execgap;
	$com.=" -b $config_file -a $o_dir/K$kmer.scafSeq -o $o_dir/K\Q$kmer\E.gcScafSeq -t $thread_num >$o_dir/gapcloser.log 2>&1\n";
	$com.="breakGCscaf.pl $o_dir/K\Q$kmer\E.gcScafSeq $o_dir/K\Q$kmer\E.gcContig\n";
    }

#************** Might be modified for different task submission system *******************
    my $job_file=$script_out."/".$s.".slurm";   #script_file
    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
	    print JOB "\#!/bin/bash\n";
    print JOB "\#SBATCH --job-name=$s\n";                          #job name
	    print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;   #queue name in the submission system
    print JOB "\#SBATCH --output=$out_file\n";                   #stdout
    print JOB "\#SBATCH --error=$err_file\n";                   #stderr
    print JOB "\#SBATCH -n $thread_num\n";                 #thread number
    print JOB "\#SBATCH --ntasks-per-node=$thread_num\n"; #use the same node
    print JOB "$com\n";                                  #commands
    close JOB;
    system("sbatch $job_file");                           #submit job
#*****************************************************************************************
}
1;
}

sub linearK{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_s $opt_g $opt_r $opt_w $opt_u $opt_c $opt_n $opt_k $opt_m $opt_q);
getopts("ht:g:s:r:w:u:c:n:k:m:q:");

my $usage="\nUsage: hupanSLURM assmble linearK [options] <fastq_data_directory> <output_directory> <soapdenovo_directory>

hupanSLURM assmble linearK is used to assemble high-quality reads on large scale.

Necessary input description:

  fastq_data_directory    <string>    This directory should contain many sub-directories
                                      named by sample names, such as Sample1, Sample2,etc.
                                      In each sub-directory, there should be several 
                                      sequencing files ended by .fastq or .fastq.gz.

  output_directory        <string>    Alignment results will be output to this directory.
                                      To avoid overwriting of existing files, we kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  sopadenovo_directory    <string>    Directory where soapdenovo2 executable files exists   

Options:
     -h                               Print this usage page.

     -t                   <int>       Threads used.
                                      Default: 1

     -g                   <int>       Genome size. Used to infer sequencing depth. 
                                      Default: 3000000000 (3G)
     
     -s                   <string>    Suffix of files within data_directory.
                                      Default: \".fastq.gz\" 

     -r                   <string>    Parameters of linear function: Kmer=2*int(0.5*(a*Depth+b))+1. 
                                      The parameter should be input as \"a,b\".
                                      Default: 0.76,20

     -w                   <int>       Step-length of Kmer change.
                                      Default: 2

     -u                   <int>       Upper limmited times of Kmer change. This parameter is set to reduce
                                      redundancy computation.
                                      Default: 10

     -c                   <string>    Parameters of soapdenovo2 config file. 8 parameters ligated by comma
                                        1)maximal read length
                                        2)average insert size
                                        3)if sequence needs to be reversed
                                        4)in which part(s) the reads are used
                                        5)use only first N bps of each read
                                        6)in which order the reads are used while scaffolding
                                        7)cutoff of pair number for a reliable connection (at least 3 for 
                                          short insert size)
                                        8)minimum aligned length to contigs for a reliable read location 
                                          (at least 32 for short insert size)
                                      Default: 80,460,0,3,80,1,3,32

     -n                   <int>       The minimum length of contigs. Contigs shorter than this length will
                                      NOT be used when calculating N50.
                                      Default: 100

     -k                   <string>    Available Kmer range. Give comma-seperated lower bound and upper bound.
                                      Default: 15,127

     -m                   <int>       The number of consecutive Ns to be broken down to contigs.This is used 
                                      in the process break gapclosed scaffolds to contigs.
                                      Default: 10.

     -q                   <string>    The queue name for job submiting. 
                                      Default: default queue
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$soapdenovo)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n");
}
#Check executable linearK.pl
my $exec="linearK.pl";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
  $p.="/".$exec;
  if(-e $p && -x $p){
     $fpflag=1;
	last;
  }
}
die("Executable linearK.pl cannot be found in your PATH!\n
") unless($fpflag);


#read threads
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}

#Adjust directory names and create output directory

$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#read samples
opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
my @sample=readdir(DATA);
closedir DATA;

#************** Might be modified for different task submission system *******************
my $job_out=$out_dir."job";       
mkdir($job_out);
my $script_out=$job_out."/scripts"; #job script directory
mkdir($script_out);
my $stderr_out=$job_out."/err";     #stdout directory
mkdir($stderr_out);
my $stdout_out=$job_out."/out";     #sdterr directory
mkdir($stdout_out);
#*****************************************************************************************

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $com="linearK.pl -t $thread_num";
    $com.=" -g $opt_g" if defined $opt_g;
    $com.=" -s $opt_s" if defined $opt_s;
    $com.=" -r $opt_r" if defined $opt_r;
    $com.=" -w $opt_w" if defined $opt_w;
    $com.=" -u $opt_u" if defined $opt_u;
    $com.=" -c $opt_c" if defined $opt_c;
    $com.=" -n $opt_n" if defined $opt_n;
    $com.=" -k $opt_k" if defined $opt_k;
    $com.=" -m $opt_m" if defined $opt_m;
    my $indir=$data_dir.$s;
    my $outdir=$out_data.$s;
    $com.=" $indir $outdir $soapdenovo";

#************** Might be modified for different task submission system *******************
    my $job_file=$script_out."/".$s.".slurm";   #script_file
    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
	    print JOB "\#!/bin/bash\n";
    print JOB "\#SBATCH --job-name=$s\n";                          #job name
print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;   #queue name in the submission system
    print JOB "\#SBATCH --output=$out_file\n";                   #stdout
    print JOB "\#SBATCH --error=$err_file\n";                   #stderr
    print JOB "\#SBATCH -n $thread_num\n";                 #thread number
    print JOB "\#SBATCH --ntasks-per-node=$thread_num\n"; #use the same node
    print JOB "$com\n";                                  #commands
    close JOB;
    system("sbatch $job_file");                           #submit job
#*****************************************************************************************
}
1;
}

1;

sub sga{
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_d $opt_q $opt_m $opt_s);
getopts("ht:d:q:m:s:");

my $usage="\nUsage: hupan assemble sga [options] <fastq_data_directory> <output_directory> <sga_directory>

hupan assemble sga is used to assemble high-quality reads on large scale.

Necessary input description:

  fastq_data_directory    <string>    This directory should contain many sub-directories
                                      named by sample names, such as Sample1, Sample2,etc.
                                      In each sub-directory, there should be several 
                                      sequencing files ended by \".fastq.gz\" or \".fq.gz\".

  output_directory        <string>    Alignment results will be output to this directory.
                                      To avoid overwriting of existing files, we kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  sga_directory           <string>    Directory where sga executable files exists   

Options:
     -h                               Print this usage page.

     -t                   <int>       Threads used.
                                      Default: 1.

     -s                   <string>    Suffix of files within data_directory.
                                      Default: \".fastq.gz\"  

     -d                   <int>       The parameter sets used for different depths of sequencing
                                      data, we obtained two optimatial paramter sets from 
                                      simulated data for 30-fold and 60-fold, respectively.
                                      Default: 30.

     -q                   <string>    The queue name for job submiting.
                                      Default: default queue

     -m                   <int>       The intermediate results are huge and we kindly suggested 
                                      delete them after finishing each steps. If you want to keep
                                      them, please set this parameter as 1.  0: delete; 1: keep.
                                      Default: 0.
";
die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$sga_dir)=@ARGV;

$data_dir=abs_path($data_dir);
$out_dir=abs_path($out_dir);
#print $data_dir."\n";
#print $out_dir."\n";
my $cmd_dir=$ENV{'PWD'};

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n");
}

#Check executable sga
$sga_dir.="/" unless $sga_dir=~/\/$/;
my $sga=$sga_dir."sga";
die("Cannot find $sga in directory $sga_dir\n") unless -e $sga;

#read threads
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}

#define depth
my $depth=30;
$depth=$opt_d if defined $opt_d;

#define file suffix
my $suffix=".fastq.gz";
$suffix=$opt_s if defined($opt_s);

#whether keep the intermediate results
my $keep=0;
$keep=$opt_m if defined $opt_m;

$data_dir.="/" unless $data_dir=~/\/$/;
$out_dir.="/" unless $out_dir=~/\/$/;
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#read samples
opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $s_dir=$data_dir.$s."/";
    print STDERR "Process sample $s ...\n";
    my $o_dir=$out_data.$s."/";
    mkdir($o_dir);
    opendir(DATA,$s_dir) || die("Error: can not read input data directory: $s_dir\n");
    my @files=readdir(DATA);
    closedir DATA;
    
    #parse fastq pairs
    my @fq_base;
    foreach my $f (@files){
        next if $f=~/^\.+$/;
        next if $f=~/^single/;
        print STDERR "Warnig: $f without suffix: $suffix\n" unless $f=~/$suffix$/;
        next unless $f=~/$suffix$/;
        my $fb=substr($f,0,length($f)-length($suffix)-1);
        #print $fb."\n";
        push @fq_base, $fb;
    }
    my %count;
    my @uniq_times=grep { ++$count{ $_ } < 2; } @fq_base;
    @fq_base=@uniq_times;
    
    my $com="";
    my $fastq1;
    my $fastq2;
    my $preprocess_fastq;
    my $length=@fq_base;
    if($length==1){
        my $b=$fq_base[0];
        $preprocess_fastq=$o_dir.$s.".fastq";
        $fastq1=$s_dir.$b."1".$suffix;
        $fastq2=$s_dir.$b."2".$suffix;
        die("Error: missed file: $fastq1\n") unless -e $fastq1;
        die("Error: missed file: $fastq2\n") unless -e $fastq2;
        $com="$sga preprocess -o $preprocess_fastq --pe-mode 1 $fastq1 $fastq2\n";
        $com.="$sga index -a ropebwt --no-reverse -t $thread_num $preprocess_fastq\n";
    }
    else{
        my @list;
        my $merge_prefix=$o_dir.$s;
        foreach my $b (@fq_base){
            $fastq1=$s_dir.$b."1".$suffix;
            $fastq2=$s_dir.$b."2".$suffix;
            die("Error: missed file: $fastq1\n") unless -e $fastq1;
            die("Error: missed file: $fastq2\n") unless -e $fastq2;
            $preprocess_fastq=$s_dir.$s.".fastq";
            $com="$sga preprocess -o $preprocess_fastq --pe-mode 1 $fastq1 $fastq2\n";
            $com.="$sga index -a ropebwt --no-reverse -t $thread_num $preprocess_fastq\n";
            push @list, $preprocess_fastq;
        }
        $com.="$sga merge -p $merge_prefix @list\n";
    }
    my $correct_fastq=$o_dir.$s.".correct.fastq";
    my $filter_fa=$o_dir.$s.".correct.filter.pass.fa";
    my $merge_fa=$o_dir.$s.".correct.filter.pass.merged.fa";
    my $rmdup_fa=$o_dir.$s.".correct.filter.pass.merged.rmdup.fa";
    my $asqg=$o_dir.$s.".correct.filter.pass.merged.rmdup.asqg.gz";
    my $assemble_prefix=$o_dir.$s.".assemble";
    $com.="$sga correct -k 55 --learn -t $thread_num -o $correct_fastq $preprocess_fastq\n";
    if($keep==0){
        my $sai_file=$o_dir.$s.".sai";
        my $bwt_file=$o_dir.$s.".bwt";
        $com.="rm $preprocess_fastq $sai_file $bwt_file\n";
    }
    $com.="$sga index -a ropebwt -t $thread_num $correct_fastq\n";
    $com.="$sga filter -x 2 -o $filter_fa -t $thread_num $correct_fastq\n";
    if($keep==0){
        my $sai_file=$o_dir.$s.".correct.sai";
        my $bwt_file=$o_dir.$s.".correct.bwt";
        my $rsai_file=$o_dir.$s.".correct.rsai";
        my $rbwt_file=$o_dir.$s.".correct.rbwt";
        $com.="rm $correct_fastq $sai_file $bwt_file $rsai_file $rbwt_file\n";
    }
    if($depth==30){
        $com.="$sga fm-merge -m 65 -t $thread_num -o $merge_fa $filter_fa\n";
        if($keep==0){
            my $sai_file=$o_dir.$s.".correct.filter.pass.sai";
            my $bwt_file=$o_dir.$s.".correct.filter.pass.bwt";
            my $rsai_file=$o_dir.$s.".correct.filter.pass.rsai";
            my $rbwt_file=$o_dir.$s.".correct.filter.pass.rbwt";
            my $discard_fa=$o_dir.$s.".correct.filter.pass.discard.fa";
            $com.="rm $filter_fa $sai_file $bwt_file $rsai_file $rbwt_file $discard_fa\n";
        }
        $com.="$sga index -d 20000000 -t $thread_num $merge_fa\n";
        $com.="$sga rmdup -t $thread_num -o $rmdup_fa $merge_fa\n";
        if($keep==0){
            my $sai_file=$o_dir.$s.".correct.filter.pass.merged.sai";
            my $bwt_file=$o_dir.$s.".correct.filter.pass.merged.bwt";
            my $rsai_file=$o_dir.$s.".correct.filter.pass.merged.rsai";
            my $rbwt_file=$o_dir.$s.".correct.filter.pass.merged.rbwt";
            $com.="rm $merge_fa $sai_file $bwt_file $rsai_file $rbwt_file\n";
        }
        $com.="$sga overlap -m 65 -t $thread_num $rmdup_fa\n";
        if($keep==0){
            my $sai_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.sai";
            my $bwt_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.bwt";
            my $rsai_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.rsai";
            my $rbwt_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.rbwt";
            my $dup_fa=$o_dir.$s.".correct.filter.pass.merged.rmdup.dups.fa";
            $com.="rm $rmdup_fa $sai_file $bwt_file $rsai_file $rbwt_file $dup_fa\n";
        }
        $com.="$sga assemble -m 91 -l 160 -o $assemble_prefix $asqg\n";
    }  
    elsif($depth==60){
        $com.="$sga fm-merge -m 85 -t $thread_num -o $merge_fa $filter_fa\n";
        if($keep==0){
            my $sai_file=$o_dir.$s.".correct.filter.pass.sai";
            my $bwt_file=$o_dir.$s.".correct.filter.pass.bwt";
            my $rsai_file=$o_dir.$s.".correct.filter.pass.rsai";
            my $rbwt_file=$o_dir.$s.".correct.filter.pass.rbwt";
            my $discard_fa=$o_dir.$s.".correct.filter.pass.discard.fa";
            $com.="rm $filter_fa $sai_file $bwt_file $rsai_file $rbwt_file $discard_fa\n";
        }
        $com.="$sga index -d 20000000 -t $thread_num $merge_fa\n";
        $com.="$sga rmdup -t $thread_num -o $rmdup_fa $merge_fa\n";
        if($keep==0){
            my $sai_file=$o_dir.$s.".correct.filter.pass.merged.sai";
            my $bwt_file=$o_dir.$s.".correct.filter.pass.merged.bwt";
            my $rsai_file=$o_dir.$s.".correct.filter.pass.merged.rsai";
            my $rbwt_file=$o_dir.$s.".correct.filter.pass.merged.rbwt";
            $com.="rm $merge_fa $sai_file $bwt_file $rsai_file $rbwt_file\n";
        }
        $com.="$sga overlap -m 85 -t $thread_num $rmdup_fa\n";
        if($keep==0){
            my $sai_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.sai";
            my $bwt_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.bwt";
            my $rsai_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.rsai";
            my $rbwt_file=$o_dir.$s.".correct.filter.pass.merged.rmdup.rbwt";
            my $dup_fa=$o_dir.$s.".correct.filter.pass.merged.rmdup.dups.fa";
            $com.="rm $rmdup_fa $sai_file $bwt_file $rsai_file $rbwt_file $dup_fa\n";
        }
        $com.="$sga assemble -m 97 -d 0.4 -g 0.1 -r 50 -l 160 -o $assemble_prefix $asqg\n";
    }
    else{
        die("Please check the depth!");
    }
    
    #remove tmp file
    #$com.="remove *bai *sai *correct.fastq *correct.filter.pass.fa *correct.filter.pass.merged.fa *.correct.filter.pass.merged.rmdup.fa\n";
    #******* Might be modified for different task submission system *******
    my $job_out=$out_dir."job";       
    mkdir($job_out);
    my $script_out=$job_out."/scripts"; #job script directory
    mkdir($script_out);
    my $stderr_out=$job_out."/err";     #stdout directory
    mkdir($stderr_out);
    my $stdout_out=$job_out."/out";     #sdterr directory
    mkdir($stdout_out);
    my $job_file=$script_out."/".$s.".slurm";   #script_file
    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
    print JOB "\#!/bin/bash\n";
    print JOB "\#SBATCH --job-name=$s\n";                   #job name
    print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;     #queue name in the submission system
    print JOB "\#SBATCH --output=$out_file\n";              #stdout
    print JOB "\#SBATCH --error=$err_file\n";               #stderr
    print JOB "\#SBATCH -n $thread_num\n";                  #thread number
    print JOB "\#SBATCH --ntasks-per-node=$thread_num\n";   #use the same node
    print JOB "\#SBATCH --time=336:00:00\n"; # set time
    print JOB "$com\n";                                     #commands
    close JOB;    
    chdir "$o_dir";             
    system("sbatch $job_file");                             #submit job                              
    chdir "$cmd_dir";
  }
1;
}
1;

