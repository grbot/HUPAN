#!/usr/bin/perl
#Created by Duan Zhongqu, 2018-7-2
package rmHigh;
sub alignContig{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_q $opt_s $opt_t);
    getopts("hq:s:t:");
    my $usage="\nUsage: hupanSLURM align [options] <data_directory> <output_directory> <MUMmer_directory> <reference.fa>

align is used to align the assembly result to reference genome sequence by numcer in MUMmer progam.

The script will call MUMmer program, so you need to tell the program where MUMmer locates.

Necessary input description:

  data_directory       <string>      This directory should contain many sub-directories
                                     named by sample names, such as Sample1, Sample2,etc.
                                     In each sub-directory, assembly results, including
                                     files \"*.scafSeq\" and \"*.contig\", should exist.

  output_directory     <string>      Both final output files and intermediate results 
                                     will be found in this directory. To avoid 
                                     overwriting of existing files, we kindly request
                                     that the output_directory should not exist. It is
                                     to say, this directory will be created by the 
                                     script itself.

  <MUMmer_directory>   <string>      MUMmer directory where nucmer and show-coords locates.

  <reference.fa>       <string>      Reference sequence file (.fa or .fa.gz).

Options:
     -h                              Print this usage page.

     -q                <string>      The queue name for job submiting.
                                     Default: default queue
     
     -t                <int>         Threads used.
                                     Default: 1. 

     -s                <string>      Suffix of assembled file. 
                                     Defult: \"contigs.fa\"
";

    die $usage if @ARGV!=4;
    die $usage if defined($opt_h);
    my ($data_dir,$out_dir,$MUMmer_dir,$ref)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files, we kindly request that the
 output directory should not exist.
");
}

    #Detect executable nucmer and show-coords
    $MUMmer_dir.="/" unless($MUMmer_dir=~/\/$/);
    my $exec_nucmer=$MUMmer_dir."nucmer";
    die("Error01: Cannot find nucmer file in directory $MUMmer_dir\n") unless(-e $exec_nucmer);
    my $exec_showcoords=$MUMmer_dir."show-coords";
    die("Error01: Cannot find show-coords file in directory $MUMmer_dir\n") unless(-e $exec_showcoords);

    #Check existence of reference sequence file
    die("Error02: Cannot find reference sequence file\n") unless(-e $ref);

    #Adjust directory names and create output directory
    unless($data_dir=~/\/$/){
	$data_dir.="/";
    }

    unless($out_dir=~/\/$/){
	$out_dir.="/";
    }
    mkdir($out_dir);

    my $suffix="contigs.fa";
    $suffix=$opt_s if(defined $opt_s);

    #read threads
    my $thread_num=1;
    if(defined($opt_t)){
        $thread_num=$opt_t;
    }
     
    #create directory for alignment result
    my $align_out=$out_dir."data/";
    mkdir($align_out);
    print STDERR "Processing samples:\n";

#******* Might be modified for different task submission system *******
    my $job_out=$out_dir."job";
    mkdir($job_out);
    my $script_out=$job_out."/scripts";
    mkdir($script_out);
    my $stderr_out=$job_out."/err";
    mkdir($stderr_out);
    my $stdout_out=$job_out."/out";
    mkdir($stdout_out);
#**********************************************************************

    #Read samples
    opendir(DATA,$data_dir) || die("Error03: cannot open data directory: $data_dir\n");
    my @sample=readdir(DATA);
    closedir DATA;

##process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $sd=$data_dir.$s."/";
    next unless(-d $sd);
    print STDERR "Process sample $s ...\n";
    #obtain *_.contig file within the sample directory
    my $contig_file="";
    opendir(ASS,$sd) || die("Error04: cannot open data directory: $sd\n");
    my @ass_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@ass_files){
        if($f=~/$suffix$/){
            $contig_file=$f;
            last;
        }
    }
    if($contig_file eq ""){
        print STDERR "Warnings: cannot find assembly file in $sd\n";
        next;
    }
    $contig_file=$sd.$contig_file;
    #create output directory for a sample
    my $sample_out=$align_out.$s."/";
    mkdir($sample_out) unless(-e $sample_out);
    my $prefix=$sample_out.$s;
    my $delta=$prefix.".delta";
    my $coords=$prefix.".coords";
    #generate command
    my $com;
    $com="$exec_nucmer --prefix=$prefix $ref $contig_file\n";
    $com.="$exec_showcoords -cl $delta >$coords";
#******* Might be modified for different task submission system *******
    my $job_file=$script_out."/".$s.".slurm";   #script_file
    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
    print JOB "\#!/bin/bash\n\n";
    print JOB "\#SBATCH -J $s","_align\n";              #job name
    print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;  #queue name in the submission system
    print JOB "\#SBATCH --output=$out_file\n";               #stdout
    print JOB "\#SBATCH --error=$err_file\n";               #stderr
    print JOB "\#SBATCH -n $thread_num\n";                  #thread number
    print JOB "\#SBATCH --ntasks-per-node=$thread_num\n\n";   #use the same node
    print JOB "\#SBATCH --time=96:00:00\n"; # set time
    print JOB "$com\n";                              #commands
    close JOB;
    system("sbatch $job_file");                       #submit job
#*********************************************************************
}
1;
}
sub extractSeq{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_i $opt_c $opt_s);
    getopts("hi:c:s:");
    my $usage="\nUsage: hupanSLURM extractSeq [options] <data_dir> <out_dir> <MUMmer_output_dir>

extractSeq is used to extract contigs that is lower similarity with reference genomes.

Necessary input description:

data_dir             <string>      This directory should contain many sub-directories
                                   named by sample names, such as Sample1, Sample2,etc.
                                   In each sub-directory, assembly results, including
                                   files *.scafSeq and *.contigs.fa, should exist.

out_dir              <string>      Results will be output to this directory.
                                   To avoid overwriting of existing files, we kindly request
                                   that the output_directory should not exist. It is
                                   to say, this directory will be created by the 
                                   script itself.

MUMmer_output_dir    <string>      This directory should contain many sub-directories
                                   named by sample names, such as Sample1, Sample2,etc.
                                   In each sub-directory, alignment results, including
                                   files *.coords, should exist.

Options:
     -h                            Print this usage page.

     -i              <float>       The theshold of identity. 
                                   Default: 0.95.

     -c              <float>       The theshold of query coverage.
                                   Default: 0.95.

     -s              <string>      Suffix of assembled file. 
                                   Defult: \"contigs.fa\" 
";
    die $usage if @ARGV!=3;
    die $usage if defined($opt_h);
    my ($contig_dir,$out_dir,$coords_dir)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n");
}

    #get identity
    my $identity=0.95;
    if(defined($opt_i)){
        $identity=$opt_i;
    }

    #get coverage
    my $coverage=0.95;
    if(defined($opt_c)){
        $coverage=$opt_c;
    }

    my $suffix="contigs.fa";
    $suffix=$opt_s if(defined $opt_s);

    #Adjust directory names and create output directory
    unless($contig_dir=~/\/$/){
	$contig_dir.="/";
    }

    unless($coords_dir=~/\/$/){
        $coords_dir.="/";
    }
    
    unless($out_dir=~/\/$/){
	$out_dir.="/";
    }
    mkdir($out_dir);

    #create directory for extractSeq result
    my $extractSeq_out=$out_dir."data/";
    mkdir($extractSeq_out);
    print STDERR "Processing samples:\n";

    #Read samples
    opendir(DATA,$contig_dir) || die("Error03: cannot open data directory: $contig_dir");
    my @sample=readdir(DATA);
    closedir DATA;

##process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $sd=$contig_dir.$s."/";
    next unless(-d $sd);
    print STDERR "Process sample $s ...\n";
    #obtain *_.contig file within the sample directory
    my $contig_file="";
    opendir(ASS,$sd) || die("Error04: cannot open data directory: $sd\n");
    my @ass_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@ass_files){
        if($f=~/$suffix$/){
            $contig_file=$f;
            last;
        }
    }
    if($contig_file eq ""){
        print STDERR "Warnings: cannot find assembly file in $sd\n";
        next;
    }
    $contig_file=$sd.$contig_file;
    #create output directory for a sample
    my $sample_out=$extractSeq_out.$s."/";
    mkdir($sample_out) unless(-e $sample_out);
    #Read the contigs name from assembly result
    open(FILE1,$contig_file);
    my %names;
    while(<FILE1>){
    chomp;
    if($_=~/^>/){
        my @t=split /\s+/,$_;
        my $name=substr($t[0],1,length($t[0])-1);
        $names{$name}=1;
        }
    }
    close FILE1;
    my $length1=keys %names;
    print "There are :".$length1." contigs in the contigs file.\n";    

    #get the coords file
    my $coords_file=$coords_dir."data/".$s."/".$s.".coords";
    
    #Read the contigs name of high similarity with reference genome from coords file
    open(FILE2,$coords_file)||die("Can't open the file: ".$coords_file);
    my %contigs;
    my $rows=1;
    while(<FILE2>){
        if($rows<6){
            $rows++;
        }
        else{
            chomp;
            $_=~s/^ +//;
            my @t=split /\s+/,$_;
            #print @t;
            my $i=$t[9];
            my $c=$t[15];
            my $contig=$t[18];
            if(($i>=$identity*100)&&($c>=$coverage*100)){
                #print $i."\t".$c."\t"."\t".$contig,"\n";
                #push @contigs, $contig;
                $contigs{$contig}=1
            }
        }
    }
    close FILE2;
    my $length2=keys %contigs;
    print "There are: ".$length2." contigs highly similarity with the reference genome with >= ".($identity*100)."% identity and >= ".($coverage*100)."% coverage.\n";

    #Remove the contigs name of high similarity with reference genome
    foreach my $key (keys %names){
        if(exists $contigs{$key}){
            delete $names{$key};
        }
    }
    my $length3=keys %names;
    print "There are :".$length3." candidate unaligned contigs.\n";

    #Collect the candidat unaligned contigs from contigs file
    my $candidate_unalign_file=$sample_out.$s.".candidate.unaligned.contig";
    open(FILE3,$contig_file)||die("Error9: cannot read file:$contig_file\n");
    open(FILE4,">$candidate_unalign_file")||die("Error10: cannot write file:$candidate_unalign_file\n");
    my $flag=0;
    while(my $line=<FILE3>){
        chomp $line;
        if($line=~/^>/){
            my @t=split /\s+/,$line;
            my $name=substr($t[0],1,length($t[0])-1);
            if(exists $names{$name}){
                print FILE4 $line."\n";
                $flag=1;
            }
            else{
                $flag=0;
            }
        }
        else{
            if($flag==1){
                print FILE4 $line."\n";            
            }
        }
    }
    close FILE3;
    close FILE4;     
}
1;
}
1;
