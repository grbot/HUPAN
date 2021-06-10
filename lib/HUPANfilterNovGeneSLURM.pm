#!/usr/bin/perl
#Created by Duan Zhongqu at 2018-11-5
package filterNovGen;

sub filterNovGen{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_s $opt_t $opt_i);
    getopts("hs:t:i:");
    my $usage="\nUsage: hupan filterNovGen [potions] <maker_result_dir> <output_dir> <reference_dir> <blast_dir> <cd-hit_dir> <RepeatMask_dir>

eupan filterNovGen is used to filter the novel precited genes

Necessary descriptions:
    
    maker_result_dir        <string>        The directory contains maker result. The novel sequence file
                                            (*.novel.fa), gff file (*.gff), protein sequence file 
                                            (*.maker.proteins.fasta) and transcription sequence file
                                            (*.maker.transcripts.fasta) should be included.

    output_dir              <string>        The output_directory.

    reference_dir           <string>        The reference genome sequence and transcripts sequence 
                                            file 'ref.genome.fa' and 'ref.transcripts.fa' should be included.

    blast_dir               <string>        The directory of blastn and makeblastdb locates.
    
    cd-hit_dir              <string>        The directory of cd-hit-est locates.
    
    RepeatMask_dir          <string>        The directory of RepeatMask locates.

Options:
  
  -h                                        Print this Usage page.

  -i                        <float>         The identity used in removing homology gene sequences.
                                            Default: 0.8
  
  -s                        <float>         The identity with referenece genome or transription.
                                            Default: 0.5

  -t                        <int>           The thread number.
                                            Default: 1
";

    die $usage if @ARGV<6;
    die $usage if $opt_h;
    my ($data_dir,$out_dir,$ref_dir,$blast_dir,$cdhit_dir,$Repeatmask_dir)=@ARGV;

#check existense of output directory
    if(-e $out_dir){
        die("Error: output directory \"$out_dir\" already exists.\nTo avoid overwriting of existing files, we kindly request that the \noutput directory should not exist.\n");
    }

#adjust directory names and create output directory
    $data_dir.="/" unless($data_dir=~/\/$/);
    $out_dir.="/" unless($out_dir=~/\/$/);
    mkdir($out_dir);
    $ref_dir.="/" unless($ref_dir=~/\/$/);
    $blast_dir.="/" unless($ref_dir=~/\/$/);
    $cdhit_dir.="/" unless($Repeatmask_dir=~/\/$/);
    $Repeatmask_dir.="/" unless($Repeatmask_dir=~/\/$/);

    my $gene_identity=0.8;
    $gene_identity=$opt_i if defined $opt_i;
    my $ref_identity=0.5;
    $ref_identity=$opt_s if defined $opt_s;
    my $thread_num=1;
    $thread_num=$opt_t if defined $opt_t;
    
#check the input file
    opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
    my @files=readdir(DATA);
    closedir DATA;
    my $seq_file;
    my $gff_file;
    my $tra_file;
    my $pro_file;
    foreach my $file (@files){
        next if $file=~/^\./;
        if($file=~/novel\.fa$/){
            $seq_file=$data_dir.$file;
        }
        elsif($file=~/gff$/){
            $gff_file=$data_dir.$file;
        }
        elsif($file=~/maker\.transcripts\.fasta$/){
            $tra_file=$data_dir.$file;
        }
        elsif($file=~/maker\.proteins\.fasta$/){
            $pro_file=$data_dir.$file;
        }
        else{
            print STDERR "Warning: $file is not used in subsequent analysis.\n";
        }
    }
    die("The gff file doesn't exist!\n") unless -e $gff_file;
    die("The sequences file doesn't exist!") unless -e $seq_file;
    die("The transcripts file doesn't exist!") unless -e $tra_file;
    die("The proteins file doesn't exist!") unless -e $pro_file;

#check the reference file
    my $genome_ref=$ref_dir."ref.genome.fa";
    die("$genome_ref doesn't exist! Please rename the reference genome file as 'ref.genome.fa'.\n") unless -e $genome_ref;
    my $tran_ref=$ref_dir."ref.transcripts.fa";
    die("$tran_ref doesn't exist! Please rename the reference genome file as 'ref.genome.fa'.\n") unless -e $tran_ref;

#check the binary file
    my $cdhit_exe=$cdhit_dir."cd-hit-est";
    die(" $cdhit_exe doesn't exist!\n") unless -e $cdhit_exe;
    my $blast_exe=$blast_dir."blastn";
    die(" $blast_exe doesn't exist!\n") unless -e $cdhit_exe;
    my $mkblastdb_exe=$blast_dir."makeblastdb";
    die(" $mkblastdb_exe doesn't exist!\n") unless -e $cdhit_exe;
    my $repeat_exe=$Repeatmask_dir."RepeatMasker";
    die(" $repeat_exe doesn't exist!\n") unless -e $repeat_exe;

    my $exefg="filterNovelgene.pl";
    my @path=split /:/,$ENV{PATH};
    my $fpflag=0;
    foreach my $p (@path){
        $p.="/".$exefg;
        if(-e $p && -x $p){
            $fpflag=1;
            last;
        }
    }
    die("Executable filterNovelgene cannot be found in your PATH!\n") unless($fpflag);

    my $data_out=$out_dir."data/";
    mkdir($data_out);
#************** Might be modified for different task submission system *******************
    my $job_out=$out_dir."job";       
    mkdir($job_out);
    my $script_out=$job_out."/scripts"; #job script directory
    mkdir($script_out);
    my $stderr_out=$job_out."/err";     #stdout directory
    mkdir($stderr_out);
    my $stdout_out=$job_out."/out";     #sdterr directory
    mkdir($stdout_out);
##****************************************************************************************

    my $com="$exefg $seq_file $gff_file $tra_file $pro_file $genome_ref $tran_ref $cdhit_exe $blast_exe $mkblastdb_exe $repeat_exe $gene_identity $ref_identity $thread_num $data_dir\n";

#generate and submit job script
#************** Might be modified for different task submission system *******************
    my $job_file=$script_out."/FilterNovelGene.slurm";   #script_file
    my $err_file=$stderr_out."/FilterNovelGene.err";   #stderr_output_file
    my $out_file=$stdout_out."/FilterNovelGene.out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error05: Unable to create job file: $job_file\n");
    print JOB "\#!/bin/bash\n\n";
    print JOB "\#SBATCH --job-name=FilterNovelGene\n";              #job name
    print JOB "\#SBATCH --output=$out_file\n";               #stdout
    print JOB "\#SBATCH --error=$err_file\n";               #stderr
    print JOB "\#SBATCH -n $thread_num\n";                  #thread number
    print JOB "\#SBATCH --ntasks-per-node=$thread_num\n\n";   #use the same node
    print JOB "$com\n";                              #commands
    close JOB;
    system("sbatch $job_file");                       #submit job
#************** Might be modified for different task submission system *******************
}
1;
