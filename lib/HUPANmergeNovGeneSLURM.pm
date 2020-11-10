#!/usr/bin/perl
##Created by Duan Zhongqu, 2018-10-24
package mergeNovGene;

sub mergeNovGene{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_p $opt_j);
    getopts("hp:j:");
    my $usage="\nhupanSLURM mergeNovGen [options] <data_directory> <output_directory> <maker_directory>

hupanSLURM mergeNovGen is used to merge maker result from multiple files.

Necessary input description:

  <data_directory>          <string>     The directory includes maker results (*.all.maker.gff, 
                                         *.all.maker.proteins.fasta and *.all.maker.transcripts.fasta) 
                                         from multiple small files.

  <output_directory>        <string>     The output directory.

  <maker_directory>         <string>     The directory maker locates.

Options:
  
          -h                             Print this usage.

          -p                <string>     Define the registered genome prefix.
                                         default: NovGene

          -j                <int>        The length of the number following the prefix.
                                         default: 6

";

    die $usage if @ARGV!=3;
    die $usage if $opt_h;
    my ($data_dir,$out_dir,$maker_dir)=@ARGV;

#check existense of output directory
    if(-e $out_dir){
        die("Error: output directory \"$out_dir\" already exists.\nTo avoid overwriting of existing files, we kindly request that the \noutput directory should not exist.\n");
    }

#adjust directory names and create output directory
    $data_dir.="/" unless($data_dir=~/\/$/);
    $out_dir.="/" unless($out_dir=~/\/$/);
    $maker_dir.="/" unless($maker_dir=~/\/$/);
    mkdir($out_dir);

    my $prefix="NovGene";
    $prefix=$opt_p if defined $opt_p;
    my $justify=6;
    $justify=$opt_j if defined $opt_j;
    
#get executable file
    my $gff_merge=$maker_dir."gff3_merge";
    my $maker_map_ids=$maker_dir."maker_map_ids";
    my $map_gff_ids=$maker_dir."map_gff_ids";
    my $map_fasta_ids=$maker_dir."map_fasta_ids";

    my $gff_file="";
    my $all_prot_file=$out_dir."combine.all.maker.proteins.fasta";
    open(PROT,">$all_prot_file")||die("Error10: cannot write the file: $all_prot_file.\n");
    my $all_tran_file=$out_dir."combine.all.maker.transcripts.fasta";
    open(TRAN,">$all_tran_file")||die("Error10: cannot write the file: $all_tran_file.\n");
#read sample
    opendir(DATA,$data_dir) || die("Error: can not open the data directory: $data_dir!\n")
;
    my @files=readdir(DATA);
    foreach my $file (@files){
        next if $file=~/^\./;
        if($file=~/\.all\.maker\.gff/){
            my $tmp_file=$data_dir.$file;
            $gff_file.=$tmp_file." ";
        }
        elsif($file=~/\.all\.maker\.proteins\.fasta/){
            my $tmp_file=$data_dir.$file;
            open(IN,$tmp_file)||die("Error09: cannot read the file: $tmp_file.\n");
            while(my $line=<IN>){
                print PROT $line;
            }
            close IN;
        }
        elsif($file=~/\.all\.maker\.transcripts\.fasta/){
            my $tmp_file=$data_dir.$file;
            open(IN,$tmp_file)||die("Error09: cannot read the file: $tmp_file.\n");
            while(my $line=<IN>){
                print TRAN $line;
            }
            close IN;
        }
        else{
            print STDERR "Warning: $file isn't a gff file or fasta file from maker result! => Not processed.\n";
        }
    }
    close PROT;
    close TRAN;
    my $all_gff_file=$out_dir."combine.all.maker.gff";
    my $com="$gff_merge -o $all_gff_file $gff_file\n";
    my $map_file=$out_dir."all.maker.map";
    $com.="$maker_map_ids --prefix $prefix --justify $justify $all_gff_file > $map_file\n";
    $com.="$map_gff_ids $map_file $all_gff_file\n";
    $com.="$map_fasta_ids $map_file $all_tran_file\n";
    $com.="$map_fasta_ids $map_file $all_prot_file\n";
    system("$com");   
}
1;
