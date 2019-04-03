################################################################################
# @Author: Ayush T Raman
# Aryee Lab, Broad Institute of MIT and Harvard
# Date: March 7th, 2019
#
# Program is used for:
# 1. post processed file per read
################################################################################

#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

## Working directory
my $file =  $ARGV[0]; ## file name
if($file =~ m/(.*)methylation_calls.tsv/){
    my $sample = $1;
#    print $sample,"\n";
    &postProcess($file, $sample);
}

sub postProcess{
    my @args = @_;
    my $input_file = $args[0];
    my $file = $args[1];
    my $result_file = $file."read-methylation-calls.tsv";
    my $meth_calls;
    my %result;
    my $read = "";
    my $count = 0;

    ## reading the file
    print "\n";
    print "Reading the input file: $input_file\n";
    open my $FILE,  "<",  $input_file or die "Cannot open the $file : $! \n";
    my $header = <$FILE>;
    while(my $line = <$FILE>){
        chomp $line;
        my @lines = split("\t", $line);
        my $chr = $lines[0];
        my $strand = $lines[1];
        my $read = $lines[4];
        my $ratio = $lines[5];
        my $call = 0 if($ratio < -2.5);
        $call = 1 if($ratio > 2.5);
        $call = "-" if($ratio >= -2.5 && $ratio <= 2.5);
        my $num_cg = $lines[9];
        $call = $call x $num_cg;

        ## read, strand and strand
        if(exists $result{$chr}{$strand}{$read}){
            my $start = $result{$chr}{$strand}{$read}{1};
            my $end = $result{$chr}{$strand}{$read}{2};
            $meth_calls = $result{$chr}{$strand}{$read}{3};
            $meth_calls = $meth_calls.$call;
            $result{$chr}{$strand}{$read}{3} = $meth_calls;
            $result{$chr}{$strand}{$read}{1} = $lines[2] if($lines[2] < $start);
            $result{$chr}{$strand}{$read}{2} = $lines[3] if($lines[3] > $end);
        }elsif(!exists $result{$chr}{$strand}{$read}){
            $result{$chr}{$strand}{$read}{1} = $lines[2];
            $result{$chr}{$strand}{$read}{2} = $lines[3];
            $result{$chr}{$strand}{$read}{3} = $call;
        }
    }
    close($FILE);

    ## writing file
    print "Output file: $result_file\n";
    open my $OUT,  ">",  $result_file or die "Cannot open the file: $! \n";
    print $OUT "chr\tstrand\tstart\tend\tread_name\taggregate_meth_calls\n";

    ## outputting the results in sorted order
    #print Dumper(\%result);
    foreach my $key1 (sort {lc $a cmp lc $b}  keys  %result){  ## $result{$chr}{$strand}{$read}
        foreach my $key2 (sort keys  %{$result{$key1}}){
            foreach my $key3 (keys %{$result{$key1}{$key2}}) {
                my $coord = $result{$key1}{$key2}{$key3}{1}."\t".$result{$key1}{$key2}{$key3}{2};
                my $geno = "chr".$key1."\t".$key2."\t".$coord;
                my $calls = $result{$key1}{$key2}{$key3}{3};
                print $OUT $geno,"\t",$key3,"\t","$calls","\n";
#                print $geno,"\t",$key3,"\t","$calls","\n";
            }
        }
    }
    close($OUT);
    print "\nDone writing the file: $result_file \n\n";
    print "**********  \n";
}
