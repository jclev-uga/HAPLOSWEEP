#!/usr/bin/env perl

use strict;
use warnings;
use threads;
no strict qw(subs refs);

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case pass_through);
use threads;
use Data::Dumper;

my $usage = <<__EOUSAGE__;

######################################################################################################################
#                                                                                                                    #
#  PIPELINE PATH:                                                                                                    #
#  --genotype_parents   Set to identify new markers between two parents                                              #
#    --> REQUIRED:                                                                                                   #
#       -vcf <string>   vcf file from SNP caller/or tab delimeted file of chromosome name\possible variant position  #
#       -b <string>     sorted/indexed bam file of alignments                                                        #
#       -r <int>        read length (default 100bp)                                                                  #
#       --single/or     only calls haplotypes from within individual reads                                           #
#       --longrange     also calls long range haplotypes using paired read information (for large insert sizes)      #
#       -o <string>     prefix name for output files                                                                 #
#       --polyploid/or  calls haplotypes with allopolyploid/heterozygous diploid in mind or inbred diploid           #
#       --diploid                                                                                                    #
#       IF DIPLOID:                                                                                                  #
#       -m <int>        maximum number of bases per haplotype to consider (optional: default is 5)                   #
#                                                                                                                    #
#  --call_population    Set to use haplotypes called to genotype a population of genotypes                           #
#    --> REQUIRED:                                                                                                   #
#       -gen <string>   .gen file of haplotype positions from --genotype_parents path                                #
#       -progeny <string> file that lists names of progeny bam files to genotype                                     #       
#                                                                                                                    #
#  EXAMPLE USAGE:                                                                                                    #
#                                                                                                                    #
#  perl Haplotyper.pl --genotype_parents --polyploid --longrange -vcf snps.vcf -b gen1.sorted.bam\                   #
#           -b gen2.sorted.bam -r 150 -o example                                                                     #
#                                                                                                                    #
#  THEN TO GENOTYPE A POPULATION WITH IDENTIFIED HAPLOTYPES:                                                         #
#                                                                                                                    #
#  perl Haplotyper.pl --call_population -gen example.gen -progeny progeny_files.txt                                  # 
#                                                                                                                    #
######################################################################################################################
######################################################################################################################


__EOUSAGE__

    ;

my @bamfile=();
my $vcf = '';
my $handle;
my $help_flag;
my $genotype_parents = 0;
my $call_population = 0;
my $polyploid = 0;
my $single = 0;
my $longrange = 0;
my $read_length = 100;
my $genfile = '';
my $progeny = '';
my $max = 5;
my $output = '';
my $diploid = 0;

&GetOptions ( 'h' => \$help_flag,
              'b=s@' => \@bamfile,
              'vcf=s' => \$vcf,
              'r=i' => \$read_length,
              'single' => \$single,
              'o=s' => \$output,
              'longrange' => \$longrange,
              'polyploid' => \$polyploid,
              'gen=s' => \$genfile,
              'progeny=s' => \$progeny,
              'genotype_parents' => \$genotype_parents,
              'call_population' => \$call_population,
              'm=i' => \$max,
              'diploid' => \$diploid
              );

if (@ARGV) {
    die "Error, don't understand arguments: @ARGV ";
}

if ($help_flag) { die $usage; }

open (STDERR, ">&STDOUT");  ## capturing stderr and stdout in a single stdout stream

main: {


     my $input_bam = join(' ', @bamfile);
     my $length = 0+@bamfile;
     my $out_hap = join($output, ".hap");

#     &process_cmd("export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:./");     

     ##Path to identify new polymorphic haplotypes between parents
     if ($genotype_parents) {

          if ($polyploid) {

               if ($single) {
                    print "Running haplotype caller:\n";
                    &process_cmd("./HAPLOSWEEP $vcf \"$input_bam\" $length $read_length $output > $out_hap");
               } 
               elsif ($longrange) {
              #      print "Running long range haplotype caller:\n");
                    &process_cmd("./HAPLOSWEEP_LONGRANGE $vcf \"$input_bam\" $length $read_length $output > $out_hap");
               }
          } 
          elsif ($diploid) {

               print "Running diploid haplotype caller:\n";
               &process_cmd("./HAPLOSWEEP_DIPLOID $vcf \"$input_bam\" $length $read_length $output $max > $out_hap");
          }
     } 
     elsif ($call_population) {

          print "Genotyping population:\n";
          &process_cmd("python CallPopulation.py $genfile $progeny $output");
     }
     
     exit(0);

}


####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $start_time = time();
    my $ret = system($cmd);
    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";

    return;
}

