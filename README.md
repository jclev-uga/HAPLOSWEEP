
# HAPLOSWEEP v1.0
# By Josh Clevenger, PhD
# Project funded by a USDA-NIFA postdoctoral fellowship grant 2016-04731

Haplotype-based genotyping designed for polyploid species

Three C++ programs that can be run individually or within a perl pipeline.

Uses Bamtools API to process aligned reads in bam files so in order to run the pipeline first type the following command in the directory where the pipeline will run:

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./

Make sure the bamtools libraries libbamtools.so and libbamtools.so.2.4.1 are in the working directory as well.

Binary files are included, but C++ code can be compiled by user if bamtools API is installed.  C++ standard 11 was used.
.cpp files can be compiled with the following command:

g++ -I ./lib/bamtools-master/include/ -L ./lib/bamtools-master/lib/ -o program program.cpp -lz -lbamtools -std=c++11

If the above is true then running the pipeline (Haplotyper.pl) is simple.

                                                                                                                  
#  PIPELINE PATH:                                                                                                    
#  --genotype_parents   Set to identify new markers between two parents                                              
#    --> REQUIRED:                                                                                                   
#       -vcf <string>   vcf file from SNP caller or tab delimeted file of chromosome name possible variant position 
#       -b <string>     sorted/indexed bam file of alignments                                                        
#       -r <int>        read length (default 100bp)                                                                  
#       --single/or     only calls haplotypes from within individual reads                                           
#       --longrange     also calls long range haplotypes using paired read information (for large insert sizes)      
#       -o <string>     prefix name for output files                                                                 
#       --polyploid/or  calls haplotypes with allopolyploid/heterozygous diploid in mind or inbred diploid           
#       --diploid                                                                                                    
#       IF DIPLOID:                                                                                                  
#       -m <int>        maximum number of bases per haplotype to consider (optional: default is 5)                   
                                                                                                                    
#  --call_population    Set to use haplotypes called to genotype a population of genotypes                           
#    --> REQUIRED:                                                                                                   
#       -gen <string>   .gen file of haplotype positions from --genotype_parents path                                
#       -progeny <string> file that lists names of progeny bam files to genotype                                     
                                                                                                                    
#  EXAMPLE USAGE:                                                                                                    
                                                                                                                   
#  perl Haplotyper.pl --genotype_parents --polyploid --longrange -vcf snps.vcf -b gen1.sorted.bam -b gen2.sorted.bam -r 150 -o example                                                                     
                                                                                                                   
#  THEN TO GENOTYPE A POPULATION WITH IDENTIFIED HAPLOTYPES:                                                         
                                                                                                                   
#  perl Haplotyper.pl --call_population -gen example.gen -progeny progeny_files.txt                                  


If you encounter the error: 
"error while loading shared libraries: libbamtools.so.2.4.1: cannot open shared object file: No such file or directory"
You must set your LD_LIBRARY_PATH so the bamtools library can be found.
Type the following command in the shell and HaploSWEEP should execute:

        export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:./
        
Module commands for stand alone usage:

./HAPLOSWEEP vcffile "gen1.bam gen2.bam" num_bams read_length out > out.hap
./HAPLOSWEEP_LONGRANGE vcffile "gen1.bam gen2.bam" num_bams read_length out > out.hap

