import sys
import pysam
from collections import defaultdict
from collections import OrderedDict

##usage: python CallPopulation.py parentfile.gen progenyfile.txt output_handle

class DefaultListOrderedDict(OrderedDict):
    def __missing__(self,k):
        self[k] = []
        return self[k]

def main():
     parentfile = open(sys.argv[1], "rU")
     progenyfile = open(sys.argv[2], "rU")
     handle = sys.argv[3]

     outfile = open(handle+".gen", "w")
     outfile2 = open(handle+".score", "w")

     counter = 1
     genotype_gen = DefaultListOrderedDict()
     genotype_score = DefaultListOrderedDict()
     progeny_array = []

##loop through progeny
     for file in progenyfile:

          file=file.strip();
          progeny_array.append(file)
          samfile = pysam.Samfile(file, "rb")

##loop through markers
          for line in parentfile:
               if line.startswith('Chr'):
                    continue

               line=line.strip()
               line=line.split('\t')

               chr=line[0]
               pos1=int(line[1])
               pos2=int(line[2])
               gen1=line[3]
               gen2=line[4]

               result=haplotype(samfile, chr, pos1, pos2)

               if not result:
                    genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append("NA")
                    genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("0")
                    continue

               if gen1 not in result and gen2 not in result:
                    genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append("NA")
                    genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("0")
                    continue

               if gen1 in result and gen2 in result:
                    genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append("FALSE")
                    genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("0")
                    continue

               for hap in result:
                    if hap == gen1:
                         if result[hap] > 0:
                              genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append(hap)
                              genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("1")
                              continue
                         else:
                              genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append("NA")
                              genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("0")
                              continue
                    if hap == gen2:
                         if result[hap] > 0:
                              genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append(hap)
                              genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("3")
                              continue
                         else:
                              genotype_gen[chr+":"+str(pos1)+":"+str(pos2)].append("NA")
                              genotype_score[chr+":"+str(pos1)+":"+str(pos2)].append("0")


          parentfile.seek(0)


     parentfile.close()

     outfile.write("Chr\tPosition1\tPosition2\t")
     outfile2.write("Chr\tPosition1\tPosition2\t")

     for progeny in progeny_array:
          outfile.write("%s\t" % (progeny))
          outfile2.write("%s\t" % (progeny))

     outfile.write("\n")
     outfile2.write("\n")


     for marker in genotype_gen:

          outfile.write("%s\t%s\t%s\t" % (marker.split(":")[0],marker.split(":")[1],marker.split(":")[2]))
          outfile2.write("%s\t%s\t%s\t" % (marker.split(":")[0],marker.split(":")[1],marker.split(":")[2]))

          for i in range(len(genotype_gen[marker])):

                outfile.write("%s\t" % genotype_gen[marker][i])
                outfile2.write("%s\t" % genotype_score[marker][i])

          outfile.write('\n')
          outfile2.write('\n')


     outfile.close()
     outfile2.close()

def haplotype(samfile, chr, position1, position2):

     haplotypes = {}
     base1 = ''
     base2 = ''
     contrast = 0
     occurence = 0
     if position1 -2 >= 0:
          for pileupcolumn in samfile.pileup(chr, position1-2, position1+2):
               if (pileupcolumn.pos + 1) == position1:
                    for pileupread in pileupcolumn.pileups:
                         if not pileupread.query_position:
                              continue
#                         base1 = pileupread.alignment.seq[pileupread.qpos]
                         base1 = pileupread.alignment.query_sequence[pileupread.query_position]
                         if 0 <= (position2-position1+pileupread.query_position) < len(pileupread.alignment.query_sequence):
                              base2 = pileupread.alignment.query_sequence[position2-position1+pileupread.query_position]
                              if (base1+base2) not in haplotypes:
                                   haplotypes[base1+base2] = 1
                              else:
                                   haplotypes[base1+base2] += 1
     return haplotypes


main()

