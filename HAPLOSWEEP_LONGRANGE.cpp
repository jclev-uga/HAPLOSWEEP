#include <time.h>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"
#include "utils/bamtools_utilities.h"


using namespace std;
using namespace BamTools;

//struct object to collect paired reads over two positions

struct paired_reads {
         paired_reads() : Ref_ID(0), Read1_Position(0), Read1_Aligned(), Read2_Position(2), Read2_Aligned(), read1_haps(), read2_haps() {}
         paired_reads(int ref, int read1, std::string aligned1, int read2, std::string aligned2, std::map<string, string> haps1, std::map<string, string> haps2)
         : Ref_ID(ref), Read1_Position(read1), Read1_Aligned(aligned1), Read2_Position(read2), Read2_Aligned(aligned2), read1_haps(haps1), read2_haps(haps2) {}
	string Name;
	int Ref_ID;
	int Read1_Position;
	string Read1_Aligned;
        int Read2_Position;
        string Read2_Aligned;
        std::map<string, string> read1_haps;
        std::map<string, string> read2_haps;
};


//function to gather haplotypes from reads structure
//Takes a map of paired_reads objects and key of read name and returns a map of counts with haplotypes

map<string, int> get_haplotypes(std::map<std::string, paired_reads> reads, int position1, int position2, int read_length, map<string, string>& bases)
{
     string base1 = "";
     string base2 = "";
     string hap = "";
     map<string, int> haps;
     int index1 = 0;
     int index2 = 0;

     //loop through read pairs
     if (reads.empty())
          return {{"NN", 0}};

     for (const auto& x : reads ) {
        base1 = "";
        base2 = "";
        hap = "";

        index1 = (position1 - reads[x.first.c_str()].Read1_Position - 1);
        index2 = (position2 - reads[x.first.c_str()].Read1_Position - 1);
        //if only one of the read pairs disregard
        if ((index1 < read_length) and (index2 < read_length)) {


           //collect positional bases form aligned reads
           base1 = reads[x.first.c_str()].Read1_Aligned[index1];

           base2 = reads[x.first.c_str()].Read1_Aligned[index2];

           if (index2 > read_length)
              continue;
           if (base1 == "N" or base2 == "N" or (base2 != "A" and base2 != "G" and base2 != "C" and base2 != "T"))
              continue;


           hap = base1+base2;

           //if the haplotype is in the map add 1 to the count, else add the haplotype to the map
           if (haps.count(hap)) {
                haps[hap] += 1;
           }
           else {
                haps[hap] = 1;
           }
        }
     }
    if (haps.empty())
        return {{"NN", 0}};
    else
        return haps;
}

/* Function to collect haplotypes from paired reads taking as input
	bamfiles and locations and returning a dictionary of haplotypes map<string, int>& Haps */

map<string, int> collect_haplotype(BamReader& reader, int ref, int position1, int position2, int read_length, int bam_number, std::map<int, std::map<string, paired_reads> >& all_reads, map<string, string>& bases, map<int, string> reference_names)
{

        string location1 = "";
        string location2 = "";
        string base1 = "";
        string base2 = "";


        std::map<std::string, paired_reads> reads;
        string mate1 = "";

        //set the region to the first position
        reader.SetRegion(ref,position1,ref,position1);

        BamAlignment al;

        location1 = reference_names[ref] + ":" + to_string(position1);
        location2 = reference_names[ref] + ":" + to_string(position2);


//      Loop while there is a valid alignment overlapping with position 1
        while (reader.GetNextAlignment(al)) {
               int check_skip = 0;

//             Check if the read has a skip in it (exome and RNAseq data) 
               for (int i=0; i < al.CigarData.size(); i++ ) {
                     if (al.CigarData[i].Type == 'I') {
                            check_skip = 1;
                            break;
                     }
               }

//             if the read has a skip in it then move to the next one
               if (check_skip == 1) {
                    continue;
               }

//            Check if reads are not in region 
               if (abs(position1-al.Position) > read_length) {
                    continue;
               }

//             if the read is the first map
	       if (al.IsFirstMate() and abs(position1-al.Position) < read_length) {
//		    add to the paired_reads struct
                    reads[al.Name].Ref_ID = al.RefID;
                    reads[al.Name].Read1_Position = al.Position;
                    reads[al.Name].Read1_Aligned = al.AlignedBases;
                    all_reads[bam_number][al.Name].Ref_ID = al.RefID;
                    all_reads[bam_number][al.Name].Read1_Position = al.Position;
                    all_reads[bam_number][al.Name].Read1_Aligned = al.AlignedBases;
                    all_reads[bam_number][al.Name].read1_haps[to_string(position1)] = al.AlignedBases[position1 - al.Position - 1];

//                  store base unformation for later checking long range haplotypes
                    base1 = al.AlignedBases[position1-al.Position-1];
                    base2 = al.AlignedBases[position2-al.Position-1];
                    if (base1 == "A" or base1 == "T" or base1 == "G" or base1 == "C"){
                    bases[location1] += base1;
                    }
                    if (base2 == "A" or base2 == "T" or base2 == "G" or base2 == "C"){
                    bases[location2] += base2;
                    }
               }
//             if the read is actually the second mate
               else if (al.IsSecondMate() and abs(position1-al.Position) < read_length) {
                    //if it is not in the map add as the first read (reads are mapped in reverse)
                    if (all_reads[bam_number].count(al.Name)) {

                        all_reads[bam_number][al.Name].Read2_Position = al.Position;

                        all_reads[bam_number][al.Name].Read2_Aligned = al.AlignedBases;
                        all_reads[bam_number][al.Name].read2_haps[to_string(position1)] = al.AlignedBases[position1 - al.Position - 1];

     //                 store base unformation for later checking long range haplotypes
                        base1 = al.AlignedBases[position1-al.Position-1];
                        base2 = al.AlignedBases[position2-al.Position-1];
                        if (base1 == "A" or base1 == "T" or base1 == "G" or base1 == "C"){
                        bases[location1] += base1;
                        }
                        if (base2 == "A" or base2 == "T" or base2 == "G" or base2 == "C"){
                        bases[location2] += base2;
                        }

                    }
                    else {
                        all_reads[bam_number][al.Name].Read2_Position = al.Position;
                        all_reads[bam_number][al.Name].Read2_Aligned = al.AlignedBases;
                        all_reads[bam_number][al.Name].read2_haps[to_string(position1)] = al.AlignedBases[position1 - al.Position - 1];

     //                 store base unformation for later checking long range haplotypes
                        base1 = al.AlignedBases[position1-al.Position-1];
                        base2 = al.AlignedBases[position2-al.Position-1];
                        if (base1 == "A" or base1 == "T" or base1 == "G" or base1 == "C"){
                        bases[location1] += base1;
                        }
                        if (base2 == "A" or base2 == "T" or base2 == "G" or base2 == "C"){
                        bases[location2] += base2;
                        }

                    }
               }

        }

//return the result of the get_haplotypes function passing the paired_reads map just constructed
if (!reads.empty()) {
     reader.Jump(ref, position1-100);
     return get_haplotypes(reads, position1, position2, read_length, bases);
}
else {
     return {{"NN",0}};
}

}


map <string, map<string, int> > get_long_range(std::map<string, paired_reads> all_reads, map<int, string> reference_names)
{
   map <string, map<string, int> > haplotypes;
   std::string hap1;
   std::string hap2;
   std::string location;
   std::string hap;

   for (const auto name : all_reads) {

      for (const auto position1 : all_reads[name.first.c_str()].read1_haps) {

         hap1 = all_reads[name.first.c_str()].read1_haps[position1.first.c_str()];
         if (hap1 != "A" and hap1 != "G" and hap1 != "C" and hap1 != "T") { 
              continue;
         }

         for (const auto position2 : all_reads[name.first.c_str()].read2_haps) {

            hap2 = all_reads[name.first.c_str()].read2_haps[position2.first.c_str()];

            if (hap2 != "A" and hap2 != "G" and hap2 != "C" and hap2 != "T") { 
               continue;
            }


            location=reference_names[all_reads[name.first.c_str()].Ref_ID] + ":" + position1.first.c_str();
            location = location + ":" + position2.first.c_str();

            hap = hap1;
            hap += hap2;

            if (haplotypes[location].find(hap) != haplotypes[location].end()) {
               haplotypes[location][hap] += 1;
            }
            else if (!haplotypes[location].count(hap)) {
               haplotypes[location][hap] = 1;
            }
         }
      }
   }
/*   for (const auto loc : haplotypes) {
      for (const auto type : haplotypes[loc.first.c_str()]) {
         cout << loc.first.c_str() << ":" << type.first.c_str() << ":" << haplotypes[loc.first.c_str()][type.first.c_str()] << endl;
      }
   }*/
   return haplotypes;
}


int main(int argc, char* argv[])
{
    string input_bams;
    int num_bams;
    int read_length;
    string variant_file;
    string handle;
    string genfile;
    map<int, string> reference_names;

//  collect parsed arguments
    if ( argc == 6 )
    {
            variant_file = argv[1];
            input_bams  = argv[2];
            num_bams = atoi(argv[3]);
            read_length = atoi(argv[4]);
            handle = argv[5];
    }
    else {
        cerr << "Wrong number of arguments." << endl;
        return 1;
    }

    //open genotype file for outputting contrasting haplotypes
    genfile = handle+".gen";

    std::ofstream outfile (genfile, std::ofstream::out);

    //split bamfiles into array of filenames
    std::istringstream buf(input_bams);
    std::istream_iterator<std::string> beg(buf), end;

    std::vector<std::string> bams(beg, end);

    //open variant file for reading
    std::ifstream infile(variant_file);

    BamReader reader[num_bams];
    BamReader reader2[num_bams];
    BamReader reader3[num_bams];
    BamReader reader4[num_bams];

    //open bam files for overlapping Bamreader objects
    for (int x=0; x < num_bams; x++) {

        if (!reader[x].Open(bams[x])) {
            cerr << "Could not open input BAM file." << endl;
            return 1;
        }
    }

    for (int x=0; x < num_bams; x++) {

        if (!reader2[x].Open(bams[x])) {
            cerr << "Could not open input BAM file." << endl;
            return 1;
        }
    }

    for (int x=0; x < num_bams; x++) {

        if (!reader3[x].Open(bams[x])) {
            cerr << "Could not open input BAM file." << endl;
            return 1;
        }
    }

    for (int x=0; x < num_bams; x++) {

        if (!reader4[x].Open(bams[x])) {
            cerr << "Could not open input BAM file." << endl;
            return 1;
        }
    }


    //collect positions
    map<int, int> positions;
    map<int, string> chromosomes;
    map<string, int> unique_chr;
    std::string line = "";
    int position = 0;
    int chr_count = 0;

    while (std::getline(infile,line)) {
        std::string prefix("#");

        if (!line.compare(0, prefix.size(), prefix)) {
            continue;
        }

        std::stringstream line_ss(line);
        std::string column = "";
        unsigned int index = 0;
        string row[9+num_bams];

        while (std::getline(line_ss,column,'\t')) {
            if(index < 3) {
                row[index] = column;
                index++;
            }
            else{
                break;
            }
        }
       positions[position] = stoi(row[1]);

       if (unique_chr.find(row[0]) == unique_chr.end()) {
          unique_chr[row[0]] = 1;
          reference_names[chr_count] = row[0];
          chr_count++;
       }

       chromosomes[position] = row[0];
       position++;
    }

    //map that will hold all of the reads
    std::map<int, std::map<string, paired_reads> > all_reads;

    //loop through positions to collect haplotypes
    int first = 0;
    string first_chr;
    int second = 0;
    string second_chr;
    //short range haplotypes
    map <string, map<string, int> > haplotypes[num_bams];

    //map to control long range haplotypes
    map <string, string> bases[num_bams];

    //long range haplotypes
    map <string, map<string, int> > long_haplotypes[num_bams];

    string physical_location;
    int size = position;
    position = 0;
    int second_position = 0;
    int reference_ID = 0;

    int read_decision = 0;

    bool gatekeeper = 1;

    for (int chr = 0; chr < chr_count; chr++) {

    while (position < size) {


//         printf("%s\n", (chromosomes[position]+":"+to_string(positions[position])+":"+to_string(positions[position+1])).c_str());

         //If change of chromosome move to next position, clear memory for previous chromosome and go to start of loop
         if (chromosomes[position+1] != chromosomes[position]) {
              position += 1;
              break;

         }

         reference_ID = reader[0].GetReferenceID(chromosomes[position]);

         //If physical distance is less than read length iterate second_position so it is > read_length and < insert_size
         if (positions[position+1] - positions[position] < read_length) {

             physical_location = chromosomes[position]+":"+to_string(positions[position])+":"+to_string(positions[position+1]);

             switch (read_decision) {
                 case 0:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader[x],reference_ID,positions[position],positions[position+1],read_length,x,all_reads,bases[x],reference_names);
                         }
                         read_decision = 1;
                         break;
                 case 1:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader2[x],reference_ID,positions[position],positions[position+1],read_length,x,all_reads,bases[x],reference_names);
                         }
                         read_decision = 2;
                         break;
                 case 2:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader3[x],reference_ID,positions[position],positions[position+1],read_length,x,all_reads,bases[x],reference_names);
                         }
                         read_decision = 3;
                         break;
                 case 3:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader4[x],reference_ID,positions[position],positions[position+1],read_length,x,all_reads,bases[x],reference_names);
                         }
                         read_decision = 0;
                         break;


              }

              position += 1;
              continue;
         }
         position += 1;
     }

    //Filter for minimum depth
    map <string, map<string, int> > min_depth_haplotypes[num_bams];

    for (const auto pos : haplotypes[0]) {
         int retain = 1;
         for (int bam = 0; bam < num_bams; bam++) {

              for (const auto hap : haplotypes[bam][pos.first.c_str()]) {

                   if (haplotypes[bam][pos.first.c_str()][hap.first.c_str()] < 2) {
                        retain = 0;
                        break;
                   }
              }
         }
         if (retain) {
              for (int x = 0; x < num_bams; x++) {
                   min_depth_haplotypes[x][pos.first.c_str()] = haplotypes[x][pos.first.c_str()];
              }
         }
    }

    //filter for ratio of minimum haplotype < 0.25 of maximum haplotype
    map <string, map<string, int> > ratio_haplotypes[num_bams];

    for (const auto pos : min_depth_haplotypes[0]) {
         int count = 0;
         int min = 0;
         int max = 0;

         for (int bam = 0; bam < num_bams; bam++) {

              for (const auto hap : min_depth_haplotypes[bam][pos.first.c_str()]) {
                   if (count > 0) {
                       if (min > min_depth_haplotypes[bam][pos.first.c_str()][hap.first.c_str()])
                           min = min_depth_haplotypes[bam][pos.first.c_str()][hap.first.c_str()];
                       if (max < min_depth_haplotypes[bam][pos.first.c_str()][hap.first.c_str()])
                           max = min_depth_haplotypes[bam][pos.first.c_str()][hap.first.c_str()];
                   }
                   else if (count == 0) {
                       min = min_depth_haplotypes[bam][pos.first.c_str()][hap.first.c_str()];
                       max = min_depth_haplotypes[bam][pos.first.c_str()][hap.first.c_str()];
                   }
                   count++;
              }
         }

         if (double(min)/double(max) > 0.25) {
              for (int x = 0; x < num_bams; x++) {
                   ratio_haplotypes[x][pos.first.c_str()] = min_depth_haplotypes[x][pos.first.c_str()];
              }
         }
    }

    //get contrasting haplotypes
    map <string, map<string, int> > contrasted_haplotypes[num_bams];
    bool keep = 0;
    bool same_size = 1;
    for (const auto pos : ratio_haplotypes[0]) {
         same_size = 1;
         keep = 0;
         string * keepArray;
         keepArray = new string[num_bams];
         for (int bam = 1; bam < num_bams; bam++) {

              if (ratio_haplotypes[bam][pos.first.c_str()].size() != ratio_haplotypes[bam-1][pos.first.c_str()].size()) {same_size = 0; break;}

              for (const auto hap : ratio_haplotypes[bam][pos.first.c_str()]) {

                   if (!ratio_haplotypes[bam-1][pos.first.c_str()].count(hap.first.c_str())) {

                        for (const auto check: ratio_haplotypes[bam-1][pos.first.c_str()]) {

                             if (hap.first.c_str()[0] != check.first.c_str()[0] and hap.first.c_str()[1] == check.first.c_str()[1]) {
                                  keep = 1;
                                  keepArray[bam-1] = check.first.c_str();
                                  keepArray[bam] = hap.first.c_str();
                             }
                             else if (hap.first.c_str()[1] != check.first.c_str()[1] and hap.first.c_str()[0] == check.first.c_str()[0]) {
                                  keep = 1;
                                  keepArray[bam-1] = check.first.c_str();
                                  keepArray[bam] = hap.first.c_str();
                             }

                        }
                   }
              if (keep) { break; }
              }
         if (keep) { break; }
         }
         if (keep and same_size) {
              for (int x = 0; x < num_bams; x++) {
                   contrasted_haplotypes[x][pos.first.c_str()] = ratio_haplotypes[x][pos.first.c_str()];
              }
              //output genotypes of contrasting haplotypes for each bamfile
              for (int bam = 0; bam < num_bams; bam++) {

                   if (bam + 1 == num_bams) {

                        for (const auto hap : ratio_haplotypes[bam][pos.first.c_str()]) {

                             if (!ratio_haplotypes[0][pos.first.c_str()].count(hap.first.c_str()) and hap.first.c_str() != "NN" and hap.first.c_str() == keepArray[bam]) {
                                  outfile << hap.first.c_str();
                                  break;
                             }
                        }
                        break;
                   }
                   for (const auto hap : ratio_haplotypes[bam][pos.first.c_str()]) {

                        if (!ratio_haplotypes[bam+1][pos.first.c_str()].count(hap.first.c_str()) and hap.first.c_str() != "NN" and (bam + 1) < num_bams) {

                                  if (bam == 0) { outfile << pos.first.c_str(); outfile << "\t"; outfile << hap.first.c_str(); outfile << "\t"; }
                                  else if (bam > 0 and bam+1 < num_bams) { outfile << hap.first.c_str(); outfile << "\t"; }
                                  break;
                        }
                        else if (hap.first.c_str() == "NN") { outfile << "NN\t"; }
                   }
              }
              outfile << "\n";
         }
         delete[] keepArray;
    }

    string location = "";
    string hap_key = "";

    //print haplotypes with counts for each genotype
    for (const auto pos : contrasted_haplotypes[0]) {
         location = pos.first.c_str();
         printf("%s\t", location.c_str());
         for (int bam = 0; bam < num_bams; bam++) {
              for (const auto hap : contrasted_haplotypes[bam][location]) {
                   hap_key = hap.first.c_str();
                   printf("%s {%i} ", hap_key.c_str(), contrasted_haplotypes[bam][location][hap_key]);
              }
              printf("\t");
         }
         printf("\n");
    }


    //process long range haplotypes
    for (int bam = 0; bam < num_bams; bam++) {
        long_haplotypes[bam] = get_long_range(all_reads[bam], reference_names);

    }

/*    PROCESS LONG RANGE HAPLOTYPES *************************************************************
************************************************************************************************/

   //filter long range haplotypes for minimum depth
    map <string, map<string, int> > longrange_min[num_bams];

    for (const auto pos : long_haplotypes[0]) {
         int retain = 1;
         for (int bam = 0; bam < num_bams; bam++) {

              for (const auto hap : long_haplotypes[bam][pos.first.c_str()]) {

                   if (long_haplotypes[bam][pos.first.c_str()][hap.first.c_str()] < 3) {
                        retain = 0;
                        break;
                   }
              }
         }
         if (retain) {
              for (int x = 0; x < num_bams; x++) {
                   longrange_min[x][pos.first.c_str()] = long_haplotypes[x][pos.first.c_str()];
              }
         }
    }


    //filter for ratio of minimum long range haplotype < 0.25 of maximum haplotype
    map <string, map<string, int> > longrange_ratio[num_bams];

    for (const auto pos : longrange_min[0]) {
         int count = 0;
         int min = 0;
         int max = 0;

         for (int bam = 0; bam < num_bams; bam++) {

              for (const auto hap : longrange_min[bam][pos.first.c_str()]) {
                   if (count > 0) {
                       if (min > longrange_min[bam][pos.first.c_str()][hap.first.c_str()])
                           min = longrange_min[bam][pos.first.c_str()][hap.first.c_str()];
                       if (max < longrange_min[bam][pos.first.c_str()][hap.first.c_str()])
                           max = longrange_min[bam][pos.first.c_str()][hap.first.c_str()];
                   }
                   else if (count == 0) {
                       min = longrange_min[bam][pos.first.c_str()][hap.first.c_str()];
                       max = longrange_min[bam][pos.first.c_str()][hap.first.c_str()];
                   }
                   count++;
              }
         }

         if (double(min)/double(max) > 0.25) {
              for (int x = 0; x < num_bams; x++) {
                   longrange_ratio[x][pos.first.c_str()] = longrange_min[x][pos.first.c_str()];
              }
         }
    }


    //get contrasting haplotypes
    map <string, map<string, int> > contrasted_longrange[num_bams];
    keep = 0;
    same_size = 1;
    string location_index = "";

    for (const auto pos : longrange_ratio[0]) {
         std::string placeholder = "";
         unsigned int split = 0;

         location_index = "";
         same_size = 1;
         keep = 0;
         string * keepArray;
         keepArray = new string[num_bams];
         for (int bam = 1; bam < num_bams; bam++) {

              if (longrange_ratio[bam][pos.first.c_str()].size() != longrange_ratio[bam-1][pos.first.c_str()].size()) {same_size = 0; break;}
              location_index = "";

              for (const auto hap : longrange_ratio[bam][pos.first.c_str()]) {

                   if (!longrange_ratio[bam-1][pos.first.c_str()].count(hap.first.c_str())) {

                        for (const auto check: longrange_ratio[bam-1][pos.first.c_str()]) {

                             if (hap.first.c_str()[0] != check.first.c_str()[0] and hap.first.c_str()[1] == check.first.c_str()[1]) {

//                                  if (longrange_ratio[bam][pos.first.c_str()].size() > 1 and longrange_ratio[bam-1][pos.first.c_str()].size() > 1) {keep = 1; break;}

                                  std::stringstream pos_ss(pos.first.c_str());
                                  while (std::getline(pos_ss,placeholder,':')) {
                                     if(split == 0) { location_index += placeholder + ":"; split++; continue; }
                                     if(split == 1) {
                            //            location_index += ":";
                                        location_index += placeholder;
                                        split++;
                                        break;
                                     }
                                  }

                                  if (bases[bam][location_index].find(check.first.c_str()[0]) != std::string::npos or bases[bam-1][location_index].find(hap.first.c_str()[0]) != std::string::npos) {
                                     keep = 0; break;
                                  }

                                  keep = 1;
                                  keepArray[bam-1] = check.first.c_str();
                                  keepArray[bam] = hap.first.c_str();
                             }
                             else if (hap.first.c_str()[1] != check.first.c_str()[1] and hap.first.c_str()[0] == check.first.c_str()[0]) {

//                                  if (longrange_ratio[bam][pos.first.c_str()].size() > 1 and longrange_ratio[bam-1][pos.first.c_str()].size() > 1) {keep = 1; break;}

                                  std::stringstream pos_ss(pos.first.c_str());
                                  while (std::getline(pos_ss,placeholder,':')) {
                                     if(split == 0) { location_index += placeholder + ":"; split++; continue;}
                                     if(split == 2) {
                            //            location_index += ":";
                                        location_index += placeholder;
                                        break;
                                     }
                                     split++;
                                  }
                                  if (bases[bam][location_index].find(check.first.c_str()[1]) != std::string::npos or bases[bam-1][location_index].find(hap.first.c_str()[1]) != std::string::npos) {

                                     keep = 0;
                                     break;
                                  }
                                  keep = 1;
                                  keepArray[bam-1] = check.first.c_str();
                                  keepArray[bam] = hap.first.c_str();
                             }

                        }
                   }
              if (keep) { break; }
              }
         if (keep) { break; }
         }
         if (keep and same_size) {
              for (int x = 0; x < num_bams; x++) {
                   contrasted_longrange[x][pos.first.c_str()] = longrange_ratio[x][pos.first.c_str()];
              }
              //output genotypes of contrasting haplotypes for each bamfile
              for (int bam = 0; bam < num_bams; bam++) {

                   if (bam + 1 == num_bams) {

                        for (const auto hap : contrasted_longrange[bam][pos.first.c_str()]) {

                             if (!contrasted_longrange[0][pos.first.c_str()].count(hap.first.c_str()) and hap.first.c_str() != "NN" and hap.first.c_str() == keepArray[bam]) {
                                  outfile << hap.first.c_str();
                                  break;
                             }
                        }
                        break;
                   }
                   for (const auto hap : contrasted_longrange[bam][pos.first.c_str()]) {

                        if (!contrasted_longrange[bam+1][pos.first.c_str()].count(hap.first.c_str()) and hap.first.c_str() != "NN" and (bam + 1) < num_bams) {

                                  if (bam == 0) { outfile << pos.first.c_str(); outfile << "\t"; outfile << hap.first.c_str(); outfile << "\t"; }
                                  else if (bam > 0 and bam+1 < num_bams) { outfile << hap.first.c_str(); outfile << "\t"; }
                                  break;
                        }
                        else if (hap.first.c_str() == "NN") { outfile << "NN\t"; }
                   }
              }
              outfile << "\n";
         }
         delete[] keepArray;
    }

/*    for (const auto pos : contrasted_longrange[0]) {
         for (int bam = 0; bam < num_bams; bam++) {
              for (const auto hap : contrasted_longrange[bam][pos.first.c_str()]) {
                 cout << bam << ":" << pos.first.c_str() << ":" << hap.first.c_str() << ":" << contrasted_longrange[bam][pos.first.c_str()][hap.first.c_str()] << endl;
              }
         }
    }
*/

/**************END OF LONG RANGE PROCESSING****************************
**********************************************************************/

    //print haplotypes with counts for each genotype
    for (const auto pos : contrasted_longrange[0]) {
         location = pos.first.c_str();
         printf("%s\t", location.c_str());
         for (int bam = 0; bam < num_bams; bam++) {
              for (const auto hap : contrasted_longrange[bam][location]) {
                   hap_key = hap.first.c_str();
                   printf("%s {%i} ", hap_key.c_str(), contrasted_longrange[bam][location][hap_key]);
              }
              printf("\t");
         }
         printf("\n");
    }

    for (int x = 0; x < num_bams; x++) {

        all_reads[x].clear();
        bases[x].clear();
        haplotypes[x].clear();
        long_haplotypes[x].clear();
        longrange_min[x].clear();
        longrange_ratio[x].clear();
        contrasted_longrange[x].clear();
        contrasted_haplotypes[x].clear();
        min_depth_haplotypes[x].clear();
        ratio_haplotypes[x].clear();
    }

    }//loop to process a chromosome at a time

    for (int x = 0; x < num_bams; x++) {
        reader[x].Close();
        reader2[x].Close();
        reader3[x].Close();
        reader4[x].Close();
    }

    return 0;
}


