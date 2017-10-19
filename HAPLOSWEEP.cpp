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
         paired_reads() : Ref_ID(0), Read1_Position(0), Read1_Aligned() {}
         paired_reads(int ref, int read1, std::string aligned1)
         : Ref_ID(ref), Read1_Position(read1), Read1_Aligned(aligned1) {}
	string Name;
	int Ref_ID;
	int Read1_Position;
	string Read1_Aligned;
};


//function to gather haplotypes from reads structure
//Takes a map of paired_reads objects and key of read name and returns a map of counts with haplotypes

map<string, int> get_haplotypes(std::map<std::string, paired_reads> reads, int position1, int position2, int read_length)
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
//        cout << x.first.c_str() << endl;
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

map<string, int> collect_haplotype(BamReader& reader, int ref, int position1, int position2, int read_length)
{

        std::map<std::string, paired_reads> reads;
        string mate1 = "";

//        reader.Rewind();

        //set the region to the first position
        reader.SetRegion(ref,position1,ref,position1);

        BamAlignment al;


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
	       if (al.IsMapped() and abs(position1-al.Position) < read_length) {
//		    add to the paired_reads struct
                    reads[al.Name] = paired_reads(al.RefID,al.Position,al.AlignedBases);

               }
//             if the read is actually the second mate
               else if (al.IsSecondMate()) {
                    //if it is not in the map add as the first read (reads are mapped in reverse)
                    reads[al.Name+"_2"] = paired_reads(al.RefID,al.Position,al.AlignedBases);
               }

        }

//return the result of the get_haplotypes function passing the paired_reads map just constructed
if (!reads.empty()) {
     reader.Jump(ref, position1-100);
     return get_haplotypes(reads, position1, position2, read_length);
}
else {
     return {{"NN",0}};
}

}


int main(int argc, char* argv[])
{
    string input_bams;
    int num_bams;
    int read_length;
    string variant_file;
    string handle;
    string genfile;

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
            cerr << "Could not open input BAM file1: " << bams[x] << endl;
            return 1;
        }
    }

    for (int x=0; x < num_bams; x++) {

        if (!reader2[x].Open(bams[x])) {
            cerr << "Could not open input BAM file2: " << bams[x] << endl;
            return 1;
        }
    }

    for (int x=0; x < num_bams; x++) {

        if (!reader3[x].Open(bams[x])) {
            cerr << "Could not open input BAM file3:" << bams[x] << endl;
            return 1;
        }
    }

    for (int x=0; x < num_bams; x++) {

        if (!reader4[x].Open(bams[x])) {
            cerr << "Could not open input BAM file4:" << bams[x] << endl;
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
          //reference_names[chr_count] = row[0];
          chr_count++;
       }

       chromosomes[position] = row[0];
       position++;
    }

    //loop through positions to collect haplotypes
    int first = 0;
    string first_chr;
    int second = 0;
    string second_chr;

    map <string, map<string, int> > haplotypes[num_bams];
    string physical_location;
    int size = position;
    position = 0;
    int second_position = 0;
    int reference_ID = 0;

    int read_decision = 0;

    for (int chr = 0; chr < chr_count; chr++) {

    while (position < size) {


//         printf("%s\n", (chromosomes[position]+":"+to_string(positions[position])+":"+to_string(positions[position+1])).c_str());

         //If change of chromosome move to next position and go to start of loop
         if (chromosomes[position+1] != chromosomes[position]) {
              position += 1;
              continue;
         }
         reference_ID = reader[0].GetReferenceID(chromosomes[position]);

         //If physical distance is less than read length iterate second_position so it is > read_length and < insert_size
         if (positions[position+1] - positions[position] < read_length) {

             physical_location = chromosomes[position]+":"+to_string(positions[position])+":"+to_string(positions[position+1]);

             switch (read_decision) {
                 case 0:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader[x],reference_ID,positions[position],positions[position+1],read_length);
                         }
                         read_decision = 1;
                         break;
                 case 1:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader2[x],reference_ID,positions[position],positions[position+1],read_length);
                         }
                         read_decision = 2;
                         break;
                 case 2:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader3[x],reference_ID,positions[position],positions[position+1],read_length);
                         }
                         read_decision = 3;
                         break;
                 case 3:
                     for (int x = 0; x < num_bams; x++) {
                         haplotypes[x][physical_location] = collect_haplotype(reader4[x],reference_ID,positions[position],positions[position+1],read_length);
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
    //loop through positions
    for (const auto pos : ratio_haplotypes[0]) {
         same_size = 1;
         keep = 0;
         string keepArray1;
         string keepArray2;
         //loop through bams
         for (int bam = 1; bam < num_bams; bam++) {
              //check if there is a difference in number of haplotypes at these positions
              if (ratio_haplotypes[bam][pos.first.c_str()].size() != ratio_haplotypes[bam-1][pos.first.c_str()].size())
              {
                   same_size = 0;
                   break;
              }
              //loop through the haplotypes at this bam and these posiitions
              for (const auto hap : ratio_haplotypes[bam][pos.first.c_str()]) {
                   //check if the haplotype is in bam-1 - because if it is then move to the next hap
                   if (!ratio_haplotypes[bam-1][pos.first.c_str()].count(hap.first.c_str())) {
                        //loop through the haps in bam-1
                        for (const auto check: ratio_haplotypes[bam-1][pos.first.c_str()]) {
                             //if the hap in quextion first positions not equal but second positions are equal
                             if (hap.first.c_str()[0] != check.first.c_str()[0] and hap.first.c_str()[1] == check.first.c_str()[1]) {
                                  //keep this haplotype and enter it into the keepArray pointer
                                  keep = 1;
                                  keepArray1 = check.first.c_str();
                                  keepArray2 = hap.first.c_str();
                             }
                             else if (hap.first.c_str()[1] != check.first.c_str()[1] and hap.first.c_str()[0] == check.first.c_str()[0]) {
                                  keep = 1;
                                  keepArray1 = check.first.c_str();
                                  keepArray2 = hap.first.c_str();
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
              outfile << pos.first.c_str(); outfile << "\t";
              int found = 0;
              for (int bam = 0; bam < num_bams; bam++) {
//                   outfile << bam << endl;

                   if (bam + 1 == num_bams) {
                        found = 0;
                        for (const auto hap : contrasted_haplotypes[bam][pos.first.c_str()]) {

                             if (hap.first.c_str() != "NN" and (hap.first.c_str() == keepArray1 or hap.first.c_str() == keepArray2)) {
                                  outfile << hap.first.c_str();
                                  found = 1;
                                  break;
                             }
                             else if (hap.first.c_str() == "NN") { outfile << "NN\t"; }
                        }
                        if (!found) { outfile << "NN\t"; }
                        break;
                   }
                   found = 0;
                   for (const auto hap : contrasted_haplotypes[bam][pos.first.c_str()]) {

                        if (hap.first.c_str() != "NN" and (bam + 1) < num_bams and (hap.first.c_str() == keepArray1 or hap.first.c_str() == keepArray2)) {

                                 outfile << hap.first.c_str();
                                 outfile << "\t";
                                 found = 1;
                                 break;
                        }
                        else if (hap.first.c_str() == "NN") { outfile << "NN\t"; }
                   }
                   if (!found) { outfile << "NN\t"; }
              }
              outfile << "\n";
         }
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

    for (int x = 0; x < num_bams; x++) {

        haplotypes[x].clear();
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

