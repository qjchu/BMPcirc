# BMPcirc
Biogenesis Mechanism Predictor of circRNA

BMPcirc is a software to predict the biogenesis mechanism of circRNA.

Version: 1.0
Last Modified: 2016-6-15
Author: Qinjie Chu (qinjiechu@zju.edu.cn)

#Prerequisites
1. perl v5.12.4 or later

2. ViennaRNA Packages (already in this package, you needn't to download or install it.)
http://www.tbi.univie.ac.at/RNA/#
[Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L. ViennaRNA Package 2.0 Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26]


#Usage: unpack the zip file, and "perl ***.pl".

#<1>#Usage of rbp.pl - predict which RBP causes the biogenesis of circRNA

###Usage: perl rbp.pl [options]

 required:
  -a	the left intron sequence of circRNA (fasta format).
  -b	the right intron sequence of circRNA (fasta format).
  -g	the species used (human/mouse drosophila or arabidopsis). ##human and mouse use the same motif file##
  -o	output file name.
 optional:
  -l	the length of flanking sequences of circRNA, which will be used to detect RBP motifs [default: input left/right intron length].


###Outputs of rbp.pl
1. **_rbpout_1.txt	===>	RBP binding site of the left intron sequence of circRNA
   **_rbpout_2.txt	===>	RBP binding site of the right intron sequence of circRNA
  |||==========================================|||
  |||format of **_rbpout_1.txt/**_rbpout_2.txt |||
  |||------------------------------------------|||
  |||>sequence_name                            |||
  |||RBP_name	binding sequence	binding site   |||
  |||==========================================|||
  
PS:
  |||==========================================|||
  |||format of RBP_name                        |||
  |||------------------------------------------|||
  |||(protein name)_(species_name)_(seq_number)|||
  |||==========================================|||

2. **_final_out.txt	===>	common RBPs that bind on both right and left introns of circRNA
  |||==========================================|||
  |||format of **_final_out.txt                |||
  |||------------------------------------------|||
  |||sequence_name	RBP_name                   |||
  |||==========================================|||


#<2># Usage of pre_db.pl - prepare files that would be used in reverse complemnt mechanism or lariat mechanism prediction

###Usage: perl basepair.pl [options]
 required:
  -a	the left intron sequence of circRNA (fasta format).  #sequence name should be same as that in right intron sequence file; the length of sequence name should shorter than 20 characters, and without special characters such as "|", and should be started with letters.
  -b	the right intron sequence of circRNA (fasta format).  #sequence name should be same as that in left intron sequence file; the length of sequence name should shorter than 20 characters, and without special characters such as "|", and should be started with letters.
  -o	output database (prepare files) name.

###Outputs of pre_db.pl
1. **_2introns.fa	===>	fasta format file of left and right intron sequence

  |||============================================================|||
  |||format of **_2introns.fa                                    |||
  |||------------------------------------------------------------|||
  |||>sequece_name|left_intron_length|right_intron_length        |||
  |||left_sequence and right sequence                            |||
  |||============================================================|||

For example
left intron sequence file:  >abc
                            ATCGATCGATCGATCGATCG
right intron sequence file: >abc
                            AAAATTTTCCCCGGGG
**_2introns.fa file:        >abc|20|16
                            ATCGATCGATCGATCGATCGAAAATTTTCCCCGGGG

2. **_RNAfold.txt	===>	the minimum free energy (MFE) structure created by command "RNAfold < **_2introns.fa > **_RNAfold.txt"

3. **_RNAfold_RNAeval.txt	===>	thermodynamic details on MFE structure created by command "RNAeval -v -d2 < **_RNAfold.txt > **_RNAfold_RNAeval.txt"

4. **_RNAfold_p.txt	===>	the minimum free energy (MFE) and thermodynamic ensemble structure created by command "RNAfold -p < **_2introns.fa > **_RNAfold_p.txt"

5. **_mountain.txt	===>	the three values (the MFE structure value, thermodynamic ensemble value of RNA structures, and the positional entropy for each position) created by mountain.pl (ViennaRNA package)
output files (2,3,4,5) are all create by ViennaRNA Packages, the detail format of which could be found on http://www.tbi.univie.ac.at/RNA/#

6. **_3values.txt	===>	txt format file that contains the minimum free energy (MFE) structure value, thermodynamic ensemble value of RNA structures, and the positional entropy for each position.
  |||=================================================================================================================|||
  |||format of **_3values.txt                                                                                         |||
  |||-----------------------------------------------------------------------------------------------------------------|||
  |||>sequece_name|left_intron_length|right_intron_length                                                             |||
  |||(six columns: the details of each column are described as followings)                                            |||
  |||# first column   ## the base of the sequence                                                                     |||
  |||# second column  ## the order number of the corresponding base                                                   |||
  |||# third column   ## order number of the base that is paired to the base that of the second column's number       |||
  |||# fourth column  ## the minimum free energy (MFE) structure value of the base that of the second column's number |||
  |||# fifth column   ## thermodynamic ensemble value of RNA structures of the base that of the second column's number|||
  |||# sixth column   ## the positional entropy for each position of the base that of the second column's number      |||
  |||=================================================================================================================|||

7. folder: **_RNAfold	===>	created by command "RNAfold < **_2introns.fa > **_RNAfold.txt"

8. folder: **_RNAfold_p	===>	created by command "RNAfold -p < **_2introns.fa > **_RNAfold.txt"


#<3>#Usage of rc.pl - predict whether reverse complementary (rc) sequence causes the biogenesis of circRNA
#####ATTENTION: pre_db.pl should be ran before running rc.pl

###Usage: perl rc.pl [options]
 required:
  -p	the name of the database, which was created by pre_db.pl
  -o	output file name.
 optional:
  -e	entropy value [default = 0.1, the shorter value the more credible the result is].
  -n	gaps between complementary sequence [default = 0, the shorter gaps the more credible the result is].

###Outputs of rc.pl
1. **_final.txt	===>	txt format file that lists the base pairs (number of the longest base pairs)
  |||======================================================================================|||
  |||format of **_final.txt                                                                |||
  |||--------------------------------------------------------------------------------------|||
  |||sequece_name|left_intron_length|right_intron_length	number of the longest base pairs |||
  |||======================================================================================|||


#<4>#Usage of lariat.pl - predict whether lariat causes the biogenesis of circRNA
#####ATTENTION: pre_db.pl should be ran before running lariat.pl

###Usage: perl lariat.pl [options]
 required:
  -p	the name of the database, which was created by pre)db.pl
  -o	output file name.
 optional:
  -e	entropy value [default = 5, the larger value the more credible the result is].
  -n	gaps between complementary sequence [default = 2, the larger gaps the more credible the result is].

###Outputs of lariat.pl
1. **_final.txt	===>	txt format file that lists the base pairs (the sum of all base pairs)
  |||======================================================================================|||
  |||format of **_final.txt                                                                |||
  |||--------------------------------------------------------------------------------------|||
  |||sequece_name|left_intron_length|right_intron_length	the sum of all base pairs        |||
  |||======================================================================================|||

