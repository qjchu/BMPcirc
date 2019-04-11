# BMPcirc
Biogenesis Mechanism Predictor of circRNA

BMPcirc is a software to predict the biogenesis mechanism of circRNA.

Version: 1.0

Last Modified: 2016-6-15

Author: Qinjie Chu (qinjiechu@zju.edu.cn)

# Prerequisites
1. perl v5.12.4 or later

2. ViennaRNA Packages (already in this package, you needn't to download or install it.)

http://www.tbi.univie.ac.at/RNA/#
[Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L. ViennaRNA Package 2.0 Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26]


# Usage
unpack the zip file, and "perl ***.pl".

#<1>#Usage of rbp.pl - predict which RBP causes the biogenesis of circRNA

###Usage: perl rbp.pl [options]

 required:
 
  -a	the left intron sequence of circRNA (fasta format).
  
  -b	the right intron sequence of circRNA (fasta format).
  
  -g	the species used (human/mouse drosophila or arabidopsis). ##human and mouse use the same motif file##
  
  -o	output file name.
  
 optional:
 
  -l	the length of flanking sequences of circRNA, which will be used to detect RBP motifs [default: input left/right intron length].


#<2># Usage of pre_db.pl - prepare files that would be used in reverse complemnt mechanism or lariat mechanism prediction

###Usage: perl pre_db.pl [options]

 required:
 
  -a	the left intron sequence of circRNA (fasta format).  #sequence name should be same as that in right intron sequence file; the length of sequence name should shorter than 20 characters, and without special characters such as "|", and should be started with letters.
  
  -b	the right intron sequence of circRNA (fasta format).  #sequence name should be same as that in left intron sequence file; the length of sequence name should shorter than 20 characters, and without special characters such as "|", and should be started with letters.
  
  -o	output database (prepare files) name.


#<3>#Usage of rc.pl - predict whether reverse complementary (rc) sequence causes the biogenesis of circRNA

#####ATTENTION: pre_db.pl should be ran before running rc.pl

###Usage: perl rc.pl [options]

 required:
 
  -p	the name of the database, which was created by pre_db.pl
  
  -o	output file name.
  
 optional:
 
  -e	entropy value [default = 0.1, the shorter value the more credible the result is].
  
  -n	gaps between complementary sequence [default = 0, the shorter gaps the more credible the result is].


#<4>#Usage of lariat.pl - predict whether lariat causes the biogenesis of circRNA

#####ATTENTION: pre_db.pl should be ran before running lariat.pl

###Usage: perl lariat.pl [options]

 required:
 
  -p	the name of the database, which was created by pre)db.pl
  
  -o	output file name.
  
 optional:
 
  -e	entropy value [default = 5, the larger value the more credible the result is].
  
  -n	gaps between complementary sequence [default = 2, the larger gaps the more credible the result is].

