#!/usr/bin/perl

#use strict;
use warnings;

use Getopt::Std;
use vars qw( $opt_a $opt_b $opt_g $opt_o $opt_l);

# Usage
my $usage = "

rbp.pl - predict which RBP causes the biogenesis of circRNA
      by  Qinjie Chu (Version: 1.0  Last Modified: 2016-6-9)

Usage: perl rbp.pl [options]
 required:
  -a	the left intron sequence of circRNA (fasta format).#sequence name should be same as that in right intron sequence file
  -b	the right intron sequence of circRNA (fasta format).#sequence name should be same as that in left intron sequence file
  -g	the species used (human/mouse drosophila or arabidopsis). ##human and mouse use the same motif file##
  -o	output file name.
 optional:
  -l	the length of flanking sequences of circRNA, which will be used to detect RBP motifs [default: input left/right intron length].

";

# command line processing.
getopts('a:b:g:o:l:');
die $usage unless $opt_a && $opt_b && $opt_g && $opt_o;

my ($fa1, $fa2, $motif, $out);
$fa1 = $opt_a if $opt_a;
$fa2 = $opt_b if $opt_b;
$motif = $opt_g if $opt_g;
$out = $opt_o if $opt_o;

my %ha_motif;
my @rbp_result1;
my @rbp_result2;

##### ##### ##### ##### #####


##### ##### ##### ##### #####
#### motif ha

open INS, $motif."\_motif.fa" or die "No $motif motif file!\n";
my $name;
while(<INS>){
	chomp;
	if(/>(.+)/){
		$name=$1;
		}
	else{
		$ha_motif{$name}=$_;
		}
	}
close INS;

#### left intron
my ($title1, @tit1, $seq1, $len1, $key, $start1, $end1);
my ($region1, $seq1_length);

open IN1, "$fa1" or die "No left intron file!\n";
while (<IN1>){
	chomp;
	if(/>/){
		$title1=$_;
		push @tit1,$title1;
		}
	else{
		
		$seq1_length=length $_;
		$region1 = $opt_l ? $opt_l : $seq1_length;
		$seq1=substr($_, -$region1);
		foreach $key(sort keys %ha_motif){
			$len1=length $ha_motif{$key};##the length of motif sequence
			while($seq1 =~ m/$ha_motif{$key}/g){
				$start1=pos($seq1)-$len1+1;
				$end1=pos($seq1);
				push @{$title1}, "$key\t$ha_motif{$key}\t$start1-$end1";
				pos($seq1)-=$len1-1;##change the start position of mapping
				}
			
			}
		}
	}
close IN1;

my ($i, $k);
foreach $i(@tit1){
	if(@{$i}){
		push @rbp_result1, $i;
		#print OUT "$i\n";
		foreach $k(@{$i}){
			push @rbp_result1, $k;
			#print OUT "$k";
			}
		}
	}

####print out the rbp result of left intron
open OUT1, ">$out\_rbpout_1.txt";
foreach (@rbp_result1){print OUT1 "$_\n";}
close OUT1;

#### right intron
my ($title2, @tit2, $seq2, $len2, $key2, $start2, $end2);
my ($region2, $seq2_length);

open IN2, "$fa2" or die "No right intron file!\n";
while (<IN2>){
	chomp;
	if(/>/){
		$title2=$_;
		push @tit2,$title2;
		@{$title2}=();##initialization is important
		}
	else{
		$seq2_length=length $_;
		$region2 = $opt_l ? $opt_l : $seq2_length;
		$seq2=substr($_, 0, $region2);

		foreach $key2(sort keys %ha_motif){
			$len2=length $ha_motif{$key2};##the length of motif sequence
			while($seq2 =~ m/$ha_motif{$key2}/g){
				$start2=pos($seq2)-$len2+1;
				$end2=pos($seq2);
				push @{$title2}, "$key2\t$ha_motif{$key2}\t$start2-$end2";
				pos($seq2)-=$len2-1;##change the start position of mapping
				}
			
			}
		}
	}
close IN2;

my ($i2, $k2);
foreach $i2(@tit2){
	if(@{$i2}){
		push @rbp_result2, $i2;
		#print OUT "$i\n";
		foreach $k2(@{$i2}){
			push @rbp_result2, $k2;
			#print OUT "$k";
			}
		}
	}

############################################################
#############example of rbp_result##########################
#############>At_ciR4414_ID=AT1G35220.1_+###################
#############BRUNOL6_ugugdkg_human_2	TGTGATG	232-238#######
#############HNRNPCL1_huuuuuk_human_1	ATTTTTG	221-227#######
#############HNRNPCL1_huuuuuk_human_1	ATTTTTG	267-273#######
#############HNRNPC_huuuuuk_human_1	ATTTTTG	221-227#########
#############HNRNPC_huuuuuk_human_1	ATTTTTG	267-273#########
#############IGF2BP2_vmahwca_human_18	CAATTCA	48-54#########
############################################################

####print out the rbp result of right intron
open OUT2, ">$out\_rbpout_2.txt";
foreach (@rbp_result2){print OUT2 "$_\n";}
close OUT2;

####common motif between two introns

my ($line, $head, %ha_a, @aa, @bb );
foreach $line(0..$#rbp_result1){
	if($rbp_result1[$line] =~ />(.+)/){
		$head=$1;
		}
	if($rbp_result1[$line] !~ />(.+)/){
		@aa=split/\t/, $rbp_result1[$line];
		@bb=split/_/, $aa[0];
		$ha_a{"$head|$bb[0]"}=1;
		}
	}

my ($line2, $head2, @cc, @dd, @common);
foreach $line2(0..$#rbp_result2){
	if($rbp_result2[$line2] =~ />(.+)/){
		$head2=$1;
		}
	if($rbp_result2[$line2] !~ />(.+)/){
		@cc=split/\t/, $rbp_result2[$line2];
		@dd=split/_/, $cc[0];
		if(exists $ha_a{"$head2|$dd[0]"}){
			push @common, "$head2\t$dd[0]";
			#print OUT "$head2\t$dd[0]\n";
			}
		}
	}


open OUT, ">$out\_final_out.txt";
my (@unique_common, @unique_name_a, %count, %count2, @ee, @name_a);
@unique_common=grep{++$count{$_}<2;}@common;

foreach (@unique_common){
	@ee=split/\t/,$_;
	push @name_a, $ee[0];
	push @{"$ee[0]"}, $ee[1];
	}

@unique_name_a=grep{++$count2{$_}<2;}@name_a;

my ($xx, $yy);
foreach $xx(@unique_name_a){
	print OUT "$xx\t";
	foreach $yy(@{$xx}){
		print OUT "$yy\t";
		}
	print OUT "\n";
	}
close OUT;


###end