#!/usr/bin/perl

#use strict;
#use warnings;

use Getopt::Std;
use vars qw( $opt_a $opt_b $opt_o);

# Usage
my $usage = "

pre_db.pl - prepare files that would be used in reverse complemnt mechanism or lariat mechanism prediction
      by  Qinjie Chu (Version: 1.0  Last Modified: 2016-6-15)

Usage: perl basepair.pl [options]
 required:
  -a	the left intron sequence of circRNA (fasta format).  #sequence name should be same as that in right intron sequence file; the length of sequence name should shorter than 20 characters, and without special characters such as "|", and should be started with letters.
  -b	the right intron sequence of circRNA (fasta format).  #sequence name should be same as that in left intron sequence file; the length of sequence name should shorter than 20 characters, and without special characters such as "|", and should be started with letters.
  -o	output database (prepare files) name.

";

# command line processing.
getopts('a:b:o:');
die $usage unless $opt_a && $opt_b && $opt_o;

my ($fa1, $fa2, $out, $en_value, $gap );
$fa1 = $opt_a if $opt_a;
$fa2 = $opt_b if $opt_b;
$out = $opt_o if $opt_o;


# join two introns into one
print "Joining the two introns (file) into one sequence (file)...";

my ($name1, $len1, %ha_len, %ha_seq);
open IN1, "$fa1" or die "No left intron sequence file!";
while (<IN1>){
	chomp;
	if(/>(.*)/){
		$name1=$1;
		}
	else{
		$len1 = length $_;
		$ha_len{$name1} = $len1;
		$ha_seq{$name1} = $_;
		}
	}
close IN1;

my ($name2, $len2, $seq );
open IN2, "$fa2" or die "No right intron sequence file!";
open OUT, ">$out\_2introns.fa";
while (<IN2>){
	chomp;
	if(/>(.*)/){
		$name2=$1;
		}
	else{
		$len2=length $_;
		$seq=$_;
		if (exists $ha_len{$name2}){
			print OUT ">$name2|$ha_len{$name2}|$len2\n$ha_seq{$name2}$seq\n";
			}
		}
	}
close IN2;
close OUT;
print "done\n";


# RNAfold RNAfold_p RNAeval
print "RNAfold progressing...";

my $dirname1="$out\_RNAfold";
system(qq(mkdir $dirname1));
system(qq(RNAfold < $out\_2introns.fa > $out\_RNAfold.txt));
system(qq(mv -f *.ps $dirname1/));
print "done\n";


print "RNAeval progressing...";
system(qq(RNAeval -v -d2 < $out\_RNAfold.txt > $out\_RNAfold_RNAeval.txt));
print "done\n";


print "RNAfold_p progressing...";
my $dirname2="$out\_RNAfold_p";
system(qq(mkdir $dirname2));
system(qq(RNAfold -p < $out\_2introns.fa > $out\_RNAfold_p.txt));
system(qq(mv -f *.ps $dirname2/));
print "done\n";


# mountain 
print "mountain.pl progressing...";
opendir D, "$dirname2";
my @arrey =grep(/dp\.ps/,readdir D);
my $list;
foreach $list(@arrey){
	system(qq(perl -e 'print "$list\n"' > out1.txt));
	$list =~ s/\|/\\\|/g;
	system(qq(ViennaRNA-2.1.9/Utils/mountain.pl $dirname2/$list > out2.txt));
	system(qq(cat out1.txt out2.txt > $list.txt));
	}
system(qq(cat *_dp.ps.txt > $out\_mountain.txt));
system(qq(rm -f *_dp.ps.txt));
system(qq(rm -f out1.txt));
system(qq(rm -f out2.txt));
closedir D;
print "done\n";


##database prepare 
print "Database preparing...";

# RNAeval to basepair
my (@bp_out, %ha_bp, $bb);

open RNAeval, "$out\_RNAfold_RNAeval.txt";
while(<RNAeval>){
	chomp;
	if(/>/){
		push @bp_out, $_;#name is >3374|74|132
		%ha_bp=();
		}
	if(/loop/){
		if(/,/){
			s/,/\t/g;
			s/\(/\t/g;
			s/\)/\t/g;
			s/ //g;
			my @aa=split/\t/, $_;
			$ha_bp{$aa[1]}=$aa[2];
			$ha_bp{$aa[4]}=$aa[5];
			$ha_bp{$aa[2]}=$aa[1];
			$ha_bp{$aa[5]}=$aa[4];
			}
		}
	if(/^[AUCG]+$/){
		my @bb=split//, $_;
		foreach $bb(0..$#bb){
			my $site=$bb+1;
			push @bp_out, "$bb[$bb]\t$site\t$ha_bp{$site}";
			}
		}
	}
close RNAeval;

####print "@bp_out\n";####
####example of @bp_out####
####>3374|74|132##########
####G	1	##################
####U	2	##################
####A	3	##################
####A	4	##################
####C	5	198###############
####A	6	197###############
####C	7	196###############
####C	8	194###############
####U	9	193###############
####U	10	192#############

# mountain to 3 values
my ($type, $count, $prelen, $t, $len, %ha_value);
my ($key, %ha_all, %ha_mfe, %ha_en);

open MOU, "$out\_mountain.txt";
while (<MOU>){
	chomp;
	if(/(.+)_dp\.ps/){###example:chr10_+_DDX21_202|1132|1771_dp.ps
		$type=">$1";
		$count=0;
		$prelen=0;
		$t=0;
		$len=0;
		%ha_value=();
		}
	if(/^\s/){###example:  19   10.445
		s/\s+/\t/g;
		my @cc=split/\t/,$_;
		$count++;
		$ha_value{$count}.="$cc[2]\t";
		$prelen=++$t/3;
		if($prelen == $len){
			foreach $key(sort {$a<=>$b} keys %ha_value){
				my @dd=split/\t/, $ha_value{$key};
				$ha_all{$type}{$key}=$dd[0];
				$ha_mfe{$type}{$key}=$dd[1];
				$ha_en{$type}{$key}=$dd[2];
				}
			}
		}
	if(/^(\d+)\s+/){###example:1030   389.86
		s/\s+/\t/g;
		my @cc=split/\t/,$_;
		$count++;
		$ha_value{$count}.="$cc[1]\t";
		$prelen=++$t/3;
		if($prelen == $len){
			foreach $key(sort {$a<=>$b} keys %ha_value){
				my @dd=split/\t/, $ha_value{$key};
				$ha_all{$type}{$key}=$dd[0];
				$ha_mfe{$type}{$key}=$dd[1];
				$ha_en{$type}{$key}=$dd[2];
				}
			}
		}
	if(/&/){###example:&
		$len=$count;
		$count=0;
		}
	}
close MOU;

my $type2;
open VAL, ">$out\_3values.txt";
foreach (@bp_out){
	if(/>/){
		$type2=$_;
		print VAL "$type2\n";
		}
	else{
		my $bp=$_;
		my @ee=split/\t/,$bp;
		print VAL "$bp\t$ha_all{$type2}{$ee[1]}\t$ha_mfe{$type2}{$ee[1]}\t$ha_en{$type2}{$ee[1]}\n";
		}
	}
close VAL;
print "done\n";

#####example of 3values.txt########
#####>3374|74|132##################
#####G	1		0	0	1.3266#############
#####U	2		0.44334	0	1.2129#######
#####A	3		0.82498	0	0.10579######
#####A	4		0.83473	0	1.1108#######
#####C	5	198	1.2118	0	1.8166#####
#####A	6	197	1.9996	1	2.0379#####
#####C	7	196	2.9076	2	2.1509#####
#####C	8	194	3.6841	3	1.9727#####
#####U	9	193	4.3937	4	1.4879#####
#####U	10	192	5.1258	5	1.2302###

