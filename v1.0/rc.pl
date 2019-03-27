#!/usr/bin/perl

#use strict;
#use warnings;

use Getopt::Std;
use vars qw( $opt_p $opt_o $opt_e $opt_n);

# Usage
my $usage = "

rc.pl - predict whether reverse complementary sequence causes the biogenesis of circRNA
      by  Qinjie Chu (Version: 1.0  Last Modified: 2016-6-15)

Usage: perl rc.pl [options]
 required:
  -p	the name of the database, which was created by pre.pl
  -o	output file name.
 optional:
  -e	entropy value [default = 0.1, the shorter value the more credible the result is].
  -n	gaps between complementary sequence [default = 0, the shorter gaps the more credible the result is].

";

# command line processing.
getopts('p:o:e:n:');
die $usage unless $opt_p && $opt_o;

my ($db, $out, $en_value, $gap );
$db = $opt_p if $opt_p;
$out = $opt_o if $opt_o;
$en_value = $opt_e ? $opt_e : 0.1;
$gap = $opt_n ? $opt_n : 0;

## important step
my ($start, $end, $name_a, @name_b);

open INVAL, "$db\_3values.txt";
while(<INVAL>){
	chomp;
	if(/>(.+)/){
		$name_a=$1;
		push @name_b, $name_a;
		my @ff=split /\|/, $name_a;
		$start = $ff[1];
		$end = $ff[1]+$ff[2];
		}
	if(/^[AUCG]+/){
		my @gg=split /\t/, $_;
		if($gg[2] =~ /\d+/ && ($gg[5] <= $en_value) && (($gg[1]<=$start && ($gg[2]>$start)) ||($gg[2]<=$start &&($gg[1]>$start)))){
			push @{$name_a}, $gg[1];
			}
		}
	}
close INVAL;

open F, ">$out\_final.txt";
my ($i, $j, $diff, @pair);
my $k=1;
foreach $i(@name_b){
	@pair=();
	my $nums=$#{$i}-1;
	foreach $j(0..$nums){
		$diff=$gap+1;
		if(${$i}[$j+1] - ${$i}[$j] <= $diff){
			$k++;
			}
		if(${$i}[$j+1] - ${$i}[$j] > $diff){
			push @pair, $k;
			$k=1;
			push @pair, $k;
			}
		if($j == $nums && (${$i}[$j+1] - ${$i}[$j] <= $diff)){
			push @pair, $k;
			}
		}
	
	use List::Util qw(max);
	my $max = max(@pair);
	print F "$i\t$max\n";
	}
close F;
