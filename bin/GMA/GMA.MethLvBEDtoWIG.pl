#!/usr/bin/perl
use strict;
use warnings;

my $in_methyl_call = shift;

my $pre_chr = "";
open(FMC, "$in_methyl_call");
#print "chr\tstart\tend\tvalue\n";
while(<FMC>){
		chomp;
		if($_ !~ /^\d/ && $_ ne "X" && $_ ne "Y"){ next; }
		my @t = split(/\t/,$_);
		my @count = split(/\W/,$t[3]);
		my $lv = sprintf("%.3f",($count[3]/($count[1]+$count[3]))*100);
		my $cur_chr = "chr$t[0]";
		if($pre_chr ne $cur_chr){
				print "variableStep chrom=$cur_chr\n";
				$pre_chr = $cur_chr;
		}
		print "$t[2]\t$lv\n";
#		print "chr$t[0]\t$t[1]\t$t[2]\t.\t$lv\t.\n";
}
close(FMC);
