#!/usr/bin/perl
use strict;
use warnings;

my $in_list = shift;
my $out_dir = shift;
`mkdir -p $out_dir`;

my @files = ();
open(FLIST, "$in_list");
while(<FLIST>){
		chomp;
		if(!-f $_){ print STDERR "check $_\n"; exit(); }
		push(@files,$_);
}
close(FLIST);
open(FU, ">$out_dir/UMRs.bed");
open(FL, ">$out_dir/LMRs.bed");
#print FU "chr\tstart\tend\tvalue1\n";
#print FL "chr\tstart\tend\tvalue1\n";

for(my $i = 0; $i <= $#files; $i++){
		open(FTXT, "$files[$i]");
		while(<FTXT>){
				chomp;
				my @t = split(/\t/,$_);
				if($t[0] eq "chr"){ next; }
				my $start = $t[1]-1;
				if($t[3] eq "UMR"){ print FU "$t[0]\t$start\t$t[2]\t$t[6]\n"; }
				if($t[3] eq "LMR"){ print FL "$t[0]\t$start\t$t[2]\t$t[6]\n"; }
		}
		close(FTXT);
}
