#!/usr/bin/perl
use strict;
use warnings;

my $in_union_txt = shift;
open(FU, "$in_union_txt");
while(<FU>){
		chomp;
		my ($chr, $pos, $t, $c) = split(/\t/,$_);
		my $start = $pos-1;
		my $u = $t-$c;
#		print "$chr\t$start\t$pos\tU:$t;M:$c\t.\t.\n";
		print "$chr\t$start\t$pos\tU:$u;M:$c\t.\t.\n";
}
close(FU);
