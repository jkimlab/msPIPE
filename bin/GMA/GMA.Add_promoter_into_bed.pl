#!/usr/bin/perl
use strict;
use warnings;

my $in_gene_bed = shift;
open(FG, "$in_gene_bed");
while(<FG>){
		chomp;
		my @t = split(/\t+/,$_);
		my ($i1, $id, $i2, $transcript_id, $i3, $gene_name) = split(/\s/,$t[-1]);
		$id = substr($id,1,-2);
		$transcript_id = substr($transcript_id,1,-2);
		my $gene_id = "$id;$transcript_id"; ## fixed at Apr 7, 2021
#		my $gene_id = $t[3]; ## fixed at Mar 7, 2021
		if($t[5] eq "+"){
				my $prom_end = $t[1];
				if($prom_end == 0){ next; }
				my $prom_start = $prom_end-1000;
				if($prom_start < 0){ $prom_start = 0; }
				print "$t[0]\t$prom_start\t$prom_end\tupstream1K:$gene_id\t$t[4]\t$t[5]\n";
				print "$t[0]\t$t[1]\t$t[2]\tgene:$gene_id\t$t[4]\t$t[5]\n"; ## fixed at Mar 7, 2021
		}
		else{
				my $prom_start = $t[2];
				my $prom_end = $prom_start+1000;
				print "$t[0]\t$t[1]\t$t[2]\tgene:$gene_id\t$t[4]\t$t[5]\n"; ## fixed at Mar 7, 2021
				print "$t[0]\t$prom_start\t$prom_end\tupstream1K:$gene_id\t$t[4]\t$t[5]\n";
		}
}
close(FG);
