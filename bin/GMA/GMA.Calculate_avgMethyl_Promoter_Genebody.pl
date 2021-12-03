#!/usr/bin/perl
use strict;
use warnings;

my $in_intersect = shift;
my $all_gene_bed = shift;

my %hs_u_count = ();
my %hs_m_count = ();
my %hs_cpg_count = ();
my %hs_methyl = ();
open(FINT, "$in_intersect");
while(<FINT>){
		chomp;
		my @t = split(/\t/,$_);
		my @c = split(/\W/,$t[9]);
		my $meth_lv = $c[3]/($c[3]+$c[1]);
		my ($type, $gene_name) = split(/:/,$t[3]);
		if(!exists($hs_cpg_count{$gene_name}{$type})){
				$hs_cpg_count{$gene_name}{$type} = 1;
		}
		else{
				$hs_cpg_count{$gene_name}{$type} += 1;
		}

		if(!exists($hs_m_count{$gene_name}{$type})){
				$hs_m_count{$gene_name}{$type} = $c[3];
				$hs_u_count{$gene_name}{$type} = $c[1];
				$hs_methyl{$gene_name}{$type} = $meth_lv;
		}
		else{
				$hs_m_count{$gene_name}{$type} += $c[3];
				$hs_u_count{$gene_name}{$type} += $c[1];
				$hs_methyl{$gene_name}{$type} += $meth_lv;
		}
}
close(FINT);

print "Gene_id\tPromoter(mC/(mC+uC))\tPromoter(methLv/totalcount)\tGene(mC/(mC+uC))\tGene(methLv/totalcount)\n";
open(FAG, "$all_gene_bed");
while(<FAG>){
		chomp;
		my @t = split(/\t/,$_);
		my ($i1, $id, $i2, $transcript_id, $i3, $gene_name) = split(/\s/,$t[-1]);
		$id = substr($id,1,-2);
		$transcript_id = substr($transcript_id,1,-2);
		$gene_name = "$id;$transcript_id"; ## fixed at Apr 7, 2021
		print "$gene_name";
		if(!exists($hs_m_count{$gene_name}{upstream1K})){
				print "\tNA\tNA";
		}
		else{
				my $total_count = $hs_m_count{$gene_name}{upstream1K} + $hs_u_count{$gene_name}{upstream1K};
				my $methyl_level = $hs_m_count{$gene_name}{upstream1K}/$total_count;
				my $av_meth = $hs_methyl{$gene_name}{upstream1K}/$hs_cpg_count{$gene_name}{upstream1K};
				print "\t$methyl_level\t$av_meth";
		}
		if(!exists($hs_m_count{$gene_name}{gene})){
				print "\tNA\tNA\n";
		}
		else{
				my $total_count = $hs_m_count{$gene_name}{gene} + $hs_u_count{$gene_name}{gene};
				my $methyl_level = $hs_m_count{$gene_name}{gene}/$total_count;
				my $av_meth = $hs_methyl{$gene_name}{gene}/$hs_cpg_count{$gene_name}{gene};
				print "\t$methyl_level\t$av_meth\n";
		}
}
close(FAG);
