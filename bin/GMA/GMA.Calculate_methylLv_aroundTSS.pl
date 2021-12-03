#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw(:all);

my $in_intersect = shift;
my $all_gene_bed = shift;
#my $down_deg = shift;
#my $up_deg = shift;
my $distance_tss = shift;
my $step_size = shift;
my $bin_size = shift;

#### calculating the num of bins
my $tss = 0;
my $iter = 0;
while(1){ #### Downstream
		my $bin_start = $tss+($iter*$step_size);
		my $bin_end = $bin_start+$bin_size;
		if($bin_end > $tss+$distance_tss){ last; }
		$iter++;
}
my $max_i = $iter-1;
print STDERR "The name of bins : 0 ~ $max_i\n";
my %hs_d_deg = ();

my %hs_u_count = ();
my %hs_m_count = ();
my %hs_cpg_count = ();

open(FINT, "$in_intersect");
while(<FINT>){
		chomp;
		my @t = split(/\t/,$_);
		my @c = split(/\W/,$t[9]);
		my ($type, $gene_name) = split(/:/,$t[3]);
		if(!exists($hs_cpg_count{$gene_name}{$type})){ $hs_cpg_count{$gene_name}{$type} = 1; }
		else{ $hs_cpg_count{$gene_name}{$type} += 1; }

		if(!exists($hs_m_count{$gene_name}{$type})){
				$hs_m_count{$gene_name}{$type} = $c[3];
				$hs_u_count{$gene_name}{$type} = $c[1];
		}
		else{
				$hs_m_count{$gene_name}{$type} += $c[3];
				$hs_u_count{$gene_name}{$type} += $c[1];
		}
}
close(FINT);

my %hs_methyl_avg = ();
my %hs_bin_count = ();
open(FAG, "$all_gene_bed");
while(<FAG>){
		chomp;
		my @t = split(/\t/,$_);
		my ($i1, $id, $i2, $transcript_id, $i3, $gene_name) = split(/\s/,$t[-1]);
		$id = substr($id,1,-2);
		$transcript_id = substr($transcript_id,1,-2);
		$gene_name = "$id;$transcript_id"; ## fixed at Apr 7, 2021
#		my $gene_name = $t[3];
#		print STDERR "$gene_name\n";
		print "$gene_name";
		for(my $i = $max_i; $i >= 0; $i--){
				my $type = "U$i";
				if(!exists($hs_m_count{$gene_name}{$type})){
						print "\tNA";
				}
				else{
								my $total_count = $hs_m_count{$gene_name}{$type} + $hs_u_count{$gene_name}{$type};
								my $methyl_level = sprintf("%.5f",$hs_m_count{$gene_name}{$type}/$total_count);
								$hs_methyl_avg{$type}{$gene_name} = $methyl_level;
								print "\t$methyl_level";
				}
		}
		for(my $i = 0; $i <= $max_i; $i++){
				my $type = "D$i";
				if(!exists($hs_m_count{$gene_name}{$type})){
						print "\tNA";
				}
				else{
								my $total_count = $hs_m_count{$gene_name}{$type} + $hs_u_count{$gene_name}{$type};
								my $methyl_level = sprintf("%.5f",$hs_m_count{$gene_name}{$type}/$total_count);
								$hs_methyl_avg{$type}{$gene_name} = $methyl_level;
								print "\t$methyl_level";
				}
		}
		print "\n";
}
close(FAG);
