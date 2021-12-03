#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $in_gene_bed;
my $distance_tss = 1500;
my $bin_size = 500;
my $step_size = 100;

my $options = GetOptions(
				"input=s" => \$in_gene_bed,
				"distance=i" => \$distance_tss,
				"bin_size=i" => \$bin_size,
				"step_size=i" => \$step_size,
);
#my $max_i = (int(2000-$bin_size)/$step_size)+1;
open(FG, "$in_gene_bed");
while(<FG>){
		chomp;
		my @t = split(/\t+/,$_);
#		my ($ann, $gene_id) = split(/:/,$t[3]);
		my ($i1, $id, $i2, $transcript_id, $i3, $gene_name) = split(/\s/,$t[-1]);
		$id = substr($id,1,-2);
		$transcript_id = substr($transcript_id,1,-2);
		my $gene_id = "$id;$transcript_id"; ## fixed at Apr 7, 2021
#		my $gene_id = $t[3]; # fixed at Apr 5, 2021
		if($t[5] eq "+"){
				my $tss = $t[1];
				my $end_of_bins = $tss-$distance_tss;
				if($end_of_bins < 0){ $end_of_bins = 0; }
				my $iter = 0;
				my $string = "";
				while(1){ #### Upstream
						my $bin_end = $tss-($iter*$step_size);
						my $bin_start = $bin_end-$bin_size;
						if($bin_start < $end_of_bins){ last; }
#						print "$t[0]\t$bin_start\t$bin_end\tU$iter:$gene_id\t$t[4]\t$t[5]\n";
						$string = "$t[0]\t$bin_start\t$bin_end\tU$iter:$gene_id\t$t[4]\t$t[5]\n".$string;
						$iter++;
				}
				print "$string";
				$end_of_bins = $tss+$distance_tss;
				$iter = 0;
				while(1){ #### Downstream
						my $bin_start = $tss+($iter*$step_size);
						my $bin_end = $bin_start+$bin_size;
						if($bin_end > $tss+$distance_tss){ last; }
						print "$t[0]\t$bin_start\t$bin_end\tD$iter:$gene_id\t$t[4]\t$t[5]\n";
						$iter++;
				}
		}
		else{
				my $tss = $t[2];
				my $end_of_bins = $tss+$distance_tss;
				my $iter = 0;
				my $string = "";
				while(1){ #### Upstream
						my $bin_start = $tss+($iter*$step_size);
						my $bin_end = $bin_start+$bin_size;
						if($bin_end > $tss+$distance_tss){ last; }
#						print "$t[0]\t$bin_start\t$bin_end\tU$iter:$gene_id\t$t[4]\t$t[5]\n";
						$string = "$t[0]\t$bin_start\t$bin_end\tU$iter:$gene_id\t$t[4]\t$t[5]\n".$string;
						$iter++;
				}
				print "$string";
				$end_of_bins = $tss-$distance_tss;
				if($end_of_bins < 0){ $end_of_bins = 0; }
				$iter = 0;
				while(1){ #### Downstream
						my $bin_end = $tss-($iter*$step_size);
						my $bin_start = $bin_end-$bin_size;
						if($bin_start < $end_of_bins){ last; }
						print "$t[0]\t$bin_start\t$bin_end\tD$iter:$gene_id\t$t[4]\t$t[5]\n";
						$iter++;
				}
		}
}
close(FG);
