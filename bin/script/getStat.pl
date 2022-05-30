#!/usr/bin/perl
use strict;
use warnings;
use Sort::Key::Natural qw(natsort);

my $in_intersection = $ARGV[0];
my $out_dir = $ARGV[1];

`mkdir -p $out_dir`;

my %hs_count = ();
my %hs_length = ();
my %hs_region = ();
my %hs_gene_direction = ();
#my %hs_dist_from_tss = ();
#my %hs_dist_between_dmbs = ();
my %hs_methyl_level = ();
my @r = ();
if($in_intersection =~ /.gz$/){ open(FINT, "gunzip -c $in_intersection|"); }
else{ open(FINT, "$in_intersection"); }
while(<FINT>){
		chomp;
		my @t = split(/\t/,$_);
		my $region_name = $t[3];
		#my ($context, $ncmd, $landrace) = split(/;/,$t[8]);
		my $ncmd = $t[9];
		my $landrace = $t[10];
		my $qval = $t[11];
		my $methyl_diff = $t[12];
		
		my $methyl_type = "";
		$hs_methyl_level{$region_name}{$t[7]} = $methyl_diff;
		if($ncmd > $landrace){ $methyl_type = "hyper"; }
		elsif($ncmd < $landrace){ $methyl_type = "hypo"; }
		if(!exists($hs_length{$region_name})){
				$hs_region{$region_name} = "$t[0]\t$t[1]\t$t[2]";
				$hs_gene_direction{$region_name} = $t[4];
				push(@r, $region_name);
				my $length = $t[2]-$t[1];
				$hs_length{$region_name} = $length;
				$hs_count{$region_name}{$methyl_type} = 1;
		}
		else{
				$hs_count{$region_name}{$methyl_type} += 1;
		}
}
close(FINT);
open(FSTAT, ">$out_dir/detailed_count.txt");
open(FTSSDIST, ">$out_dir/distance_from_TSS.txt");
open(FDIST, ">$out_dir/distance_between_DMBs.txt");
open(FMETH, ">$out_dir/methylated_position.txt");
print FSTAT "#CHR\tSTART\tEND\tREGION_TYPE\tGENEID\tLENGTH\tC_HYPER\tC_HYPO\tC_TOTAL\tNormC_HYPER\tNormC_HYPO\tNormC_TOTAL\n";
print FTSSDIST "#avg_dist_from_TSS\tnum_total_DMBs\tregion_type\tgene_name\tgene_direction\tdistances\n";
print FDIST "#avg_dist_between_DBMs\tnum_total_DMBs\tregion_type\tgene_name\tgene_direction\tdistances\n";
print FMETH "#num_total_DMBs\tregion_type\tgene_name\tgene_direction\tpositions\n";
for(my $i = 0; $i <= $#r; $i++){
		if(!exists($hs_count{$r[$i]}{"hyper"})){ $hs_count{$r[$i]}{"hyper"} = 0; }
		if(!exists($hs_count{$r[$i]}{"hypo"})){ $hs_count{$r[$i]}{"hypo"} = 0; }
		my $hyper_count = $hs_count{$r[$i]}{"hyper"};
		my $hypo_count = $hs_count{$r[$i]}{"hypo"};
		my $total_count = $hyper_count + $hypo_count;
		my ($region_type, $gene_name) = split(/;/,$r[$i]);
		print FSTAT "$hs_region{$r[$i]}\t$region_type\t$gene_name\t$hs_length{$r[$i]}\t$hyper_count\t$hypo_count\t$total_count";
		my $norm_hyper = $hyper_count/$hs_length{$r[$i]};
		my $norm_hypo = $hypo_count/$hs_length{$r[$i]};
		my $norm_total = $norm_hyper + $norm_hypo;
		print FSTAT "\t$norm_hyper\t$norm_hypo\t$norm_total\n";
		my @dmb_position = keys(%{$hs_methyl_level{$r[$i]}});
#		@dmb_position = natsort(@dmb_position);
		@dmb_position = sort {$a<=>$b} (@dmb_position);
		if($hs_gene_direction{$r[$i]} eq "-"){
				@dmb_position = reverse(@dmb_position);
		}
		my ($chr, $start, $end) = split(/\t/,$hs_region{$r[$i]});
		my $type = "";
		if($region_type =~ /upstream/){ $type = "U"; }
		elsif($region_type =~ /down/){ $type = "D"; }
		else{ $type = "G"; }
		my $string_tss = "";
		my $string_dmb = "";
		my $total_dist_tss = 0;
		my $total_dist_dmb = 0;
		my $prev_dmb = "";
		my $string_methyl = "$type\t$gene_name\t$hs_gene_direction{$r[$i]}";
##################################### Distance between DMBs
		for(my $dp = 0; $dp <= $#dmb_position; $dp++){
				$string_methyl .= "\t$chr:$dmb_position[$dp];$hs_methyl_level{$r[$i]}{$dmb_position[$dp]}";
				if($prev_dmb ne ""){
						my $dist_between_dmb = abs($dmb_position[$dp]-$prev_dmb);
						$string_dmb .= "\t$dist_between_dmb";
						$total_dist_dmb += $dist_between_dmb;
				}
				$prev_dmb = $dmb_position[$dp];
		}
##################################### Distance from TSS
		if($type ne "D"){
				my $tss = 0;
				if($hs_gene_direction{$r[$i]} eq "+"){ 
						if($type eq "G"){ $tss = $start+1; }
						else{ $tss = $end; }
				}
				else{
						if($type eq "G"){ $tss = $end; }
						else{  $tss = $start+1; }
				}
				for(my $dp = 0; $dp <= $#dmb_position; $dp++){
						my $dist_from_tss = abs($dmb_position[$dp]-$tss);
						$string_tss .= "\t$dist_from_tss";
						$total_dist_tss += $dist_from_tss;
				}
		}
#################################### Print output
		my $num_of_dmbs =  scalar(@dmb_position);
		$total_dist_tss /= $num_of_dmbs;
		if($string_dmb eq ""){
				$total_dist_dmb = "-";
				$string_dmb = "-";
		}
		else{
				$total_dist_dmb /= $num_of_dmbs-1;
		}
		print FTSSDIST "$total_dist_tss\t$num_of_dmbs\t$type\t$gene_name\t$hs_gene_direction{$r[$i]}\t$string_tss\n";
		print FDIST "$total_dist_dmb\t$num_of_dmbs\t$type\t$gene_name\t$hs_gene_direction{$r[$i]}\t$string_dmb\n";
		print FMETH "$num_of_dmbs\t$string_methyl\n";
}
close(FSTAT);
close(FTSSDIST);
close(FDIST);
`awk '\$8==0' $out_dir/detailed_count.txt > $out_dir/hyperDMC_detailed_count_methyl.txt`;
`awk '\$7==0' $out_dir/detailed_count.txt > $out_dir/hypoDMC_detailed_count_methyl.txt`;
