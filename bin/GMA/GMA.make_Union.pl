#!/usr/bin/perl
use strict;
use warnings;

#my $in_methyl_pos = shift;
my $ref_name = $ARGV[0];
my $cpu = $ARGV[1];
my $in_chr_list = $ARGV[2];
my $out_dir = $ARGV[-1];

`mkdir -p $out_dir`;
my @arr_all_chr = ();
open(FCHR, "$in_chr_list");
while(<FCHR>){
		chomp;
		my @t = split(/\t/,$_);
		push(@arr_all_chr, $t[0]);
}
close(FCHR);
my @out_path = split(/\//,$out_dir);
my @msr_dir = @out_path;
pop(@msr_dir);
my $methylseekr_dir = join("/",@msr_dir);
$methylseekr_dir .= "/MethylSeekR";
open(FPAOUT, ">$out_dir/params.txt");

open(FALLOUT, ">$out_dir/all_methylCalls.txt");
foreach my $chr (@arr_all_chr){
		my %hs_call = ();
		for(my $i= 3; $i < $#ARGV; $i++){
				my $cur_dir = $ARGV[$i];
				my @dir_path = split(/\//,$cur_dir);
				my $cur_file = "$cur_dir/$dir_path[-3].$chr.txt";
				print STDERR "FILE: $cur_file\n";
				if($cur_file =~ /.gz$/){
						open(FCALL, "gunzip -c $cur_file |");
				}
				else{
						open(FCALL, "$cur_file"); 
				}
				while(<FCALL>){
						chomp;
						if($_ =~ /^chrBase/){ next; }
						my ($t_chr, $base, $end, $perc, $c_count, $t_count) = split(/\s+/,$_);
						my $cc = $c_count;
						my $ct = $t_count;
						if(!exists($hs_call{$chr}{$base}{C})){
								$hs_call{$chr}{$base}{C} = $cc;
								$hs_call{$chr}{$base}{T} = $ct;
						}
						else{
								$hs_call{$chr}{$base}{C} += $cc;
								$hs_call{$chr}{$base}{T} += $ct;
						}
				}
				close(FCALL);
		}
		open(FOUT, ">$out_dir/$chr.txt");
		if($chr !~ /_/){ ## Fixed at Apr 2, 2021
				print FPAOUT "$ref_name\t$cpu\t$chr\t$out_dir/$chr.txt\t$methylseekr_dir\n";
		}
		foreach my $base (sort {$a<=>$b} keys(%{$hs_call{$chr}})){
				if(exists($hs_call{$chr}{$base}{C})){
						my $total_cnt = $hs_call{$chr}{$base}{T} + $hs_call{$chr}{$base}{C};
						if($total_cnt >= 10){
								print FOUT "$chr\t$base\t$total_cnt\t$hs_call{$chr}{$base}{C}\n";
								print FALLOUT "$chr\t$base\t$total_cnt\t$hs_call{$chr}{$base}{C}\n";
						}
				}
		}
		close(FOUT);
}
close(FPAOUT);
close(FALLOUT);
