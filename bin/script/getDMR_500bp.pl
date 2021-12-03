#!/usr/bin/perl
#
#
use strict;
use warnings;
use Sort::Key::Natural 'natsort';

my $f_in = $ARGV[0];
my $qval_cutoff = $ARGV[1];

my $cnt_500 = 0;
my $DMR_st = 0;
my $DMR_st_trim = 0;
my %hs_DMR = ();#{startidx}
my %hs_idx_info = ();
my $chr_bf = "init";
my $pos_bf = 0;
my $md_bf = 1;

open(F,"$f_in");
while(<F>){
	chomp;
	if ($_ =~ /^"chr"/){ #header
		next;
	}else{
		my $trimmed = $_;
		$trimmed =~ tr/\"\\\///d;
		my ($idx, $chr,$st, $ed, $str, $pv, $qv, $md) = split(/\s+/,$trimmed);
		$hs_idx_info{"$chr:$st"} = $idx;
		if ($DMR_st == 0){
			$cnt_500 = 1;
			$DMR_st = $st;
			$DMR_st_trim = $trimmed;
		
		}else{
			if ($chr_bf eq $chr){
				$cnt_500 += $st - $pos_bf;
			}else{
				$DMR_st = $st;
				$DMR_st_trim = $trimmed;
				$cnt_500 = 1;
				$chr_bf = $chr;
				$pos_bf = $st;
				$md_bf = $st;
				next;
			}
			if($cnt_500 <= 500 && $chr_bf eq $chr && $md_bf * $md > 0){ #edited for +/- meth.diff level 
				#print STDERR "find DMR! $DMR_st, $st, distance: $cnt_500\n";
				if (!exists($hs_DMR{$hs_idx_info{"$chr:$DMR_st"}})){
					($idx, $chr,$st, $ed, $str, $pv, $qv, $md) = split(/\s+/,$DMR_st_trim);
					if ($qv<= $qval_cutoff){
						$hs_DMR{$hs_idx_info{"$chr:$DMR_st"}} = $DMR_st_trim."\n";
					}
					($idx, $chr,$st, $ed, $str, $pv, $qv, $md) = split(/\s+/,$trimmed);
					if ($qv<= $qval_cutoff){
						$hs_DMR{$hs_idx_info{"$chr:$DMR_st"}} .= $trimmed."\n";
					}

				}else{
					if ($qv<= $qval_cutoff){
						$hs_DMR{$hs_idx_info{"$chr:$DMR_st"}} .= $trimmed."\n";
					}
				}
				$cnt_500 = 1;
			}else{
				$DMR_st = $st;
				$DMR_st_trim = $trimmed;
				$cnt_500 = 1;
			}

		}
		$chr_bf = $chr;
		$pos_bf = $st;
		$md_bf = $md;
	}
}
close(F);
	

print "chr\tstart\tend\tcase-control\tdmc_count\tavg_qvalue\n";

foreach my $pos (natsort keys %hs_DMR){
	#print STDERR $pos."\n";
	#print STDERR $hs_DMR{$pos}."\n";
	my $chr = "";
	my $from = 999999999;
	my $to = 0;
	my $mdpn = ""; #md positive-negative 
	my $dmc_cnt = 0	;
	
	my $qv =0;
	my $qv_sum =0;
	my $qv_avg=0;

	foreach my $record (natsort split(/\n/,$hs_DMR{$pos})){
		my @tmp_ar = split(/\s+/, $record);
		if(@tmp_ar==0){next;}
		$chr = $tmp_ar[1];
		$from = $tmp_ar[2] < $from ? $tmp_ar[2] : $from;
		$to = $tmp_ar[2] > $to ? $tmp_ar[2]:$to;
		$qv =  $tmp_ar[6];
		if ( $tmp_ar[7]>0){
			$mdpn = "positive";
		}else{
			$mdpn = "negative";
		}		
		$dmc_cnt ++;
		$qv_sum += $qv;
	}
	$qv_avg = $qv_sum/$dmc_cnt;
	print $chr."\t".($from-1)."\t".$to."\t".$mdpn."\t".$dmc_cnt."\t".$qv_avg."\n";
	
}





