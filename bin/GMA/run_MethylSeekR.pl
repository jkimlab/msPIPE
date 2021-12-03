#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $in_param = $ARGV[0];
my $out_dir = $ARGV[-1];

`mkdir -p $out_dir`;
$out_dir = abs_path("$out_dir");
open(FPARAM, "$in_param");
while(<FPARAM>){
		chomp;
		my ($ref_name, $cpu, $prefix, $file) = split(/\t/,$_);
		if(!-f $file){
				print STDERR "Please check file [ $file ] \n";
				exit();
		}
		$file = abs_path($file);
		`Rscript $Bin/run_MethylSeekR.R $ref_name $cpu $prefix $file $out_dir`;
		`mv Rplots.pdf $out_dir/$prefix.plots.pdf`;
}
close(FPARAM);
