#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use Cwd 'abs_path';


my $in_param;
my $out_dir = "./output";
my $range = 1500;
my $bin_size = 500;
my $step_size = 100;
my $cpu = 1;
my $help;

my $convert2bed = "convert2bed";
my $bedtools = "bedtools";
my $annotation = $Bin."/GMA.annotation.py";

my $options = GetOptions(
		"param|p=s" => \$in_param,
		"outdir|o=s" => \$out_dir,
		"range_tss|r=i" => \$range,
		"bin=i" => \$bin_size,
		"step=i" => \$step_size,
		"cpu|t=i" => \$cpu,
		"h|help" => \$help,
);

if(!defined($in_param)){ HELP(); }
if(defined($help)){ HELP(); }
`mkdir -p $out_dir/annotations`;
my $annot_dir = abs_path("$out_dir/annotations");
$out_dir = abs_path("$out_dir");
open(FLOG, ">$out_dir/log.gma.txt");
print FLOG "[GMA] output directory: $out_dir\n";
print FLOG "[GMA] annotation fils: $annot_dir\n";

##### Checking parameters
print STDERR "[GMA] Preprocessing\n";
print FLOG "[GMA] Preprocessing\n";
my $data_type = "";
my $ref_name = "";
my $ref_file = "";
my %hs_input = ();
open(FPARAM, "$in_param");
while(<FPARAM>){
		chomp;
		if($_ =~ /^@(\S+)/){
				$data_type = $1;
		}
		elsif($_ eq ""){ next; }
		elsif($_ =~ /^#/){ next; }
		else{
				my $file = $_;
				if($data_type eq "GENE_GTF"){
						if(!-f $file){
								print STDERR "[GMA] Please check gene annotation file\n($file)\n";
								print FLOG "[GMA] Please check gene annotation file\n($file)\n";
								exit();
						}
						print STDERR "[GMA] processing annotation fils \n";
						print "$annotation $file $annot_dir\n";
						`$annotation $file $annot_dir`;
						print"... done\n";

						## merge genes
						`$bedtools sort -i $annot_dir/gene.bed > $annot_dir/sorted.gene.bed`;
						`mv $annot_dir/sorted.gene.bed $annot_dir/gene.bed`;
						`$bedtools merge -i $annot_dir/gene.bed > $annot_dir/merge.gene.bed`;
						
						## merge exons
						`$bedtools sort -i $annot_dir/exon.bed > $annot_dir/sorted.exon.bed`;
						`mv $annot_dir/sorted.exon.bed $annot_dir/exon.bed`;
						`$bedtools merge -i $annot_dir/exon.bed > $annot_dir/merge.exon.bed`;
						
						## make intron regions
						`$bedtools subtract -a $annot_dir/merge.gene.bed -b $annot_dir/merge.exon.bed > $annot_dir/intron.bed`;

				}
				elsif($data_type eq "REF_FA"){
						$ref_file = $file;
				}
				elsif($data_type eq "REF_NAME"){
						$ref_name = $file;
				}
				else{
						my $context_type = "";
						my @tmp_file_path = split(/\//,$file);
						pop(@tmp_file_path);
						my $dir_path = join("/",@tmp_file_path);
						if($file =~ /CpG.cov.txt$/){
								$context_type = "CpG";
								$dir_path .= "/methylcontext/CpG_chr";
						}
						elsif($file =~ /CHG.cov.txt$/){
								$context_type = "CHG";
								$dir_path .= "/methylcontext/CHG_chr";
						}
						elsif($file =~ /CHH.cov.txt$/){
								$context_type = "CHH";
								$dir_path .= "/methylcontext/CHH_chr";
						}
						if(!exists($hs_input{$data_type}{$context_type})){
								$hs_input{$data_type}{$context_type} = "$dir_path";
						}
						else{
								$hs_input{$data_type}{$context_type} .= " $dir_path";
						}
				}
		}
}
close(FPARAM);
`faSize -detailed $ref_file > $out_dir/ref.size`;
`grep -v "_" $out_dir/ref.size | awk '{print \$1 \"\t0\t\" \$2 \"\t.\tgvar\"}' - > $out_dir/Ideogram.bed`;
`grep -v "_" $out_dir/ref.size | awk '{print \$1}' - > $out_dir/chr_list.txt`;
`$bedtools subtract -a $out_dir/Ideogram.bed -b $annot_dir/merge.gene.bed | cut -f1,2,3 - > $annot_dir/intergenic.bed`;

`$Bin/GMA.Add_promoter_into_bed.pl $annot_dir/gene.bed > $annot_dir/tmp.Promoter_Genebody.bed`;
`bedtools sort -i $annot_dir/tmp.Promoter_Genebody.bed > $annot_dir/Promoter_Genebody.bed`;
`rm -f $annot_dir/tmp.Promoter_Genebody.bed`;

`grep "upstream" $annot_dir/Promoter_Genebody.bed > $annot_dir/tmp.promoter.bed`;
`$bedtools sort -i $annot_dir/tmp.promoter.bed > $annot_dir/promoter.bed`;
`rm -f $annot_dir/tmp.promoter.bed`;
`$bedtools merge -i $annot_dir/promoter.bed > $annot_dir/merge.promoter.bed`;

`awk '{print \$1 \"\t\" \$2 \"\t\" \$3  \"\tpromoter\t.\t.\"}' $annot_dir/merge.promoter.bed > $annot_dir/Genomic_context.bed`;
`awk '{print \$1 \"\t\" \$2 \"\t\" \$3  \"\tgene\t.\t.\"}' $annot_dir/merge.gene.bed >> $annot_dir/Genomic_context.bed`;
`awk '{print \$1 \"\t\" \$2 \"\t\" \$3  \"\texon\t.\t.\"}' $annot_dir/merge.exon.bed >> $annot_dir/Genomic_context.bed`;
`awk '{print \$1 \"\t\" \$2 \"\t\" \$3  \"\tintron\t.\t.\"}' $annot_dir/intron.bed >> $annot_dir/Genomic_context.bed`;
`awk '{print \$1 \"\t\" \$2 \"\t\" \$3  \"\tintergenic\t.\t.\"}' $annot_dir/intergenic.bed >> $annot_dir/Genomic_context.bed`;
`$bedtools sort -i $annot_dir/Genomic_context.bed > $annot_dir/sort.Genomic_context.bed`;
my $chr_list = "$out_dir/chr_list.txt";

print STDERR "[GMA] Methylation analysis\n";
print FLOG "[GMA] Methylation analysis\n";
print STDERR "Rscript $Bin/check_MethylSeekR_avail.R $ref_name 2> $out_dir/log.check_MethylSeekR_avail.txt\n";
my $check_msr = `Rscript $Bin/check_MethylSeekR_avail.R $ref_name 2> $out_dir/log.check_MethylSeekR_avail.txt`;
my ($tmp1, $c_msr) = split(/\s+/,$check_msr);
my @arr_sample = (sort keys(%hs_input));
print STDERR "\tall samples: @arr_sample\n";
foreach my $data_type (sort keys(%hs_input)){
		print STDERR "\t$data_type\n";
		print FLOG "\t$data_type\n";
		`mkdir -p $out_dir/$data_type`;

		print STDERR "\t$data_type - CpG context analysis \n";
		print FLOG "\t$data_type - CpG context analysis \n";
		print STDERR "\t$Bin/GMA.make_Union.pl $ref_name $cpu $chr_list $hs_input{$data_type}{CpG} $out_dir/$data_type/union_CpG\n";
		`$Bin/GMA.make_Union.pl $ref_name $cpu $chr_list $hs_input{$data_type}{CpG} $out_dir/$data_type/union_CpG`;

		print STDERR "\t$Bin/GMA.Union_txt2bed.pl $out_dir/$data_type/union_CpG/all_methylCalls.txt > $out_dir/$data_type/CpG_methylCalls.bed\n";
		`$Bin/GMA.Union_txt2bed.pl $out_dir/$data_type/union_CpG/all_methylCalls.txt > $out_dir/$data_type/CpG_methylCalls.bed`;
		
		`$bedtools intersect -wo -a $annot_dir/Genomic_context.bed -b $out_dir/$data_type/CpG_methylCalls.bed > $out_dir/$data_type/methylation.Genomic_Context.CpG.txt`;


		print STDERR "\t$bedtools intersect -wo -a $annot_dir/Promoter_Genebody.bed -b $out_dir/$data_type/CpG_methylCalls.bed > $out_dir/$data_type/methylation_in_Promoter_Genebody.CpG.bed\n";
		`$bedtools intersect -wo -a $annot_dir/Promoter_Genebody.bed -b $out_dir/$data_type/CpG_methylCalls.bed > $out_dir/$data_type/methylation_in_Promoter_Genebody.CpG.bed`;
		
		print STDERR "\t$Bin/GMA.Calculate_avgMethyl_Promoter_Genebody.pl $out_dir/$data_type/methylation_in_Promoter_Genebody.CpG.bed  $annot_dir/gene.bed > $out_dir/$data_type/Average_methyl_lv.txt\n";
		`$Bin/GMA.Calculate_avgMethyl_Promoter_Genebody.pl $out_dir/$data_type/methylation_in_Promoter_Genebody.CpG.bed  $annot_dir/gene.bed > $out_dir/$data_type/Average_methyl_lv.txt`;
		
		`mkdir -p $out_dir/$data_type/AroundTSS`;
		print STDERR "\tTSS +,- $range | bin size $bin_size | step size $step_size\n";
		
		print STDERR "\t$Bin/GMA.make_bins.pl -input $annot_dir/gene.bed -distance $range -bin_size $bin_size -step_size $step_size > $out_dir/$data_type/AroundTSS/tmp.gene.bin.bed\n";
		`$Bin/GMA.make_bins.pl -input $annot_dir/gene.bed -distance $range -bin_size $bin_size -step_size $step_size > $out_dir/$data_type/AroundTSS/tmp.gene.bin.bed`;
		
		print STDERR "\t$bedtools sort -i $out_dir/$data_type/AroundTSS/tmp.gene.bin.bed > $out_dir/$data_type/AroundTSS/gene.bin.bed\n";
		`$bedtools sort -i $out_dir/$data_type/AroundTSS/tmp.gene.bin.bed > $out_dir/$data_type/AroundTSS/gene.bin.bed`;
		`rm -f $out_dir/$data_type/AroundTSS/tmp.gene.bin.bed`;
		
		print STDERR "\t$bedtools intersect -wo -a $out_dir/$data_type/AroundTSS/gene.bin.bed -b $out_dir/$data_type/CpG_methylCalls.bed > $out_dir/$data_type/AroundTSS/methylation_in_bin.bed\n";
		`$bedtools intersect -wo -a $out_dir/$data_type/AroundTSS/gene.bin.bed -b $out_dir/$data_type/CpG_methylCalls.bed > $out_dir/$data_type/AroundTSS/methylation_in_bin.bed`;
		
		print STDERR "\t$Bin/GMA.Calculate_methylLv_aroundTSS.pl $out_dir/$data_type/AroundTSS/methylation_in_bin.bed $annot_dir/gene.bed $range $step_size $bin_size > $out_dir/$data_type/AroundTSS/meth_lv.$data_type.txt\n";
		`$Bin/GMA.Calculate_methylLv_aroundTSS.pl $out_dir/$data_type/AroundTSS/methylation_in_bin.bed $annot_dir/gene.bed $range $step_size $bin_size > $out_dir/$data_type/AroundTSS/meth_lv.$data_type.txt`;
		
		if($c_msr == 1){
				print STDERR "[GMA] UMRs and LMRs analysis\n";
				print FLOG "[GMA] UMRs and LMRs analysis\n";
				print FLOG "\t$Bin/run_MethylSeekR.pl $out_dir/$data_type/union_CpG/params.txt $out_dir/$data_type/MethylSeekR\n";
				print STDERR "\t$Bin/run_MethylSeekR.pl $out_dir/$data_type/union_CpG/params.txt $out_dir/$data_type/MethylSeekR\n";
				`$Bin/run_MethylSeekR.pl $out_dir/$data_type/union_CpG/params.txt $out_dir/$data_type/MethylSeekR`;

				`ls $out_dir/$data_type/MethylSeekR/*.tsv > $out_dir/$data_type/UMRsLMRs.list.txt`;
				print STDERR "\t$Bin/GMA.convert_umrlmr2bed.pl $out_dir/$data_type/UMRsLMRs.list.txt $out_dir/$data_type\n";
				`$Bin/GMA.convert_umrlmr2bed.pl $out_dir/$data_type/UMRsLMRs.list.txt $out_dir/$data_type`;
				`bedtools sort -i $out_dir/$data_type/UMRs.bed > $out_dir/$data_type/sort.UMRs.bed`;
				`bedtools sort -i $out_dir/$data_type/LMRs.bed > $out_dir/$data_type/sort.LMRs.bed`;
				`bedtools intersect -wa -c -a $annot_dir/promoter.bed -b $out_dir/$data_type/sort.UMRs.bed > $out_dir/$data_type/UMR-Promoter.cnt.bed`;
				`bedtools intersect -wo -a $annot_dir/promoter.bed -b $out_dir/$data_type/sort.UMRs.bed > $out_dir/$data_type/UMR-Promoter.pos.bed`;
		}
		else{
				print STDERR "[GMA] (skipped) UMRs and LMRs analysis\n";
				print FLOG "[GMA] (skipped) UMRs and LMRs analysis\n";
				print STDERR "$ref_name is not available for MethylSeekR. (UMR, LMR analysis are skipped)\n";
		}
	
		print STDERR "[GMA] $data_type - CHG context analysis\n";
		print FLOG "[GMA] $data_type - CHG context analysis\n";
		print STDERR "\t$Bin/GMA.make_Union.pl $ref_name $cpu $chr_list $hs_input{$data_type}{CHG} $out_dir/$data_type/union_CHG\n";
		`$Bin/GMA.make_Union.pl $ref_name $cpu $chr_list $hs_input{$data_type}{CHG} $out_dir/$data_type/union_CHG`;

		print STDERR "\t$Bin/GMA.Union_txt2bed.pl $out_dir/$data_type/union_CHG/all_methylCalls.txt > $out_dir/$data_type/CHG_methylCalls.bed\n";
		`$Bin/GMA.Union_txt2bed.pl $out_dir/$data_type/union_CHG/all_methylCalls.txt > $out_dir/$data_type/CHG_methylCalls.bed`;
		
		`$bedtools intersect -wo -a $annot_dir/Genomic_context.bed -b $out_dir/$data_type/CHG_methylCalls.bed > $out_dir/$data_type/methylation.Genomic_Context.CHG.txt`;

		print STDERR "[GMA] $data_type - CHH context analysis\n";
		print FLOG "[GMA] $data_type - CHH context analysis\n";
		print STDERR "\t$Bin/GMA.make_Union.pl $ref_name $cpu $chr_list $hs_input{$data_type}{CHH} $out_dir/$data_type/union_CHH\n";
		`$Bin/GMA.make_Union.pl $ref_name $cpu $chr_list $hs_input{$data_type}{CHH} $out_dir/$data_type/union_CHH`;

		print STDERR "\t$Bin/GMA.Union_txt2bed.pl $out_dir/$data_type/union_CHH/all_methylCalls.txt > $out_dir/$data_type/CHH_methylCalls.bed\n";
		`$Bin/GMA.Union_txt2bed.pl $out_dir/$data_type/union_CHH/all_methylCalls.txt > $out_dir/$data_type/CHH_methylCalls.bed`;
		
		`$bedtools intersect -wo -a $annot_dir/Genomic_context.bed -b $out_dir/$data_type/CHH_methylCalls.bed > $out_dir/$data_type/methylation.Genomic_Context.CHH.txt`;
}

close(FLOG);

sub HELP{
		my $src = basename($0);
		print "\nUsage: $src -p <params.txt> -o <output dir> -t cpu\n";
		print "options\n";
		print "== Mendatory\n";
		print "-param | -p : input parameter files\n";
		print "-outdir | -o : output directory name or path\n";
		print "== Optional\n";
		print "-cpu | -t : number of threads (default: 1)\n";
		print "-range_tss | -r : distance from TSS where you want to see the methylaion level (default: 1500)\n";
		print "-bin : the window size of each bin (default: 500)\n";
		print "-step | the step size of sliding window (default: 100)\n";
		exit();
}
