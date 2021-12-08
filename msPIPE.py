#! /usr/bin/env python3.6

import sys
import os
import subprocess as sub
import configparser
import time

pipe_path = os.path.abspath(sys.argv[0])
pipe_path = os.path.split( pipe_path )[0]
binD = pipe_path + "/bin"
sys.path.append(pipe_path + '/lib')
sys.path.append(pipe_path + '/bin')

from func import *
from methClass import *


def make_argparser():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('--param','-p', type=str, help= 'config format parameter file', metavar = 'params.conf', required=True)
    parser.add_argument('--out','-o', type=str, help = 'output directory', metavar = 'PATH', required=True)
    parser.add_argument('--core', '-c', type=int, metavar='int', help='core (default:5)', default=5)
    parser.add_argument('--qvalue', '-q', type=float, help='q-value cutoff (default:0.5)', default=0.5, metavar = 'float' )
    parser.add_argument('--skip_trimming', help = 'skip the trimgalore trimming', action='store_true')
    parser.add_argument('--skip_mapping', help = 'skip the bismark mapping', action='store_true')
    parser.add_argument('--skip_calling', help = 'skip the methylation calling', action='store_true')
    parser.add_argument('--skip_HMR', help = 'skip the HMR analysis', action='store_true')
    parser.add_argument('--skip_DMR', help = 'skip the DMR analysis', action='store_true')

    
    return parser



parser = make_argparser()
args = parser.parse_args()
paramF = args.param
outD = os.path.abspath(args.out)

## ===== result space =====
curD = os.getcwd()
timestamp = time.strftime('%y%m%d%H%M', time.localtime(time.time()))
log_proc = outD + f'/log.Process_{timestamp}.txt'
if os.path.exists(log_proc):
    O = open(log_proc, 'w')
    O.write('\n')
    O.close()



sub.call(f'mkdir -p {outD}', shell=True)
sub.call(f'cp {paramF} {outD}/', shell=True)
if not os.path.exists(paramF):
    sys.exit()



## =====    get parameters ====

config = configparser.ConfigParser()
config.read(paramF)
config_keys = config.keys()

ANALYSISs = []
LIBs = []
for key in config_keys:
    if key.find("LIB") != -1 :
        LIBs.append( Parameters(config, key) ) 
    
    DMR = Parameters(config, 'DMR')
    REF = Parameters(config, 'REFERENCE')


SAMPLE = {}
## ======   programs    =====

class Programs():
    def __init__(self):
        binD = pipe_path + "/bin"
        self.multiqc = 'multiqc'
        
       # Requirment 
        self.trim_galore = 'trim_galore'
        self.samtools = 'samtools'
        self.bismark = 'bismark'
        self.bismark_indexing = 'bismark_genome_preparation'
        self.bismark_methylation_extractor = 'bismark_methylation_extractor'
        self.bismark_bedgraph = 'bismark2bedGraph'
        self.cutadapt = 'cutadapt'
        
        scriptD = binD + '/script' 
        self.split = scriptD+'/splitF_bychr.py'
        self.getDMC = scriptD + '/getDMC.py'
        self.getStat = scriptD + '/getStat.pl'
        self.dmr = scriptD + '/getDMR_500bp.pl'

        self.gene = binD + '/GMA/Gene-Methyl-Analysis.pl'
        self.bed2wig = scriptD + '/MethLvBEDtoWIG.pl'
        
        visD = f"{pipe_path}/bin/vis_script/"
        self.analCpG = visD + 'visualization_parallel.R'
        self.window = visD + "win100kb_methylLevel.R"
        self.circos = visD + 'GMA.Circos_100kb.R'
        self.contextLev = visD + 'genomic_context_levels.R'


prog = Programs()


## ===== check input file  =====

printlog(log_proc, '- command: '+' '.join(sys.argv) + '\n')
printlog(log_proc, '-Check all input')

for i in range(len(LIBs)):
    lib = LIBs[i]
    exit= False
    check_line =  f'  LIB{i+1}'
    if not lib.lib_type :
        check_line += f'\n*** LIB_TYPE is necessary(P or S)'
        exit = True
    else:
        if not os.path.exists(lib.file_1):
            check_line += f'\n*** FILE_1({lib.file_1}) is not existed'
            exit = True    
        if lib.lib_type == 'S': continue # pass to check file2 when libtype is single end
        if not os.path.exists(lib.file_2):
            check_line += f'\n*** FILE_2({lib.file_2}) is not existed'
            exit = True 
    
    if exit :
        printlog(log_proc, check_line)
        sys.exit()
    else :
        printlog(log_proc, check_line+' .. ok')

## =====    check reference =====
printlog(log_proc, f'\n- Check reference {REF.ucsc_name}')

if not 'ucsc_name' in REF.__dict__ or REF.ucsc_name == "":
    print(f'*** "ucsc_NAME" is mandatory option (UCSC genome name). Please check the parameter file')
    sys.exit()


refD = pipe_path + f'/reference/{REF.ucsc_name}'
if not( os.path.exists(refD)):
    sub.call(f'mkdir -p {refD}', shell=True)        

if 'fasta' in REF.__dict__:
    if not REF.fasta.replace('.gz','').split(".")[-1] in ('fa','fasta'):
        printlog(log_proc, '*** FASTA format must be (.fa, .fa.gz, .fasta or .fasta.gz) file extensions')
        sys.exit()

    fastaF = f'{refD}/{REF.fasta.split("/")[-1]}'
    if not (os.path.exists(fastaF.replace('.gz',''))):
        sub.call(f'cp {REF.fasta} {refD}/', shell=True)
        if fastaF.split('.')[-1] == 'gz':
            sub.call(f'gzip -d {fastaF}', shell=True)
            fastaF = fastaF.replace('.gz','')
        #cmd = f'grep ">" {fastaF} '
        #chr_list = sub.call(cmd, stdout=pipe, shell=True)
    else:
        printlog(log_proc, f'{fastaF} exist..')

else:
    fastaF = f'{refD}/{REF.ucsc_name}.fa'
    if (os.path.exists(fastaF+'.gz')):
        os.remove(f'{fastaF}.gz')
    if not (os.path.exists(f'{fastaF}')):
        printlog(log_proc, f"-- REFERENCe DOWNLOAD ... (fasta)")
        ucsc_return = UCSC_download(REF.ucsc_name, 'fa', refD)
        if ucsc_return !=0:
            printlog(log_proc, f'\n*** reference download ERR')
            printlog(log_proc, f"*** check logfile > {refD}/log.download_fa.txt\n")
            sys.exit()
    else:
        printlog(log_proc, f'{fastaF} exist..')


if 'gtf' in REF.__dict__:
    annotF = f'{refD}/{REF.gtf.split("/")[-1]}'
    if not (os.path.exists(annotF.replace('.gz',''))):
        sub.call(f'cp {REF.gtf} {refD}/', shell=True)
        if fastaF.split('.')[-1] == 'gz':
            sub.call(f'gzip -d {annotF}', shell=True)
            annotF = annotF.replace('.gz','')
    else:
        printlog(log_proc, f'{annotF.replace(".gz","")} exist..')
else:
    annotF = f'{refD}/{REF.ucsc_name}.ncbiRefSeq.gtf'
    if (os.path.exists(annotF+'.gz')):
        os.remove(f'{annotF}.gz')
    if not (os.path.exists(f'{annotF}')):
        printlog(log_proc, f"-- REFERENCe DOWNLOAD ... (gtf)")
        ucsc_return = UCSC_download(REF.ucsc_name, 'gtf', refD)
        
        if ucsc_return !=0:
            printlog(log_proc, f'\n*** reference download ERR')
            printlog(log_proc, f"*** check logfile > {refD}/log.download_gtf.txt\n")
            sys.exit()
    else:
        printlog(log_proc, f'{annotF} exist..')




## Bisulfite genome
if not ( os.path.exists(refD + '/Bisulfite_Genome') ):
    printlog(log_proc, f"\n-- Building the bismark reference genome")
    
    cmd = f'{prog.bismark_indexing} --parallel {args.core} --verbose {refD}'
    subproc(cmd,f'{refD}/log.indexing.txt',log_proc)
    printlog(log_proc, "\tDone.\n")
else:
    printlog(log_proc, f"\n-- Building the bismark reference genome .. skipped")


## =====    processing  & calling====
printlog(log_proc, '\n'+'='*60)
callD = outD + '/methylCALL'
sub.call(f'mkdir -p {callD}', shell=True)
preproc_complete = True

for i in range(len(LIBs)):
    LIB = LIBs[i]
    cur_lib = f'LIB{i+1}'
    printlog(log_proc, f'\nProcessing .. {cur_lib}: {LIB.lib_name}\n')
    
    # result space
    libD = callD +'/'+ LIB.lib_name
    dataD = libD + '/data'
    methylD = libD + '/methylcontext'
    
    chr_CpGD = methylD + '/CpG_chr'
    chr_CHGD = methylD + '/CHG_chr'
    chr_CHHD = methylD + '/CHH_chr'

    sub.call(f'mkdir -p {libD}/logs', shell=True )
    sub.call(f'mkdir -p {dataD}', shell=True )
    sub.call(f'mkdir -p {chr_CpGD}', shell=True )
    sub.call(f'mkdir -p {chr_CHGD}', shell=True )
    sub.call(f'mkdir -p {chr_CHHD}', shell=True )
    
    fastqc_args = f'"-f fastq -o {dataD} -d {libD} -t 6"'
    cmd_trimGalor = f'{prog.trim_galore} --fastqc --fastqc_args {fastqc_args} --phred33 --gzip --length 20 -o {dataD} --cores {args.core}'
    cmd_bismark_align = f'{prog.bismark} --score_min L,0,-0.6 -N 0 -L 20 --parallel {args.core} --temp_dir {dataD}/temp {refD} -o {dataD}'
    
    if LIB.lib_type == 'P':
        cmd_trimGalor += f' --paired {LIB.file_1} {LIB.file_2}'
        cmd_bismark_align += f' -1 {dataD}/*_val_1.fq.gz -2 {dataD}/*_val_2.fq.gz'
         
    elif LIB.lib_type == 'S':
        cmd_trimGalor += f' {LIB.file_1}'
        cmd_bismark_align += f' {dataD}/*_trimmed.fq.gz'
        
    ## fastq file processing
    if args.skip_trimming or args.skip_mapping:
        printlog(log_proc, "1. Input read QC (TrimGalore!) .. skip")
    else:
        ##1. TrimGalore
        printlog(log_proc, "1. Input read QC (TrimGalore!) .. start")
        call_return = subproc(cmd_trimGalor, f'{libD}/logs/log.TrimGalore.txt',log_proc)
    
    if args.skip_mapping :
        printlog(log_proc, '2. Mapping (bismark) .. skip')
    else: 
        ##2. bismark alignment
        printlog(log_proc, '2. Mapping (bismark) .. start')
        call_return = subproc(cmd_bismark_align,f'{libD}/logs/log.bismark.txt',log_proc)
        if call_return != 0:
            preproc_complete = False
            continue

        bamF = f'{dataD}/{LIB.lib_name}_bismark_mapping.bam'
        sub.call(f'mv {dataD}/*.bam {bamF}', shell=True)
        sub.call(f'mv {dataD}/*bismark*_report.txt {dataD}/{LIB.lib_name}_bismark_mapping_report.txt', shell=True)
    
    bamF = f'{dataD}/{LIB.lib_name}_bismark_mapping.bam'

    if not args.skip_calling:
        ##2.5 multiQC
        cmd = f'{prog.multiqc} {dataD} -o {libD}'
        call_return = subproc(cmd, 'stdout',log_proc)
        
        
        ##3.calling:  Methylation extractor
        printlog(log_proc, '\n3. Methylation extractor .. start')
        cmd = f'{prog.bismark_methylation_extractor} -{LIB.lib_type.lower()} --no_overlap --comprehensive --gzip --CX --cytosine_report --genome_folder {refD}/Bisulfite_Genome -o {methylD} {bamF}'
        call_return = subproc(cmd, f'{libD}/logs/log.methylcontext.txt',log_proc)
        if call_return != 0:
            preproc_complete = False
            continue

        ##4. bedGraph
        printlog(log_proc,'4. bedGraph (CpG context) .. start')
        #CpG
        cmd = f'{prog.bismark_bedgraph} -o CpG.cov --dir {methylD} {methylD}/CpG_context_*'
        call_return= subproc(cmd,f'{libD}/logs/log.CpG_bismark2bedGraph.txt',log_proc)
        sub.call( f'gzip -dc {methylD}/CpG.cov.gz.bismark.cov.gz > {libD}/{LIB.lib_name}_CpG.cov.txt', shell=True)
        sub.call( [prog.split, f'{libD}/{LIB.lib_name}_CpG.cov.txt', f'{chr_CpGD}/{LIB.lib_name}'] )

        #CHG
        cmd = f'{prog.bismark_bedgraph} -o CHG.cov --dir {methylD} --CX {methylD}/CHG_context_*'
        call_return= subproc(cmd,f'{libD}/logs/log.CHG_bismark2bedGraph.txt',log_proc)
        sub.call( f'gzip -dc {methylD}/CHG.cov.gz.bismark.cov.gz > {libD}/{LIB.lib_name}_CHG.cov.txt', shell=True)
        sub.call( [prog.split, f'{libD}/{LIB.lib_name}_CHG.cov.txt', f'{chr_CHGD}/{LIB.lib_name}'] )

        #CHH
        cmd = f'{prog.bismark_bedgraph} -o CHH.cov --dir {methylD} --CX {methylD}/CHH_context_*'
        call_return= subproc(cmd,f'{libD}/logs/log.CHH_bismark2bedGraph.txt',log_proc)
        sub.call( f'gzip -dc {methylD}/CHH.cov.gz.bismark.cov.gz > {libD}/{LIB.lib_name}_CHH.cov.txt', shell=True)
        sub.call( [prog.split, f'{libD}/{LIB.lib_name}_CHH.cov.txt', f'{chr_CHHD}/{LIB.lib_name}'] )
        

        # keep bismark CpG result
    if not LIB.sample_name in SAMPLE:
        SAMPLE[LIB.sample_name] = {}
    SAMPLE[LIB.sample_name][LIB.lib_name] = [ f'{libD}/{LIB.lib_name}_CpG.cov.txt', f'{libD}/{LIB.lib_name}_CHG.cov.txt', f'{libD}/{LIB.lib_name}_CHH.cov.txt' ] 




if preproc_complete :
    printlog(log_proc, f'\nPre-processing and calling is done')
    ## multiQC for all data
    cmd = f'{prog.multiqc} {dataD} -o {callD}'
    call_return = subproc(cmd, 'stdout',log_proc)
else :
    printlog("ERR: problem in preprocessing & calling step")
    sys.exit()


#-------------------------------------------------------------------------------------------------------
sample_names = list(SAMPLE.keys())

if len(sample_names) == 0:
    printlog('stdout', 'Error. There is no samples to compare (only one sample existed).\n')
    sys.exit()

## sample_level result

resultD = outD + '/Analysis'
result_logD = resultD + '/logs'
sub.call(f'mkdir -p {result_logD}', shell=True)

printlog(log_proc,'\n\nGene & Methyl Analysis ...'+'-'*60)

##make parameter file
param_file = f'{resultD}/params.txt'

p_line = f'@REF_NAME\n{REF.ucsc_name}\n\n@REF_FA\n{fastaF}\n\n@GENE_GTF\n{annotF}\n\n'
for sample in sample_names:
    p_line += f'@{sample}\n'
    for lib in SAMPLE[sample].keys():
        for F in SAMPLE[sample][lib]:
            p_line += f'{F}\n'

    p_line += '\n'

paramO = open(param_file,'w')
paramO.write(p_line)
paramO.close()



## GMA 
cmd = f'{prog.gene} -p {param_file} -o {resultD} -cpu {args.core}'
subproc(f'{prog.gene} -p {param_file} -o {resultD} -cpu {args.core}', f'{result_logD}/log.genemethylAnal.txt',log_proc)


## VISUALIZATION

## vis_param
vis_param = 'sample\tcontext\tfile\n'

for sample in sample_names:
    for cx in ("CpG", "CHG", "CHH"):
        call_file =  f'{resultD}/{sample}/union_{cx}/all_methylCalls.txt'
        vis_param += f'{sample}\t{cx}\t{call_file}\n'

vis_param_file = f'{resultD}/vis_params.txt'
paramO = open(vis_param_file,'w')
paramO.write(vis_param)
paramO.close()


##make wig, window, histogram
wig_cmd =  f'Rscript {binD}/vis_script/visualization_parallel2.R {vis_param_file} {resultD}/Ideogram.bed {resultD}/ {args.core}'
subproc(wig_cmd, result_logD + f'/log.bed2wigANDwindow', log_proc)

for sample in sample_names:
    sampleD = resultD + f'/{sample}'

    # wig2bigwig
    bw_cmd = f'wigToBigWig {sampleD}/CpG_methylLev.wig {resultD}/ref.size {sampleD}/CpG_methylLev.bw'
    subproc(bw_cmd, result_logD + f'/log.{sample}.wig2bw', log_proc)

    #genomic context visualization
    context_cmd = f'Rscript {prog.contextLev} {sampleD}/methylation.Genomic_Context.CpG.txt {sample} {sampleD}/Genomic_Context_CpG.pdf 1> {sampleD}/Avg_Genomic_Context_CpG.txt'    
    subproc(context_cmd, result_logD+f'/log.{sample}.Genomic_Context', log_proc)

    #circos
    circos_cmd = f'{prog.circos} {resultD}/Ideogram.bed {sampleD}/ {sampleD}/Circos.CpG_UMRs_LMRs.pdf'
    subproc(circos_cmd, result_logD+f'/log.{sample}.Circos', log_proc)
    


## get Promoter
promoterF = f'{resultD}/annotations/promoter.bed'



#------------------------------------------------------------------------------------------------------------


## DMR analysis
printlog(log_proc,'\n\nDMR Analysis ...')
printlog(log_proc, '-'*60)

if DMR.param_values() == [] or DMR.param_values() == [""]:
    print("NO DMR analysis group set")
else:
    dmrD = resultD + '/DMR'
    sub.call(f'mkdir -p {dmrD}', shell=True)

    for AN in DMR.param_values():
        control, case = map(lambda x: x.strip(), AN.split(','))
        file_list = []
        treatment = []
        sample_lib = []
        printlog(log_proc, f'{AN}')
        for lib in SAMPLE[control].keys():
            treatment.append('0')
            file_name = SAMPLE[control][lib][0].split('/')[-1]
            file_list.append(f'"{SAMPLE[control][lib][0]}"' )
            sample_lib.append(f'"{lib}"')
    
        for lib in SAMPLE[case].keys():
            treatment.append('1')
            file_name = SAMPLE[case][lib][0].split('/')[-1]
            file_list.append( f'"{SAMPLE[case][lib][0]}"' )
            sample_lib.append(f'"{lib}"')


        ## methylKit Rscript        ====================================
            analD = f'{dmrD}/{control}.{case}'
            anal_logD = f'{dmrD}/{control}.{case}/logs'
            methylkitD = f'{analD}/methylkit'
            sub.call(f'mkdir -p {methylkitD}', shell=True)
            sub.call(f'mkdir -p {anal_logD}', shell=True)
            
            methylkit_script = f'{methylkitD}/run_methylKit.R'
            
            methF = f'{methylkitD}/CpG_united_filtered_meth.txt'
            rawF = f'{methylkitD}/CpG_united_filtered_raw_DMBs.txt'
    
        
            sample_ids = ','.join(sample_lib)
            treatments = ','.join(treatment)
            methylKO = open(methylkit_script, 'w')
            ## write methylKit script
            line =( f'#!/usr/bin/env Rscript\n' + 
                    f'library(methylKit);\n' +
                    f'file.list=list( {",".join(file_list)});')

            line += f'myobj=methRead(file.list, sample.id=list({sample_ids}), assembly="{REF.ucsc_name}",\n\n'
            line += f'treatment=c({treatments}), context="CpG", pipeline="bismarkCoverage");\n'
            line += f'filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9);\n\n'
            line += f'meth=unite(filtered.myobj, destrand=FALSE);\n'
            line += f'write.table(meth, file=paste("{methF}", sep="\\t"));\n'
            line += f'myDiff=calculateDiffMeth(meth);\n\n'
            line += f'write.table(myDiff, file=paste("{rawF}", sep="\\t"));\n'
            methylKO.write(line)
            methylKO.close()
            subproc(f'Rscript {methylkit_script}', f'{anal_logD}/log.methylkit.txt',log_proc)

        if True:
        ## Methylation Level
            methylLevelF = f"{analD}/DMC_q{args.qvalue}.bed"
            cmd = f'{prog.getDMC} {methF} {rawF} {args.qvalue} {treatments} {methylLevelF} {control} {case}'
            subproc(cmd, f'{anal_logD}/log.DMC.txt',log_proc)
        
        ##get methylation Levels on promoter region
            intersect_promoterF = f'{analD}/intersection.DMC2Promoter.txt'
            cmd_intersect = f'bedtools intersect -wa -wb -a {promoterF} -b {analD}/reform.DMC_q{args.qvalue}.bed > {intersect_promoterF}'
            subproc(cmd_intersect, "stdout",log_proc)
        
            statF = f'{analD}/Stat.CpG_qvalue{args.qvalue}.out'
            cmd_stat = f'{prog.getStat} {intersect_promoterF} {analD}'
            subproc(cmd_stat, f"{anal_logD}/log.getStat.txt",log_proc)
        
        
            if not os.path.isfile(methylLevelF):
                printlog(log_proc,"There is no significantly differentially methylated base.\n")
                printlog(log_proc, "->DMR analysis is skipped.\n")
            else:
                pass
                ## getDMR
                subproc(f'{prog.dmr} {rawF} {args.qvalue} 1> {analD}/DMR.bed', f'{anal_logD}/log.getDMR.txt',log_proc)
            





