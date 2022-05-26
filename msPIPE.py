#! /usr/bin/env python3

import sys
import os
import subprocess as sub
import configparser
import argparse
import time
from multiprocessing import Pool

pipe_path = os.path.split( os.path.abspath(sys.argv[0]) ) [0]
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
    parser.add_argument('--qvalue', '-q', type=float, metavar = 'float', help='q-value cutoff (default:0.5)', default=0.5)
    parser.add_argument('--skip_trimming', help = 'skip the trimgalore trimming', action='store_true')
    
    parser.add_argument('--program', type=str, metavar = 'bismark or bs2', help = 'program option for mapping & calling', default='bismark')
    parser.add_argument('--bsmooth', help = 'use bsmooth for DMR analysis', action='store_true')
    
    
    parser.add_argument('--skip_mapping', help = 'skip the bismark mapping', action='store_true')
    parser.add_argument('--skip_calling', help = 'skip the methylation calling', action='store_true')
    parser.add_argument('--calling_data','-m', type=str, help = 'methylCALL directory', metavar = 'PATH')

    parser.add_argument('--skip_GMA', help = 'skip the Gene-Methyl analysis', action='store_true')
    parser.add_argument('--skip_DMR', help = 'skip the DMR analysis', action='store_true')
    
    
    return parser


## ===== result space =====
def make_result_space(outD):
    global log_proc
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
        print(f'no {paramF}')
        sys.exit()

    #set call path
    if args.__dict__['calling_data'] == None:
        args.__dict__['calling_data'] = outD + '/methylCALL'
    

## =====    get parameters ====
def get_parameters(paramF):
    global LIBs
    global DMR
    global REF


    ## parsing config file
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


##  ====    check input and reference files ====
def check_inputs(refD):
    global LIBs
    global REF
    global pipe_path
    global args
    

    ## check program option
    if not args.program  in ['bismark', 'bs2']:
        printlog(log_proc, f'\n*** --program option should be "bismark" or "bs2" (default: bismark))')
        sys.exit()

    
    outD = os.path.abspath(args.out)

    printlog(log_proc, '- command: '+' '.join(sys.argv) + '\n')
    printlog(log_proc, '-Check all input')
    

    ## check input file
    exit= False
    check_line = ''
    for i in range(len(LIBs)):
        lib = LIBs[i]
        if not lib.lib_type :
            check_line += f'\n*** LIB_TYPE is necessary(P or S)'
            exit = True

        else:
            
            if not (args.skip_trimming or args.skip_mapping or args.skip_calling): # input= LIB.FILE1 (fastq)
                if not os.path.exists(lib.file_1) :
                    check_line += f'\n*** FILE_1({lib.file_1}) is not existed'
                    exit = True    
                if lib.lib_type == 'S': continue # pass to check file2 when libtype is single end
                if not os.path.exists(lib.file_2):
                    check_line += f'\n*** FILE_2({lib.file_2}) is not existed'
                    exit = True 
            else:
                dataD = f'{args.calling_data}/{lib.lib_name}/data'

                # ----------------------
                if args.skip_trimming: # skip_trimming, input = LIB.TRIMMED_FILE1,2:
                    if not 'trimmed_file_1' in lib.param_keys():    # When there's no input trimming file, search for pipeline output
                        if lib.lib_type == 'P':
                            lib.trimmed_file_1 = f'{dataD}/{lib.lib_name}_val_1.fq.gz'
                            lib.trimmed_file_2 = f'{dataD}/{lib.lib_name}_val_2.fq.gz'
                        
                        else:
                            lib.trimmed_file_1 = f'{dataD}/{lib.lib_name}_trimmed.fq.gz'
                
                    if not os.path.exists(lib.trimmed_file_1) :
                        ## Use input in parameter file
                        lib.trimmed_file_1 = lib.file_1

                        if not os.path.exists(lib.trimmed_file_1) :
                            check_line += f'\n*** TRIMMED_FILE_1({lib.trimmed_file_1}) is not existed'
                            exit = True    

                    if lib.lib_type == 'S': continue # pass to check file2 when libtype is single end
                    
                    if not os.path.exists(lib.trimmed_file_2):
                        lib.trimmed_file_2 = lib.file_2

                        if not os.path.exists(lib.trimmed_file_2) :
                            check_line += f'\n*** TRIMMED_FILE_2({lib.trimmed_file_2}) is not existed'
                            exit=True

                if args.skip_calling: # skip_calling, input = methylCALL
                    ## check output path
                    if not os.path.exists(dataD):
                        check_line += f'\n*** data directory({dataD}) is not existed'
                        exit = True 
                        continue
                    else:            
                        datadir_list = os.listdir( dataD )
                    
                    if not os.path.exists(args.calling_data):
                        check_line += f'\n*** calling data ({args.calling_data}) is not existed (skip_calling)'
                        exit = True
                    else:
                        continue
                
                elif args.skip_mapping: # skip_mapping, input = LIB.BAM_FILE
                    ## check output path
                    if not os.path.exists(dataD):
                        check_line += f'\n*** data directory({dataD}) is not existed'
                        exit = True 
                        continue
                    else:            
                        datadir_list = os.listdir( dataD )
                    
                    if not 'bam_file' in lib.param_keys():    # When there's no input bam file, search for pipeline output
                        lib.bam_file = f'{dataD}/{lib.lib_name}_*.bam'
                        for F1 in datadir_list:
                            if F1.endswith('bam'):
                                lib.bam_file = f'{dataD}/{F1}'

                    if not os.path.exists(lib.bam_file):
                        print('no bam')
                        check_line += f'\n*** BAM_FILE({lib.bam_file}) is not existed (skip_mapping)'
                        exit = True 
                    
                    continue 
                

    
    if exit :
        printlog(log_proc, check_line)
        sys.exit()
    else :
        printlog(log_proc, check_line+' .. ok')
        

    ## check reference DIR
    printlog(log_proc, f'\n- Check reference {REF.ucsc_name}')

    if not 'ucsc_name' in REF.__dict__ or REF.ucsc_name == "":
        print(f'*** "UCSC_NAME" is mandatory option (UCSC genome name). Please check the parameter file')
        sys.exit()


    if not( os.path.exists(refD)):
        sub.call(f'mkdir -p {refD}', shell=True)        

    ## check fasta file
    if 'fasta' in REF.__dict__:
        ## file format
        if not REF.fasta.replace('.gz','').split(".")[-1] in ('fa','fasta'):
            printlog(log_proc, '*** FASTA format must be (.fa, .fa.gz, .fasta or .fasta.gz) file extensions')
            sys.exit()

        fastaF = f'{refD}/{REF.fasta.split("/")[-1]}'
        if not (os.path.exists(fastaF.replace('.gz',''))):
            sub.call(f'cp {REF.fasta} {refD}/', shell=True)
            if fastaF.split('.')[-1] == 'gz':
                sub.call(f'gzip -d {fastaF}', shell=True)
                fastaF = fastaF.replace('.gz','')
        else:
            printlog(log_proc, f'{fastaF} exist..')

    else:
        fastaF = f'{refD}/{REF.ucsc_name}.fa'
        if (os.path.exists(fastaF+'.gz')):
            os.remove(f'{fastaF}.gz')
        if not (os.path.exists(f'{fastaF}')):
            printlog(log_proc, f"-- REFERENCE DOWNLOAD ... (fasta)")
            ucsc_return = UCSC_download(REF.ucsc_name, 'fa', refD)
            if ucsc_return !=0:
                printlog(log_proc, f'\n*** reference download ERR')
                printlog(log_proc, f"*** check logfile > {refD}/log.download_fa.txt\n")
                sys.exit()
        else:
            printlog(log_proc, f'{fastaF} exist..')

    ## check annotation file
    if 'gtf' in REF.__dict__:
        annotF = f'{refD}/{REF.gtf.split("/")[-1]}'
        if not (os.path.exists(annotF.replace('.gz',''))):
            sub.call(f'cp {REF.gtf} {refD}/', shell=True)
            if annotF.split('.')[-1] == 'gz':
                sub.call(f'gzip -d {annotF}', shell=True)
                annotF = annotF.replace('.gz','')
        else:
            printlog(log_proc, f'{annotF.replace(".gz","")} exist..')
    else:
        annotF = f'{refD}/{REF.ucsc_name}.ncbiRefSeq.gtf'
        if (os.path.exists(annotF+'.gz')):
            os.remove(f'{annotF}.gz')
        if not (os.path.exists(f'{annotF}')):
            printlog(log_proc, f"-- REFERENCE DOWNLOAD ... (gtf)")
            ucsc_return = UCSC_download(REF.ucsc_name, 'gtf', refD)
        
            if ucsc_return !=0:
                printlog(log_proc, f'\n*** reference download ERR')
                printlog(log_proc, f"*** check logfile > {refD}/log.download_gtf.txt\n")
                sys.exit()
        else:
            printlog(log_proc, f'{annotF} exist..')
    
    REF.fastaF = fastaF
    REF.annotF = annotF

    return fastaF, annotF


## Bismark indexing -> make bisulfite genome
def REF_bisulfite_genome(refD):
    global prog
    global args
    global REF

    if args.program == 'bismark':
        ref_genome_dir = refD + '/Bisulfite_Genome' 
        ref_building_cmd = f'{prog.bismark_indexing} --parallel {args.core} --verbose {refD}'

    elif args.program == 'bs2':
        ref_genome_dir = refD + f'/{REF.ucsc_name}.fa_bowtie' 
        ref_building_cmd = f'python2 {prog.bs2_build} --aligner=bowtie --db {refD} -f {REF.fastaF} 1> {refD}/log.indexing.txt '
        
    if not ( os.path.exists( ref_genome_dir ) ):
        printlog(log_proc, f"\n-- Building the {args.program} reference genome")
        subproc(ref_building_cmd,f'{refD}/log.indexing.txt',log_proc)

    else:
        printlog(log_proc, f"\n-- Building the {args.program} reference genome .. skipped")



##  ==== preprocess the reads ====
def Reads_preprocessing(callD):
    global prog
    global args
    global LIBs

    preproc_complete = True
    
    if args.skip_trimming or args.skip_mapping or args.skip_calling:
        printlog(log_proc, "1. Input read QC (TrimGalore!) .. skip")
        return preproc_complete

    else:
        ##1. TrimGalore
        printlog(log_proc, "1. Input read QC (TrimGalore!) .. start")
    
    ## set multi job and core
    job_number, core = set_job_number(len(LIBs), args.core)
    CMD_list = []
    for i in range(len(LIBs)):
        LIB = LIBs[i]
        cur_lib = f'LIB{i+1}'

        # result space
        libD = callD +'/'+ LIB.lib_name
        dataD = libD + '/data'
        sub.call(f'mkdir -p {libD}/logs', shell=True )
        sub.call(f'mkdir -p {dataD}', shell=True )
    
    
        fastqc_args = f'"-f fastq -o {dataD} -d {libD} -t 6"'
        cmd_trimGalor = f'{prog.trim_galore} --fastqc --fastqc_args {fastqc_args} --phred33 --gzip --length 20 -o {dataD} --cores {core}'

        if LIB.lib_type == 'P':
            cmd_trimGalor += f' --paired {LIB.file_1} {LIB.file_2}'
        elif LIB.lib_type == 'S':
            cmd_trimGalor += f' {LIB.file_1}'
        
        CMD_list.append( [cmd_trimGalor, f'{libD}/logs/log.TrimGalore.txt',log_proc ] )
        

    with Pool(int(job_number)) as pool:
        multi_returns = pool.map(multi_run_wrapper, CMD_list)
        if not sum(multi_returns) == 0:
            printlog(log_proc, f'READS PREPROCESSING(trimming) IS NOT COMPLETED.')
            return False
        else:
            printlog(log_proc, '\n\t.. done\n')
            for i in range(len(LIBs)):
                LIB = LIBs[i]
                if LIB.lib_type == 'P':
                    LIB.trimmed_file_1 = f'{dataD}/*_val_1.fq.gz'
                    LIB.trimmed_file_2 = f'{dataD}*_val_2.fq.gz'
                else:
                    LIB.trimmed_file_1 = f'{dataD}/*_trimmed.fq.gz'
            
            return True


##2. bismark alignment  ----------------------------------------------------------------------------
def Bismark_mapping(callD, refD):
    global prog
    global args
    global LIBs
    global REF

    if args.skip_mapping or args.skip_calling:
        printlog(log_proc, f'2. Mapping ({args.program}) .. skip')
        return True
    else: 
        printlog(log_proc, f'2. Mapping ({args.program}) .. start')
    ## set multi job and core
    job_number, core = set_job_number(len(LIBs), args.core)
    CMD_list = []
    sortCMD_list = []
    
    for i in range(len(LIBs)):
        LIB = LIBs[i]
        cur_lib = f'LIB{i+1}'
        
        # result space
        libD = callD +'/'+ LIB.lib_name
        dataD = libD + '/data'
        tmpD = dataD + '/temp'
        sub.call(f'mkdir -p {libD}/logs', shell=True )
        sub.call(f'mkdir -p {tmpD}', shell=True )
        
        bamF = f'{dataD}/{LIB.lib_name}_{args.program}.bam'
        sorted_bamF = f'{dataD}/{LIB.lib_name}_{args.program}.sorted.bam'
        sortCMD_list.append( [f'{prog.samtools} sort -o {sorted_bamF} -@ {core} {bamF}', f'{libD}/logs/log.sortbam.txt',log_proc ] )
    
        # 1)bismark
        if args.program == 'bismark':
            cmd_align = f'{prog.bismark} --score_min L,0,-0.6 -N 0 -L 20 --parallel {core} --temp_dir {tmpD} {refD} -o {dataD}'
        
            if LIB.lib_type == 'P':
                cmd_align += f' -1 {LIB.trimmed_file_1} -2 {LIB.trimmed_file_2}'
            elif LIB.lib_type == 'S':
                cmd_align += f' {LIB.trimmed_file_1}'
        
            cmd_align += f' 1> {libD}/logs/log.Bismark.txt'
            CMD_list.append( [cmd_align, f'{libD}/logs/log.Bismark.txt',log_proc ] )
        
        # 2) bs2
        elif args.program == 'bs2':
            cmd_align = f'python2 {prog.bs2_align}'
            
            if LIB.lib_type == 'P':
                cmd_align += f' -1 {LIB.trimmed_file_1} -2 {LIB.trimmed_file_2}'
            elif LIB.lib_type == 'S':
                cmd_align += f'-i {LIB.trimmed_file_1}'

            cmd_align += f' -m 0 --aligner=bowtie -g {REF.fastaF} --db {refD} --temp_dir={tmpD} -o {bamF}'
            CMD_list.append( [cmd_align, f'{libD}/logs/log.bs2_align.txt',log_proc ] )
    

    with Pool(int(job_number)) as pool:
        multi_returns = pool.map(multi_run_wrapper, CMD_list)
        if not sum(multi_returns) == 0:
            printlog(log_proc, f'READS MAPPING IS NOT COMPLETED.')
            return False
        else: 
            printlog(log_proc, '\n\t.. done\n')
            for i in range(len(LIBs)):
                LIB = LIBs[i]
                cur_lib = f'LIB{i+1}'
                libD = callD +'/'+ LIB.lib_name
                dataD = libD + '/data'
                
                sorted_bamF = f'{dataD}/{LIB.lib_name}_{args.program}.sorted.bam'

                if args.program == 'bismark':
                    sub.call(f'mv {dataD}/*.bam {sorted_bamF}', shell=True)
                    sub.call(f'mv {dataD}/*bismark*_report.txt {dataD}/{LIB.lib_name}_bismark_mapping_report.txt', shell=True)
                
                multi_returns = pool.map(multi_run_wrapper, sortCMD_list)
                LIB.bam_file = sorted_bamF
                
                ##2.5 multiQC
                cmd = f'{prog.multiqc} {dataD} -o {libD}'
                call_return = subproc(cmd, 'stdout',log_proc)
            return True
    

# ----------------------------------------------------------------------------

def Bismark_calling(callD):
    global LIBs
    global args
    global prog
    global SAMPLE
    
    ## set job_number
    job_number, core = set_job_number(len(LIBs), args.core)
    callingCMD_list = []
    bedgraphCMD_list = []
    gzipCMD_list = []
    splitCMD_list = []
    
    call_complete = True
    
    for i in range(len(LIBs)):
        LIB = LIBs[i]
        cur_lib = f'LIB{i+1}'
        
        ## result space
        libD = callD +'/'+ LIB.lib_name
        dataD = libD + '/data'
        methylD = libD + '/methylcontext'
    
        chr_CpGD = methylD + '/CpG_chr'
        chr_CHGD = methylD + '/CHG_chr'
        chr_CHHD = methylD + '/CHH_chr'

        ##3.calling:  Methylation extractor
        if not args.skip_calling:
            cmd = f'{prog.bismark_methylation_extractor} -{LIB.lib_type.lower()} --no_overlap --comprehensive --gzip --CX --cytosine_report --genome_folder {refD}/Bisulfite_Genome -o {methylD} {LIB.bam_file}'
            callingCMD_list.append( [cmd, f'{libD}/logs/log.Bismark_call.txt',log_proc ] )
        
        ##4.bedGraph
            cpgcmd = f'{prog.bismark_bedgraph} -o CpG.cov --dir {methylD} {methylD}/CpG_context_*'
            bedgraphCMD_list.append( [cpgcmd, f'{libD}/logs/log.CpG_bismark2bedGraph.txt',log_proc ] )
        
            chgcmd = f'{prog.bismark_bedgraph} -o CHG.cov --dir {methylD} --CX {methylD}/CHG_context_*'
            bedgraphCMD_list.append( [chgcmd, f'{libD}/logs/log.CHG_bismark2bedGraph.txt',log_proc ] )
        
            chhcmd = f'{prog.bismark_bedgraph} -o CHH.cov --dir {methylD} --CX {methylD}/CHH_context_*'
            bedgraphCMD_list.append( [chhcmd, f'{libD}/logs/log.CHH_bismark2bedGraph.txt',log_proc ] )
    
            ## split files
            gzipCMD_list.append( [f'gzip -dc {methylD}/CpG.cov.gz.bismark.cov.gz 1> {libD}/{LIB.lib_name}_CpG.cov.txt', 'stdout', log_proc])
            gzipCMD_list.append( [f'gzip -dc {methylD}/CHG.cov.gz.bismark.cov.gz 1> {libD}/{LIB.lib_name}_CHG.cov.txt', 'stdout', log_proc])
            gzipCMD_list.append( [f'gzip -dc {methylD}/CHH.cov.gz.bismark.cov.gz 1> {libD}/{LIB.lib_name}_CHH.cov.txt', 'stdout', log_proc])
    
            splitCMD_list.append( [f'{prog.split} {libD}/{LIB.lib_name}_CpG.cov.txt {chr_CpGD}/{LIB.lib_name}', 'stdout', log_proc] )
            splitCMD_list.append( [f'{prog.split} {libD}/{LIB.lib_name}_CHG.cov.txt {chr_CHGD}/{LIB.lib_name}', 'stdout', log_proc] )
            splitCMD_list.append( [f'{prog.split} {libD}/{LIB.lib_name}_CHH.cov.txt {chr_CHHD}/{LIB.lib_name}', 'stdout', log_proc] )
        
    
    ## CALLING
    if args.skip_calling :
        printlog(log_proc, '3. Calling (bismark) .. skip')
    else: 
        printlog(log_proc, '3. Calling (bismark) .. start')
    
        with Pool(int(job_number)) as pool:
            multi_returns = pool.map(multi_run_wrapper, callingCMD_list)
            if not sum(multi_returns) == 0:
                printlog(log_proc, f'METHYLATION CALLING IS NOT COMPLETED.')
                call_complete =  False
                return False
            else:
                printlog(log_proc, '\n\tMethylation extractor .. done')
        
        ## bedgreph
        job_number, core = set_job_number(len(LIBs)*3, args.core)
        printlog(log_proc,'4. bedGraph .. start')
        
        with Pool(int(job_number)) as pool:
            multi_returns = pool.map(multi_run_wrapper, bedgraphCMD_list)
            if not sum(multi_returns) == 0:
                printlog(log_proc, f'BEDGRAPH TRANSFORMING IS NOT COMPLETED.')
                call_complete =  False
            else:
                printlog(log_proc,'\n\tbedGraph .. done')
    

        printlog(log_proc,'4-2. bedGraph file processing .. start')
        processing_returns = 0
        with Pool(int(job_number)) as pool:
            multi_returns = pool.map(multi_run_wrapper, gzipCMD_list)
            processing_returns += sum(multi_returns)
    
            multi_returns = pool.map(multi_run_wrapper, splitCMD_list)
            processing_returns += sum(multi_returns)
    
            if not processing_returns == 0:
                printlog(log_proc, f'BEDGRAPH PROCESSING IS NOT COMPLETED.')
                call_complete =  False
            else:
                printlog(log_proc,'\n\tbedGraph file processing .. done')

    return call_complete
        

# ----------------------------------------------------------------------------
def BS2_calling(callD):
    global LIBs
    global args
    global prog
    global SAMPLE
    ## set job_number
    job_number, core = set_job_number(len(LIBs), args.core)
    callingCMD_list = []
    convertCMD_list = []
    splitCMD_list = []
    call_complete =  True
    
    for i in range(len(LIBs)):
        LIB = LIBs[i]
        cur_lib = f'LIB{i+1}'
        
        ## result space
        libD = callD +'/'+ LIB.lib_name
        dataD = libD + '/data'
        methylD = libD + '/methylcontext'

        ##3.calling:  BSseeker2 bs_call_methylation
        if not args.skip_calling:
            cmd = f'python2 {prog.bs2_call} -i {LIB.bam_file} -o {methylD}/bs2_call --sorted --rm-overlap -d {REF.refD}/{REF.ucsc_name}.fa_bowtie 1> log.bs2_call.txt'
            callingCMD_list.append( [cmd, f'{libD}/logs/log.bs2_call.txt',log_proc ] )
            convertCMD_list.append( [f'{prog.bs22bismark} {methylD}/bs2_call.CGmap {libD}/{LIB.lib_name}','stdout',log_proc ])
            
            splitCMD_list.append( [f'{prog.split} {libD}/{LIB.lib_name}_CpG.cov.txt {methylD}/CpG_chr/{LIB.lib_name}', 'stdout', log_proc] )
            splitCMD_list.append( [f'{prog.split} {libD}/{LIB.lib_name}_CHG.cov.txt {methylD}/CHG_chr/{LIB.lib_name}', 'stdout', log_proc] )
            splitCMD_list.append( [f'{prog.split} {libD}/{LIB.lib_name}_CHH.cov.txt {methylD}/CHH_chr/{LIB.lib_name}', 'stdout', log_proc] )
    
    
    ## CALLING
    if args.skip_calling :
        printlog(log_proc, '3. Calling (bs2) .. skip')
    else: 
        printlog(log_proc, '3. Calling (bs2) .. start')
    
        with Pool(int(job_number)) as pool:
            multi_returns = pool.map(multi_run_wrapper, callingCMD_list)
            if not sum(multi_returns) == 0:
                printlog(log_proc, f'METHYLATION CALLING IS NOT COMPLETED.')
                call_complete =  False
                return False
            else:
                printlog(log_proc, '\n\tBSseeker2 bs2_call_methylation .. done')

        ## -------
        sub.call(f'gzip -d {callD}/*/methylcontext/*gz',shell=True)

        ## convert bs2 result to bismark format
        with Pool(int(job_number)) as pool:
            multi_returns = pool.map(multi_run_wrapper, convertCMD_list)
            multi_returns = pool.map(multi_run_wrapper, splitCMD_list)
            if not sum(multi_returns) == 0:
                printlog(log_proc, f'bs2 result -> bismark result IS NOT COMPLETED.')
                call_complete =  False
            else:
                printlog(log_proc, '\n\tBSseeker2 result converting .. done')
    
    return call_complete
        


def Methylation_calling(callD):
    global args
    global SAMPLE
    global LIB
    

    SAMPLE = {}
    ## sample name and output file
    for LIB in LIBs:
        name = LIB.sample_name
        libD = callD +'/'+ LIB.lib_name
        if not name in SAMPLE:
            SAMPLE[name] = {}
        SAMPLE[name][LIB.lib_name] = [ f'{libD}/{LIB.lib_name}_CpG.cov.txt', f'{libD}/{LIB.lib_name}_CHG.cov.txt', f'{libD}/{LIB.lib_name}_CHH.cov.txt' ]
    
    for i in range(len(LIBs)):
        LIB = LIBs[i]
        cur_lib = f'LIB{i+1}'
        
        ## result space
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
        

    ## calling
    if args.program == 'bismark':
        value = Bismark_calling(callD)
    
    if args.program == 'bs2':
        value = BS2_calling(callD)
    
    return value

# ----------------------------------------------------------------------------

def GMA(outD, refD):
    global REF
    global LIBs
    global args
    global prog
    global SAMPLE

    ## sample_level result
    resultD = outD + '/Analysis'
    result_logD = resultD + '/logs'
    sub.call(f'mkdir -p {result_logD}', shell=True)

    printlog(log_proc,'\n\nGene & Methyl Analysis ...'+'-'*60)
    
    if args.skip_GMA :
        printlog(log_proc, 'GMA .. skip')
        return 0
    else: 
        printlog(log_proc, 'GMA .. start')
    
        ##make parameter file
        param_file = f'{resultD}/params.txt'
        sample_names = list(SAMPLE.keys())

        p_line = f'@REF_NAME\n{REF.ucsc_name}\n\n@REF_FA\n{REF.fastaF}\n\n@GENE_GTF\n{REF.annotF}\n\n'
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
        gma_return = subproc(f'{prog.gene} -p {param_file} -o {resultD} -cpu {args.core}', f'{resultD}/log.genemethylAnal.txt',log_proc)
        return gma_return


def VISandDMR(outD):
    global REF
    global LIBs
    global args
    global prog
    global SAMPLE

    ## VISUALIZATION ------------------------------------------------------------------------
    printlog(log_proc,'\n\nVisualization ... '+'-'*60)
    resultD = outD + '/Analysis'
    result_logD = resultD + '/logs'
    sample_names = list(SAMPLE.keys())

    ## make visualization parameter file
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
    wig_cmd =  f'{prog.Rscript} {prog.analCpG_vis_paral} {vis_param_file} {resultD}/Ideogram.bed {resultD}/ {args.core}'
    subproc(wig_cmd, result_logD + f'/log.bed2wigANDwindow', log_proc)
    
    visCMD_list = []
    cmd_number = 0
    for sample in sample_names:
        sampleD = resultD + f'/{sample}'

        # wig2bigwig
        bw_cmd = f'wigToBigWig {sampleD}/CpG_methylLev.wig {resultD}/ref.size {sampleD}/CpG_methylLev.bw'
        visCMD_list.append( [bw_cmd, result_logD + f'/log.{sample}.wig2bw', log_proc] )

        #genomic context visualization
        context_cmd = f'{prog.Rscript} {prog.contextLev} {sampleD}/methylation.Genomic_Context.CpG.txt {sample} {sampleD}/Genomic_Context_CpG.pdf 1> {sampleD}/Avg_Genomic_Context_CpG.txt'    
        visCMD_list.append( [context_cmd, result_logD+f'/log.{sample}.Genomic_Context', log_proc] )

        #circos
        circos_cmd = f'{prog.circos} {resultD}/Ideogram.bed {sampleD}/ {sampleD}/Circos_{sample}.CpG_UMRs_LMRs.pdf'
        visCMD_list.append( [circos_cmd, result_logD+f'/log.{sample}.Circos', log_proc] )
        cmd_number += 3
    
    job_number, core = set_job_number(cmd_number, args.core)
    with Pool(int(job_number)) as pool:
        multi_returns = pool.map(multi_run_wrapper, visCMD_list)
        if not sum(multi_returns) == 0:
            printlog(log_proc, f'UPPER STEPS ARE NOT COMPLETED.')
    
    
    ## get Promoter
    promoterF = f'{resultD}/annotations/promoter.bed'


    ## DMR analysis   ------------------------------------------------------------------------------
    printlog(log_proc,'\n\nDMR Analysis ...'+ '-'*60)
    if args.skip_DMR :
        printlog(log_proc,'DMR .. skip\n')
        return 0

    elif DMR.param_values() == [] or DMR.param_values() == [""]:
        printlog(log_proc, "NO DMR analysis group set\n")
        return 1

    else:
        printlog(log_proc,'DMR .. start')
        dmrD = resultD + '/DMR'
        sub.call(f'mkdir -p {dmrD}', shell=True)
        

    #----------------------------
        complete = 0

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

            analD = f'{dmrD}/{control}.{case}'
            anal_logD = f'{dmrD}/{control}.{case}/logs'
            sub.call(f'mkdir -p {anal_logD}', shell=True)
           
            ## BSmooth Rscript        ====================================
            if args.bsmooth:
                ## make bsmooth_param
                bsmooth_paramF = f'{analD}/bsmooth_param.txt'
                bsmooth_line = ''
                for group in [control, case]:
                    for lib in SAMPLE[group].keys():
                        bsmooth_line += f'{group}\t{lib}\t{SAMPLE[group][lib][0]}\n'
                
                open(bsmooth_paramF, 'w').write(bsmooth_line.rstrip())

                ## BSmooth
                bsmooth_cmd = f'{prog.Rscript} {prog.bsmooth} {analD}/bsmooth_param.txt {control},{case} {resultD} {args.qvalue} 1> {anal_logD}/log.DMR_BSmooth.txt'
                bsmooth_return = subproc(bsmooth_cmd, f'{anal_logD}/log.DMR_BSmooth.txt',log_proc)
            

                methylLevelF = f"{analD}/DMR_q{args.qvalue}.bed"
            
            # convert 
                reformF = f'{analD}/reform.DMR_q{args.qvalue}.bed'
                cmd_bsreform = f'{prog.bsmoothreform} {methylLevelF} 1> {reformF}' 
                subproc(cmd_bsreform, "stdout",log_proc)

            ##get methylation Levels on promoter region
                intersect_promoterF = f'{analD}/intersection.DMR2Promoter.txt'
                cmd_intersect = f'bedtools intersect -wa -wb -a {promoterF} -b {reformF} > {intersect_promoterF}'
                subproc(cmd_intersect, "stdout",log_proc)


            ## gprofiler
                dmc_genelistF = f'{analD}/DMR_genelist.txt'
                cmd_genename= f'cut -f 4 {intersect_promoterF} | cut -d \';\' -f1| cut -d \':\' -f2 |sort -u > {dmc_genelistF}'
                subproc(cmd_genename, "stdout",log_proc)

                subproc(f'{prog.Rscript} {prog.gprofiler} {dmc_genelistF} {REF.ucsc_name} {analD}/DMR_gene', f'{anal_logD}/log.gprofiler.txt',log_proc)
            
                return bsmooth_return

            else:
                methylkitD = f'{analD}/methylkit'
                sub.call(f'mkdir -p {methylkitD}', shell=True)

                ## methylKit Rscript        ====================================
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
                subproc(f'{prog.Rscript} {methylkit_script}', f'{anal_logD}/log.methylkit.txt',log_proc)

            ## Methylation Level
                methylLevelF = f"{analD}/DMC_q{args.qvalue}.bed"
                cmd = f'{prog.getDMC} {methF} {rawF} {args.qvalue} {treatments} {methylLevelF} {control} {case}'
                subproc(cmd, f'{anal_logD}/log.DMC.txt',log_proc)
        
            ##get methylation Levels on promoter region
                intersect_promoterF = f'{analD}/intersection.DMC2Promoter.txt'
                cmd_intersect = f'bedtools intersect -wa -wb -a {promoterF} -b {analD}/reform.DMC_q{args.qvalue}.bed > {intersect_promoterF}'
                subproc(cmd_intersect, "stdout",log_proc)


            ## gprofiler
                dmc_genelistF = f'{analD}/DMC_genelist.txt'
                cmd_genename= f'cut -f 4 {intersect_promoterF} | cut -d \';\' -f1| cut -d \':\' -f2 |sort -u > {dmc_genelistF}'
                subproc(cmd_genename, "stdout",log_proc)

                subproc(f'{prog.Rscript} {prog.gprofiler} {dmc_genelistF} {REF.ucsc_name} {analD}/DMC_gene', f'{anal_logD}/log.gprofiler.txt',log_proc)

            ## get stat
                statF = f'{analD}/Stat.CpG_qvalue{args.qvalue}.out'
                cmd_stat = f'{prog.getStat} {intersect_promoterF} {analD}'
                subproc(cmd_stat, f"{anal_logD}/log.getStat.txt",log_proc)
        
            ## DMR
                if not os.path.isfile(methylLevelF):
                    printlog(log_proc,"There is no significantly differentially methylated base.")
                    printlog(log_proc, "->DMR analysis is skipped.\n")
                    return 0
                else:
                    ## getDMR
                    completed = subproc(f'{prog.dmr} {rawF} {args.qvalue} 1> {analD}/DMR_{args.qvalue}.bed', f'{anal_logD}/log.getDMR.txt',log_proc)
                    return completed
            


if __name__ == "__main__":

    parser = make_argparser()
    args = parser.parse_args()
    paramF = args.param
    outD = os.path.abspath(args.out)
    
    
    make_result_space(outD)
    get_parameters(paramF)
    refD = pipe_path + f'/reference/{REF.ucsc_name}'
    REF.refD = refD
    
    fastaF, annoF = check_inputs(refD) ## check input files
    
    callD = args.calling_data.rstrip("/")
    prog = Programs(pipe_path)

    REF_bisulfite_genome(refD)
    
    completed = Reads_preprocessing(callD)
    if not completed:
        printlog(log_proc,"ERR: problem in preprocessing step")
        sys.exit()
    
    completed = Bismark_mapping(callD, refD)
    if not completed:
        printlog(log_proc,"ERR: problem in mapping step")
        sys.exit()

    completed = Methylation_calling(callD)
    if not completed:
        printlog(log_proc,"ERR: problem in methylation calling step")
        sys.exit()
    else:
        printlog(log_proc, f'\nPre-processing and calling is done')

        if not args.skip_calling :
            ## multiQC for all data
            cmd = f'{prog.multiqc} {callD} -o {callD}'
            call_return = subproc(cmd, 'stdout',log_proc)

    completed = GMA(outD, refD)
    if completed != 0:
        printlog(log_proc,"ERR: problem in Gene Methylation Analysis\n")
    else:
        printlog(log_proc, f'\n Gene Methylation Analysis is done\n')
    
    completed = VISandDMR(outD)
    if completed != 0:
        printlog(log_proc,"ERR: problem in Visualization and Differentail methylation analysis\n")

    printlog(log_proc, f'\n\n ALL PIPELINE IS COMPLETED \n\n')


# ====================================================================================================================================================

