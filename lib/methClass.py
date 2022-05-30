import sys

class Parameters():
    def __init__(self, config, key1):
        for key2 in config[key1]:
            self.__dict__[key2] = config[key1][key2]

    def param_keys(self):
        return list(self.__dict__.keys())

    def param_values(self):
        return list(self.__dict__.values())

class Programs():
    def __init__(self, pipe_path):
        binD = pipe_path + "/bin"
        self.multiqc = 'multiqc'

       # Requirment
        self.trim_galore = 'trim_galore'
        self.samtools = 'samtools'
<<<<<<< HEAD

        bismarkD = f'{binD}/Bismark'
        self.bismark = f'{bismarkD}/bismark'
        self.bismark_indexing = f'{bismarkD}/bismark_genome_preparation'
        self.bismark_methylation_extractor = f'{bismarkD}/bismark_methylation_extractor'
        self.bismark_bedgraph = f'{bismarkD}/bismark2bedGraph'
=======
        self.bismark = 'bismark'
        self.bismark_indexing = 'bismark_genome_preparation'
        self.bismark_methylation_extractor = 'bismark_methylation_extractor'
        self.bismark_bedgraph = 'bismark2bedGraph'
>>>>>>> 8bd99478acddf4ab6f72ca7b810f24915461cd8a
        
        bs2D = f'{binD}/BSseeker2'
        self.bs2_build = f'{bs2D}/bs_seeker2-build.py'
        self.bs2_align = f'{bs2D}/bs_seeker2-align.py'
        self.bs2_call = f'{bs2D}/bs_seeker2-call_methylation.py'
        self.Rscript = 'Rscript'
        
        self.cutadapt = 'cutadapt'

        scriptD = binD + '/script'
        self.split = scriptD+'/splitF_bychr.py'
        self.getDMC = scriptD + '/getDMC.py'
        self.getStat = scriptD + '/getStat.pl'
        self.dmr = scriptD + '/getDMR_500bp.pl'
        self.bsmooth = scriptD + '/BSmooth_DMR.R'
        self.bsmoothreform = scriptD + '/BSmooth_reform.py'
        self.bs22bismark = scriptD + '/bs22bismark.py'

        self.gene = binD + '/GMA/Gene-Methyl-Analysis.pl'
        self.gprofiler = binD + '/GMA/gProfiler.R'

        visD = f"{pipe_path}/bin/vis_script/"
        self.analCpG_vis_paral = visD + 'visualization_parallel.R'
        self.window = visD + "win100kb_methylLevel.R"
        self.circos = visD + 'GMA.Circos_100kb.R'
        self.contextLev = visD + 'genomic_context_levels.R'


