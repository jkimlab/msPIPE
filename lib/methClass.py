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

        visD = f"{pipe_path}/bin/vis_script/"
        self.analCpG = visD + 'visualization_parallel.R'
        self.window = visD + "win100kb_methylLevel.R"
        self.circos = visD + 'GMA.Circos_100kb.R'
        self.contextLev = visD + 'genomic_context_levels.R'


