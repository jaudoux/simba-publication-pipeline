import os 
from snakemake.utils import R

__author__ = "Jérôme Audoux (jerome.audoux@inserm.fr)"

OLD_DATASETS    = ["GRCh38-100bp-160M-normal", "GRCh38-100bp-160M-somatic", "GRCh38-200bp-160M-somatic"]
DATASETS    = ["GRCh38-101bp-160M-normal", "GRCh38-101bp-160M-somatic","GRCh38-150bp-160M-somatic", "GRCh38-150bp-160M-normal"]

# DIRECTORIES
ABS_DIR       = os.getcwd()
DOWNLOAD_DIR  = "download"
DATASET_DIR   = "dataset"
MAPPING_DIR   = "mapping"
CALLING_DIR   = "calling"
FIGURES_DIR   = "figures"
BENCHCT_DIR   = "benchCT"
SCRIPTS_DIR   = "scripts"
BENCHCT_FUSION_DIR    = BENCHCT_DIR + "/fusions"
BENCHCT_MUTATION_DIR  = BENCHCT_DIR + "/mutations"
BENCHCT_MAPPING_DIR  = BENCHCT_DIR + "/mapping"
TMP_DIR       = "/data/scratch/audoux"

# PARAMETERS
MAX_SPLICE_LENGTH = 300000
MIN_CHIM_OVERHANG = 15
NB_THREADS_MAPPERS = 10

# FILES
GENOME_DIR         = "/data/genomes/GRCh38/chr"
REFERENCE_BASENAME = "/data/genomes/GRCh38/GRCh38_with_MT"
ANNOTATIONS        = DOWNLOAD_DIR + "/Homo_sapiens.GRCh38.86.gtf.gz"
POLYMORPHISMS      = DOWNLOAD_DIR + "/20M-human-SNPs.vcf.gz"
FLUX_ERROR_MODEL   = ABS_DIR + "/" + DATASET_DIR  + "/illumina-hiseq2500-error-model.err"

# DATASET CONDITIONS
CONDITIONS = {}
CONDITIONS['normal']  = { 'sub_rate' : '0.0014', 'indel_rate' : '0.0001', 'nb_fusions' : '0', 'k' : '-0.7', 'x0' : '15000', 'x1' : '225000000', 'nb_molecules' : 10000000 }
CONDITIONS['somatic']   = { 'sub_rate' : '0.0016', 'indel_rate' : '0.0001', 'nb_fusions' : '100', 'k' : '-0.7', 'x0' : '25000', 'x1' : '625000000', 'nb_molecules' : 10000000 }

CHIMSEGMENT_VALUES = [10,12,14,16,18,20,22,24,26,28,30]

# MAPPERS
STAR_NAME          = "star-2.5.2b"
HISAT2_NAME        = "hisat2-2.0.4"
HISAT2_2PASS_NAME  = HISAT2_NAME + "_2pass"
CRAC_NAME          = "crac-2.5.0"

# CALLERS
GATK_NAME       = "gatk"
FREEBAYES_NAME  = "freebayes"
SAMTOOLS_NAME   = "samtools"
CRACTOOLS_NAME  = "cractools"

# MAPPERS DIRECTORIES
STAR_DIR          = MAPPING_DIR + "/" + STAR_NAME
CRAC_DIR          = MAPPING_DIR + "/" + CRAC_NAME
HISAT2_DIR        = MAPPING_DIR + "/" + HISAT2_NAME
HISAT2_2PASS_DIR  = MAPPING_DIR + "/" + HISAT2_2PASS_NAME

MAPPERS = [STAR_NAME, HISAT2_NAME, HISAT2_2PASS_NAME]
CALLERS = [GATK_NAME, FREEBAYES_NAME, SAMTOOLS_NAME]

MAPPER_CALLER_PIPELINES = expand("{mapper}_{caller}", 
                        mapper = MAPPERS,
                        caller = CALLERS)
MAPPER_CALLER_PIPELINES.append(CRAC_NAME + "_" + STAR_NAME + "_" + GATK_NAME)
MAPPER_CALLER_PIPELINES.append(HISAT2_2PASS_NAME + "_" + GATK_NAME + "_" + STAR_NAME + "_" + GATK_NAME)
CALLING_PIPELINES = list(MAPPER_CALLER_PIPELINES)
CALLING_PIPELINES.append(CRAC_NAME)

MAPPERS.append(CRAC_NAME)

MUTATION_TYPES = ["snp","insertion","deletion"]

# BINARIES
SIMCT                   = "simCT"
STAR                    = "STAR"
CRAC                    = "crac"
MARKDUPLICATES          = "PicardCommandLine MarkDuplicates"
ADDORREPLACEREADGROUPS  = "PicardCommandLine AddOrReplaceReadGroups"
CREATESEQUENCEDICT      = "PicardCommandLine CreateSequenceDictionary"
GATK                    = "java -Djava.io.tmpdir=" + TMP_DIR + " -jar /data/share/GATK/GenomeAnalysisTK.jar"
SPLITNCIGARREADS        = GATK + " -T SplitNCigarReads"
HAPLOTYPECALLER         = GATK + " -T HaplotypeCaller"
HISAT2                  = "hisat2"
FREEBAYES               = "freebayes"
MPILEUP                 = "samtools mpileup"
BCFTOOLSCALL            = "bcftools call"
CRACTOOLSEXTRACT        = "cractools extract"
REMOVESIMCTFNCHIM       = SCRIPTS_DIR + "/removeSimCTfNchimeras.pl"

ruleorder: benchct_mutations > benchct_generic > vcf_sort > cractools > HaplotypeCaller > SplitNCigarReads > MarkDuplicates > addRG > mapping_stats > hisat2_2pass > hisat2 > star > crac

rule all:
  input: 
    figures_fusion = expand("{dir}/{sample}/{type}_accuracy_sensitivity.pdf",
                     dir = FIGURES_DIR,
                     #sample = ["GRCh38-101bp-160M-somatic","GRCh38-150bp-160M-somatic"],
                     sample = DATASETS,
                     type = ["fusions"]),
    figures = expand("{dir}/{sample}/{type}_accuracy_sensitivity.pdf",
                     dir = FIGURES_DIR,
                     sample = DATASETS,
                     type = ["mutations", "mapping"]),
    tp_figs = expand("{dir}/{sample}/{event}-true-positives.pdf",
                     dir = FIGURES_DIR,
                     sample = DATASETS,
                     event = MUTATION_TYPES),
    mapping_time_fig = FIGURES_DIR + "/mapping-time.pdf",

rule flux_par:
  input:
    error_model = FLUX_ERROR_MODEL
  output: DATASET_DIR + "/{params}-{condition}/flux.par"
  params:
    k   = lambda wildcards: CONDITIONS[wildcards.condition]["k"],
    x0  = lambda wildcards: CONDITIONS[wildcards.condition]["x0"],
    x1  = lambda wildcards: CONDITIONS[wildcards.condition]["x1"],
    tmp_dir = TMP_DIR
  run:
    shell("echo 'EXPRESSION_K {params.k}' > {output}")
    shell("echo 'EXPRESSION_X0 {params.x0}' >> {output}")
    shell("echo 'EXPRESSION_X1 {params.x1}' >> {output}")
    shell("echo 'TMP_DIR {params.tmp_dir}' >> {output}")
    shell("echo 'ERR_FILE {input.error_model}' >> {output}")

rule simct:
  input:
    genome = GENOME_DIR,
    annot  = ANNOTATIONS,
    flux_param = DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}/flux.par"
  params:
    nb_reads   = "{nb_reads}",
    sub_rate   = lambda wildcards: CONDITIONS[wildcards.condition]["sub_rate"],
    indel_rate = lambda wildcards: CONDITIONS[wildcards.condition]["indel_rate"],
    nb_fusions = lambda wildcards: CONDITIONS[wildcards.condition]["nb_fusions"],
    nb_molecules = lambda wildcards: CONDITIONS[wildcards.condition]["nb_molecules"], 
    out_dir = DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}",
  output:
    info = DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}/info.txt",
    chimeras = DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}/chimeras.tsv.gz",
    r1 =  DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}/reads_1.fastq.gz",
    r2 =  DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}/reads_2.fastq.gz",
  log: DATASET_DIR + "/{genome}-{read_length}bp-{nb_reads}M-{condition}/stderr.log",
  shell: """{SIMCT} -g {input.genome} -a {input.annot} -o {params.out_dir} \
            -s {params.sub_rate} -i {params.indel_rate} \
            -d {params.indel_rate} -f {params.nb_fusions} --nb-reads {params.nb_reads}000000 \
            --nb-molecules {params.nb_molecules} \
            --fragment-length 250 --fragment-sd 50 --uniq-ids --vcf-file {POLYMORPHISMS} \
            --vcf-ratio 0.95 --flux-par {input.flux_param} 2> {log}"""

rule bam_index:
  input: "{file}.bam"
  output: "{file}.bam.bai"
  shell: "samtools index {input}"

rule fa_dict:
  input: "{file}.fa"
  output: "{file}.dict"
  shell: "{CREATESEQUENCEDICT} R={input} O={output}"

rule vcf_sort:
  input: "{file}.vcf"
  output: temp("{file}_sorted.vcf")
  shell: "grep '^#' {input} > {output} && grep -v '^#' {input} | sort -t$'\t' -k1,1 -k2,2n >> {output}"

rule vcf_compress:
  input: "{file}_sorted.vcf"
  output: "{file}.vcf.gz"
  shell: "bgzip -c {input} > {output}"

rule vcf_index:
  input: "{file}.vcf.gz"
  output: "{file}.vcf.gz.csi"
  shell: "bcftools index {input}"

rule addRG:
  input: MAPPING_DIR + "/{mapper}/{sample}.bam"
  output: temp(MAPPING_DIR + "/{mapper}/{sample}_with_RG.bam")
  shell: "{ADDORREPLACEREADGROUPS} I={input} O={output} SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample"

rule MarkDuplicates:
  input: MAPPING_DIR + "/{mapper}/{sample}_with_RG.bam"
  output: 
    bam = temp(MAPPING_DIR + "/{mapper}/{sample}_marked_duplicates.bam"),
    log = MAPPING_DIR + "/{mapper}/{sample}_marked_dup_metrics.txt"
  shell: "{MARKDUPLICATES} I={input} O={output.bam} M={output.log}"

rule SplitNCigarReads:
  input: 
    bam = MAPPING_DIR + "/{mapper}/{sample}_marked_duplicates.bam",
    bai = MAPPING_DIR + "/{mapper}/{sample}_marked_duplicates.bam.bai",
    ref = REFERENCE_BASENAME + ".fa",
    dic = REFERENCE_BASENAME + ".dict"
  output: MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam"
  shell: """{SPLITNCIGARREADS} -R {input.ref} -I {input.bam} -o {output} \
            -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
            -U ALLOW_N_CIGAR_READS"""

rule HaplotypeCaller:
  input:
    bam = MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam",
    bai = MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam.bai",
    ref = REFERENCE_BASENAME + ".fa",
    dic = REFERENCE_BASENAME + ".dict"
  output: CALLING_DIR + "/{mapper}_" + GATK_NAME + "/{sample}.vcf"
  shell: "{HAPLOTYPECALLER} -R {input.ref} -I {input.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o {output}"

rule star:
  input:
    r1 = DATASET_DIR + "/{sample}/reads_1.fastq.gz",
    r2 = DATASET_DIR + "/{sample}/reads_2.fastq.gz",
  params:
    index = "/data/indexes/STAR/GRCh38_with_MT"
  output:
    dir   = STAR_DIR + "/{sample}/",
    bam   = STAR_DIR + "/{sample}.bam",
    time  = STAR_DIR + "/{sample}-time.txt"
  log: STAR_DIR + "/{sample}-star.log"
  threads: NB_THREADS_MAPPERS
  run:
    if not os.path.exists(output.dir):
      os.makedirs(output.dir)
    shell("""/usr/bin/time -v -o {output.time} \
            {STAR} --genomeDir {params.index} --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat --twopassMode Basic --outStd SAM \
            --outFileNamePrefix {output.dir} --alignMatesGapMax {MAX_SPLICE_LENGTH} \
            --alignIntronMax {MAX_SPLICE_LENGTH} --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 --chimSegmentReadGapMax parameter 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --runThreadN {threads} 2> {log} | samtools view -@5 -bS - \
            2> {log} | samtools sort -@5 -m 5G - -o {output.bam}""")

rule crac:
  input:
    r1 = DATASET_DIR + "/{sample}/reads_1.fastq.gz",
    r2 = DATASET_DIR + "/{sample}/reads_2.fastq.gz",
  params:
    index = "/data/indexes/crac/GRCh38_with_MT"
  output:
    bam   = CRAC_DIR + "/{sample}.bam",
    time  = CRAC_DIR + "/{sample}-time.txt"
  log: CRAC_DIR + "/{sample}-crac.log"
  threads: NB_THREADS_MAPPERS
  shell: """/usr/bin/time -v -o {output.time} \
            {CRAC} -k 22 -o - --detailed-sam --no-ambiguity --deep-snv \
            -r {input.r1} {input.r2} --nb-threads {threads} -i {params.index} \
            2> {log} | samtools view -@5 -bS - | samtools sort -@5 -m 5G - -o {output.bam}"""
    

rule hisat2:
  input: 
    r1 = DATASET_DIR + "/{sample}/reads_1.fastq.gz",
    r2 = DATASET_DIR + "/{sample}/reads_2.fastq.gz",
  params:
    index = "/data/indexes/hisat2/GRCh38_with_MT/GRCh38_with_MT"
  output:
    bam           = HISAT2_DIR + "/{sample}.bam",
    novel_splice  = HISAT2_DIR + "/{sample}_novel_splice.bed",
    time          = HISAT2_DIR + "/{sample}-time.txt"
  log: HISAT2_DIR + "/{sample}-hisat.log"
  threads: NB_THREADS_MAPPERS
  shell: """/usr/bin/time -v -o {output.time} \
            {HISAT2} -x {params.index} -1 {input.r1} -2 {input.r2} \
            --max-intronlen {MAX_SPLICE_LENGTH} \
            --novel-splicesite-outfile {output.novel_splice} \
            -p {threads} 2> {log} \
            | samtools view -bS - | samtools sort - -o {output.bam}"""

rule hisat2_2pass:
  input: 
    r1 = DATASET_DIR + "/{sample}/reads_1.fastq.gz",
    r2 = DATASET_DIR + "/{sample}/reads_2.fastq.gz",
    novel_splice = HISAT2_DIR + "/{sample}_novel_splice.bed",
  params:
    index = "/data/indexes/hisat2/GRCh38_with_MT/GRCh38_with_MT"
  output:
    bam           = HISAT2_2PASS_DIR + "/{sample}.bam",
    novel_splice  = HISAT2_2PASS_DIR + "/{sample}_novel_splice.bed",
    time          = HISAT2_2PASS_DIR + "/{sample}-time.txt"
  log: HISAT2_DIR + "/{sample}-hisat.log"
  threads: NB_THREADS_MAPPERS
  shell: """/usr/bin/time -v -o {output.time} \
            {HISAT2} -x {params.index} -1 {input.r1} -2 {input.r2} \
            --max-intronlen {MAX_SPLICE_LENGTH} \
            --novel-splicesite-infile {input.novel_splice} \
            --novel-splicesite-outfile {output.novel_splice} \
            -p {threads} 2> {log} \
            | samtools view -bS - | samtools sort - -o {output.bam}"""

rule mapping_stats:
  input: expand("{dir}/{mapper}/{sample}-time.txt", dir = MAPPING_DIR, mapper = MAPPERS, sample = DATASETS)
  output: MAPPING_DIR + "/time.tsv"
  run:
    fo = open(output[0],'w')
    for i,software in enumerate(MAPPERS):
      for j,sample in enumerate(DATASETS):
        f = open(MAPPING_DIR + "/" + software + "/" + sample + "-time.txt")
        for line in f:
          fields = line.rstrip().split(":")
          value = fields[1].lstrip()
          stat_name = fields[0].lstrip()
          stat_short_name = ''
          if stat_name == 'User time (seconds)':
            stat_short_name = "time"
            value = str(int(float(value)) / 60)
          elif stat_name == 'Maximum resident set size (kbytes)':
            stat_short_name = "memory"
            value = str(int(float(value)) / 1000000)
          if stat_short_name != '':
            fo.write("\t".join([software,sample,stat_short_name,value]) + "\n")

rule mapping_stats_figure:
  input: MAPPING_DIR + "/time.tsv"
  output: FIGURES_DIR + "/mapping-time.pdf"
  version: "0.01"
  run:
    R("""
    library(ggplot2)
    dat <- read.table("{input}")
    ggplot(dat, aes(x=V1,y=V4,fill=V2)) + geom_bar(stat="identity",position="dodge") + facet_grid(V3~.,scales="free")
    ggsave("{output}",width=10,height=3)
    """)

rule freebayes:
  input:
    bam = MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam",
    bai = MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam.bai",
    ref = REFERENCE_BASENAME + ".fa"
  output:
    vcf   = CALLING_DIR + "/{mapper}_" + FREEBAYES_NAME + "/{sample}.vcf",
    time  = CALLING_DIR + "/{mapper}_" + FREEBAYES_NAME + "/{sample}-time.txt"
  shell: """/usr/bin/time -v -o {output.time} \
            {FREEBAYES} -f {input.ref} {input.bam} > {output.vcf}"""

rule samtools:
  input:
    bam = MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam",
    bai = MAPPING_DIR + "/{mapper}/{sample}_splited_N_cigar_reads.bam.bai",
    ref = REFERENCE_BASENAME + ".fa"
  output:
    vcf   = CALLING_DIR + "/{mapper}_" + SAMTOOLS_NAME + "/{sample}.vcf",
    time  = CALLING_DIR + "/{mapper}_" + SAMTOOLS_NAME + "/{sample}-time.txt"
  shell: """/usr/bin/time -v -o {output.time} \
            {MPILEUP} -ugf {input.ref} {input.bam} |
            {BCFTOOLSCALL} -vcO v > {output.vcf}"""

rule cractools:
  input:
    bam = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}.bam",
    bai = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}.bam.bai",
    ref = REFERENCE_BASENAME + ".fa"
  output:
    vcf     = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}.vcf",
    chimera = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}-chimera.tsv",
    splice  = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}-splice.bed",
    time    = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}-cractools-time.txt",
  threads: 10
  shell: """/usr/bin/time -v -o {output.time} \
            {CRACTOOLSEXTRACT} -p {threads} -s {output.splice} \
            -c {output.chimera} -m {output.vcf} --coverless-splices \
            -r {input.ref} {input.bam}"""

rule crac_star_calling_merge:
  input:
    crac_vcf = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}.vcf.gz",
    crac_vcf_index = MAPPING_DIR + "/" + CRAC_NAME + "/{sample}.vcf.gz.csi",
    star_vcf = CALLING_DIR + "/" + STAR_NAME + "_" + GATK_NAME + "/{sample}.vcf.gz",
    star_vcf_index = CALLING_DIR + "/" + STAR_NAME + "_" + GATK_NAME + "/{sample}.vcf.gz.csi",
  output:
    CALLING_DIR + "/" + CRAC_NAME + "_" + STAR_NAME + "_" + GATK_NAME + "/{sample}.vcf"
  shell: "bcftools merge {input.crac_vcf} {input.star_vcf} > {output}"

rule star_hisat_calling_merge:
  input:
    star_vcf = CALLING_DIR + "/" + STAR_NAME + "_" + GATK_NAME + "/{sample}.vcf.gz",
    star_vcf_index = CALLING_DIR + "/" + STAR_NAME + "_" + GATK_NAME + "/{sample}.vcf.gz.csi",
    hisat_vcf = CALLING_DIR + "/" + HISAT2_2PASS_NAME + "_" + GATK_NAME + "/{sample}.vcf.gz",
    hisat_vcf_index = CALLING_DIR + "/" + HISAT2_2PASS_NAME + "_" + GATK_NAME + "/{sample}.vcf.gz.csi",
  output:
    CALLING_DIR + "/" + HISAT2_2PASS_NAME + "_" + GATK_NAME + "_" + STAR_NAME + "_" + GATK_NAME + "/{sample}.vcf"
  shell: "bcftools merge --force-samples {input.hisat_vcf} {input.star_vcf} > {output}"


rule benchct_configfile_mutation:
  input: 
    infos =      DATASET_DIR + "/{sample}/info.txt",
    mutations =  DATASET_DIR + "/{sample}/mutations.vcf.gz",
    splices =    DATASET_DIR + "/{sample}/splices.bed.gz",
    chimeras =   DATASET_DIR + "/{sample}/chimeras.tsv.gz",
  output: BENCHCT_MUTATION_DIR + "/{sample}.yaml"
  version: "1.47"
  params:
    calling_pipelines = MAPPER_CALLER_PIPELINES,
    tp_dir = BENCHCT_MUTATION_DIR + "/{sample}/true-positives",
    sample = "{sample}"
  run:
    if not os.path.exists(params.tp_dir):
      os.makedirs(params.tp_dir)
    f = open(output[0], 'w')
    f.write("---\n")
    f.write("checker:\n")
    f.write("  files:\n")
    f.write("    infos: " + input.infos + "\n")
    f.write("    mutations: " + input.mutations + "\n")
    f.write("output:\n")
    f.write("  statistics: [Accuracy, Sensitivity, true-negatives, false-positives, false-negatives, true-positives, nb-elements]\n")
    f.write("softwares:\n")
    for x in params.calling_pipelines:
      f.write("  - name: " + x + "\n")
      f.write("    files:\n")
      f.write("      - name: " + CALLING_DIR + "/" + x + "/" + params.sample + ".vcf\n")
      f.write("        type: VCF\n")
      f.write("        true_positives: " + params.tp_dir + "/" + x + "\n")
      f.write("        check: all\n")
    f.write("  - name: " + CRAC_NAME + "\n")
    f.write("    files:\n")
    f.write("      - name: " + CRAC_DIR + "/" + params.sample + ".vcf\n")
    f.write("        true_positives: " + params.tp_dir + "/" + CRAC_NAME + "\n")
    f.write("        type: VCF\n")
    f.write("        check: all\n")
    f.close()

rule benchct_mutations:
  input:
    conf = BENCHCT_MUTATION_DIR + "/{sample}.yaml",
    vcf = expand("{calling_dir}/{pipeline}/{{sample}}.vcf", 
                        calling_dir = CALLING_DIR,
                        pipeline = MAPPER_CALLER_PIPELINES),
    crac_vcf = expand("{mapping_dir}/{crac_dir}/{{sample}}.vcf",
                        mapping_dir = MAPPING_DIR,
                        crac_dir    = CRAC_NAME),
  output: 
    bench = BENCHCT_MUTATION_DIR + "/{sample}.tsv",
    tp_log = expand("{bench_dir}/{{sample}}/true-positives/{pipeline}-{event}.log",
                    bench_dir = BENCHCT_MUTATION_DIR,
                    pipeline = CALLING_PIPELINES,
                    event = MUTATION_TYPES)
  shell: "benchCT -v {input.conf} > {output.bench}"

rule benchct_configfile_mapping:
  input: 
    infos =      DATASET_DIR + "/{sample}/info.txt",
    splices =    DATASET_DIR + "/{sample}/splices.bed.gz",
    bam_files = expand("{dir}/{mapper}/{{sample}}.bam",dir=MAPPING_DIR,mapper=MAPPERS),
    bai_files = expand("{dir}/{mapper}/{{sample}}.bam.bai",dir=MAPPING_DIR,mapper=MAPPERS),
    #cractools_junc = expand("{dir}/{mapper}/{{sample}}-splice.bed",dir=MAPPING_DIR,mapper=[k for k in MAPPERS if 'crac' in k]),
    cractools_junc = expand("{dir}/{mapper}/{{sample}}-splice.bed",dir=MAPPING_DIR,mapper=CRAC_NAME),
  output: BENCHCT_MAPPING_DIR + "/{sample}.yaml"
  version: "0.03"
  params:
    mapping_pipelines = MAPPERS,
    tp_dir = BENCHCT_MAPPING_DIR + "/{sample}/true-positives",
    sample = "{sample}"
  run:
    if not os.path.exists(params.tp_dir):
      os.makedirs(params.tp_dir)
    f = open(output[0], 'w')
    f.write("---\n")
    f.write("checker:\n")
    f.write("  files:\n")
    f.write("    infos: " + input.infos + "\n")
    f.write("    splices: " + input.splices + "\n")
    f.write("analyzers:\n")
    f.write("  SAM:\n")
    f.write("    options:\n")
    f.write("      max_hits: 1\n")
    f.write("      sampling_rate: 0.001\n")
    f.write("output:\n")
    f.write("  statistics: [Accuracy, Sensitivity, true-negatives, false-positives, false-negatives, true-positives, nb-elements]\n")
    f.write("softwares:\n")
    for i, name in enumerate(params.mapping_pipelines):
      f.write("  - name: " + name + "\n")
      f.write("    files:\n")
      f.write("      - name: " + input.bam_files[i] + "\n")
      f.write("        type: SAM\n")
      f.write("        check:\n")
      f.write("          - mapping\n")
      if name.find("star"):
        f.write("      - name: " + MAPPING_DIR + "/" + name + "/" + params.sample + "/SJ.out.tab\n")
        f.write("        type: STAR::Junction\n")
      elif name.find("crac"):
        f.write("      - name: " + MAPPING_DIR + "/" + name + "-splice.bed\n")
        f.write("        type: BED::Junction\n")
      elif name.find("hisat"):
        f.write("      - name: " + MAPPING_DIR + "/" + name + "_novel_splice.bed\n")
        f.write("        type: Hisat::Splice\n")
      else:
        # Otherwise we check splices from SAM file
        f.write("          - splice\n")
    f.close()

rule benchct_generic:
  input:
    conf = "{bench_dir}/{sample}.yaml",
  output: 
    bench = "{bench_dir}/{sample}.tsv",
  threads: 5
  shell: "benchCT -p {threads} -v {input.conf} > {output.bench}"

rule merge_true_positives:
  input: 
    bench = expand("{bench_dir}/{{sample}}/true-positives/{pipeline}-{{event}}.log",
                bench_dir = BENCHCT_MUTATION_DIR,
                pipeline = CALLING_PIPELINES)
  output: BENCHCT_MUTATION_DIR + "/{sample}/true-positives/{event}.tsv"
  version: "0.06"
  run: 
    shell("rm -f {output}")
    for i, name in enumerate(CALLING_PIPELINES):
      shell("cut -f1 " + input.bench[i] + "| awk '{{print \"" + name + "\",$1}}' >> {output}")

rule mutation_intersection_plot:
  input:  BENCHCT_MUTATION_DIR + "/{sample}/true-positives/{event}.tsv"
  output: FIGURES_DIR + "/{sample}/{event}-true-positives.pdf"
  version: "0.07"
  run:
    R("""
    library(UpSetR)
    dat <- read.table("{input}")
    mut <- as.data.frame.matrix(t(table(dat)))
    pdf(file = "{output}", width = 10, height = 5, onefile=FALSE)
    upset(mut, sets.bar.color = "#56B4E9", order.by = "freq", sets = colnames(mut))
    dev.off(which = dev.cur())
    """)

rule precision_recall_plot:
  input:
    bench = BENCHCT_DIR + "/{type}/{sample}.tsv",
  output:
    acc_sen = FIGURES_DIR + "/{sample}/{type}_accuracy_sensitivity.pdf",
  version: "0.06"
  run:
    R("""
    library(ggplot2)
    dat <- read.table("{input.bench}", header = TRUE)
    dat2 <- subset(dat, variable == "Sensitivity" | variable == "Accuracy")
    ggplot(dat2, aes(x=software,y=value,color=variable)) + 
      geom_point() + 
      facet_grid(~event, scales = "free") + 
      coord_flip() + 
      #ylim(0,1) + 
      theme_bw() + 
      scale_color_manual(values=c("#4E77AA", "#E69F00"))
    ggsave("{output.acc_sen}", width = 9, height = 3)
    """)


rule star_fusion:
  input:
    r1 = DATASET_DIR + "/{sample}/reads_1.fastq.gz",
    r2 = DATASET_DIR + "/{sample}/reads_2.fastq.gz",
  params:
    index = "/data/indexes/STAR/GRCh38_with_MT",
    chimSegmentMin = "{nb}"
  output:
    dir   = STAR_DIR + "_fusion/{sample}-{nb}chimSegmentMin/",
    chim  = STAR_DIR + "_fusion/{sample}-{nb}chimSegmentMin/Chimeric.out.junction",
    time  = STAR_DIR + "_fusion/{sample}-{nb}chimSegmentMin/time.txt"
  threads: NB_THREADS_MAPPERS
  run:
    if not os.path.exists(output.dir):
      os.makedirs(output.dir)
    shell("""/usr/bin/time -v -o {output.time} \
            {STAR} --genomeDir {params.index} --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat --twopassMode Basic --outSAMtype None \
            --outFileNamePrefix {output.dir} \
            --alignMatesGapMax {MAX_SPLICE_LENGTH} \
            --alignIntronMax {MAX_SPLICE_LENGTH} \
            --chimSegmentMin {params.chimSegmentMin} \
            --chimJunctionOverhangMin 12 \
            --chimSegmentReadGapMax 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --runThreadN {threads}""")


rule benchct_configfile_fusion:
  input: 
    infos =      DATASET_DIR + "/{sample}/info.txt",
    chimeras =   DATASET_DIR + "/{sample}/chimeras.tsv.gz",
    star_fusion = expand("{dir}_fusion/{sample}-{nb}chimSegmentMin/Chimeric.out.junction.post",
                         dir = STAR_DIR,
                         sample = "{sample}",
                         nb = CHIMSEGMENT_VALUES),
    crac_fusion = CRAC_DIR + "/{sample}-chimera.tsv"
  output: BENCHCT_FUSION_DIR + "/{sample}.yaml"
  version: "0.01"
  params:
    pipelines = expand("{nb}chimSegmentMin", nb = CHIMSEGMENT_VALUES),
    sample = "{sample}"
  version: "0.01"
  run:
    f = open(output[0], 'w')
    f.write("---\n")
    f.write("checker:\n")
    f.write("  files:\n")
    f.write("    infos: " + input.infos + "\n")
    f.write("    chimeras: " + input.chimeras + "\n")
    f.write("softwares:\n")
    for x in params.pipelines:
      f.write("  - name: " + x + "\n")
      f.write("    files:\n")
      f.write("      - name: " + STAR_DIR + "_fusion/" + params.sample + "-" + x + "/Chimeric.out.junction\n")
      f.write("        type: STAR::Chimera\n")
      f.write("        check: all\n")
    f.write("  - name: " + CRAC_NAME + "\n")
    f.write("    files:\n")
    f.write("      - name: " + input.crac_fusion + "\n")
    f.write("        type: CRAC::Chimera\n")
    f.write("        check: all\n")
    f.close()

rule star_fusion_post:
  input: "{sample}/Chimeric.out.junction"
  output: "{sample}/Chimeric.out.junction.post"
  params:
    min_recurrence = 2
  shell: """cat {input} | awk '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' | \
            sort | uniq -c | sort -k1,1rn | \
            awk 'BEGIN {{ OFS = "\t"}} $1 >= {params.min_recurrence} \
            {{print $2,$3,$4,$5,$6,$7,$8,$9,$10}}' > {output}"""
