############################################################
#                                                          #
#                         BPNS                             #
#                       Seasonal                           #
#                         MTX                              #
#                                                          #
############################################################

# This snakemake workflow starts from raw metatranscriptomic 
# samples and yields an annotated de novo metatranscriptome

# ------------------------- GOALS -------------------------#
# This script performs:
# - Quality control (fastQC, multiQC)
# - Merge sequences (for if downstream tools require those)
# - Trim the adapters from the raw PE sequences (trimmomatic)
# - Remove rRNA reads (ribodetector)
# - Run assemblies per sample for PE reads (rnaSPAdes)
# - Merge assemblies and cluster at 95% identity (MMseqs2)
# - Clean up assemblies (remove ERCC transcripts)
# - Map trimmed PE reads to the assembly using Kallisto
# - Predict proteins from long enough transcripts in clean 
#   assembly
# - Annotate final metatranscriptome for taxonomic and func-
#   tional information

# ------------------------ REMARKS ------------------------#

# ------------------------ SET-UP -------------------------#
# Load the config file with all the file paths to be loaded 
# into global variables below
configfile: 'config.yaml'

import os
import pandas as pd

# ------------------- GLOBAL VARIABLES --------------------#
input_dir = config['raw_read_directory']
scratch_dir = config['scratch']
output_dir = config['output_directory']
conda_dir = config['conda_environment_directory']
scripts_dir = config['scripts_directory']
adapters = config['adapters']
ERCC_folder = config['ERCC_folder']
sample_info = pd.read_csv(config['metaT_sample_metadata'], sep = ';')
samples = sample_info['sample'].tolist()

# ------------------------- RULES -------------------------#

rule all:
    input:
        fastqc_rawPE = expand('{base}qc/fastqc/{sample}_R{num}_fastqc.zip', 
                                base = output_dir,
                                sample= samples,
                                num = [1,2]),
        trimmed_seqs = expand('{base}trimmed/{sample}_R{num}.trimmed.fastq.gz', 
                                base = scratch_dir,
                                sample= samples,
                                num = [1,2]),
        merged_seqs = expand('{base}merged/{sample}_merged.fastq.gz', 
                                base = scratch_dir,
                                sample= samples),
        fastqc_merged = expand('{base}qc/fastqc/merged_reads/{sample}_fastqc.zip',
                                base = output_dir,
                                sample = samples),
        rRNA_removed_seqs = expand('{base}cleanup/{sample}_R{num}.rRNA_removed.fastq.gz',
                                base = scratch_dir,
                                sample = samples,
                                num = [1,2]),
        assembly_SPAdes = expand('{base}assembly/SPAdes/{sample}/transcripts.fasta',
                                base = output_dir,
                                sample = samples),
        metatranscriptome = output_dir + "assembly/SPAdes/metatranscriptome_renamed.fasta",
        proteins = output_dir + "assembly/protein/metatranscriptome.pep",
        kallisto = expand("{base}kallisto/{sample}",
                                base = output_dir,
                                sample = samples),
        kallisto_ERCC = expand("{base}ERCC92/kallisto/{sample}",
                                base = output_dir,
                                sample = samples),
        # For the below rules: let the whole snakefile run, then uncomment the merge rules!
        #merge_kallisto = output_dir + "kallisto/merged",
        #merge_kallisto_ERCC = output_dir + "ERCC92/kallisto/merged",

rule fastqc:
        input:
                input_dir + "{sample}_R{num}.fastq.gz"
        output:
                html = output_dir + 'qc/fastqc/{sample}_R{num}_fastqc.html',
                zip = output_dir + 'qc/fastqc/{sample}_R{num}_fastqc.zip'
        wrapper:
                "0.27.1/bio/fastqc"

rule trimmomatic:
        input:
                r1 = input_dir + "{sample}_R1.fastq.gz",
                r2 = input_dir + "{sample}_R2.fastq.gz",
        output:
                r1 = os.path.join(scratch_dir, "trimmed", "{sample}_R1.trimmed.fastq.gz"),
                r2 = os.path.join(scratch_dir, "trimmed", "{sample}_R2.trimmed.fastq.gz"),
                # reads where trimming entirely removed the mate
                r1_unpaired = os.path.join(scratch_dir, "trimmed", "{sample}_R1.unpaired.fastq.gz"),
                r2_unpaired = os.path.join(scratch_dir, "trimmed", "{sample}_R2.unpaired.fastq.gz")
        params:
                trimmer=["ILLUMINACLIP:{}:2:30:7".format(adapters), 
                        "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:50"]
        #conda: conda_dir + "trimmomatic.yaml"
        shell:'''
        module load  Trimmomatic/0.39-Java-11
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
                {input.r1} {input.r2} \
                {output.r1} {output.r1_unpaired} \
                {output.r2} {output.r2_unpaired} \
                {params.trimmer}
        '''

rule merge_reads:
        input:
                r1 = input_dir + "{sample}_R1.fastq.gz",
                r2 = input_dir + "{sample}_R2.fastq.gz"
        output:
                r_merged = scratch_dir + "merged/{sample}_merged.fastq.gz"
        shell:'''
        module load BBMap
        bbmerge-auto.sh in1={input.r1} in2={input.r2} \
            out={output.r_merged} adapters={adapters} \
            rem k=62 extend2=50 ecct -Xmx40g
        '''

rule fastqc_merged:
        input:
                scratch_dir + "merged/{sample}_merged.fastq.gz"
        output:
                html = output_dir + 'qc/fastqc/merged_reads/{sample}_fastqc.html',
                zip = output_dir + 'qc/fastqc/merged_reads/{sample}_fastqc.zip'
        wrapper:
                "0.27.1/bio/fastqc"

rule sample_rRNA_cleanup:
    '''
    This step removes ribosomal RNA reads using ribodetector.
    The 
    '''
    input:
        r1 = scratch_dir + 'trimmed/{sample}_R1.trimmed.fastq.gz',
        r2 = scratch_dir + 'trimmed/{sample}_R2.trimmed.fastq.gz'
    output:
        r1 = scratch_dir + 'cleanup/{sample}_R1.rRNA_removed.fastq.gz',
        r2 = scratch_dir + 'cleanup/{sample}_R2.rRNA_removed.fastq.gz',
        rna1 = scratch_dir + 'rRNA/{sample}_R1.rRNA.fastq.gz',
        rna2 = scratch_dir + 'rRNA/{sample}_R2.rRNA.fastq.gz'
    conda: conda_dir + 'ribodetector.yaml'
    shell:'''
    ribodetector_cpu -t 20 \
        -i {input.r1} {input.r2} \
        -o {output.r1} {output.r2} \
        -r {output.rna1} {output.rna2} \
        -e rrna -l 140 --chunk_size 256
    '''

rule rnaSPAdes:
        input:
                r1 = scratch_dir + "cleanup/{sample}_R1.rRNA_removed.fastq.gz",
                r2 = scratch_dir + "cleanup/{sample}_R2.rRNA_removed.fastq.gz" 
        output:
                dir = directory(output_dir + "assembly/SPAdes/{sample}"),
                transcripts = output_dir + "assembly/SPAdes/{sample}/transcripts.fasta"
        # conda: conda_dir + "spades.yaml"
        params:
                k = '75,99,127'
        shell:'''
        module load SPAdes/3.15.3-GCC-11.2.0
        spades.py --rna \
                -1 {input.r1} \
                -2 {input.r2} \
                -o {output.dir} \
                -k {params.k}
        '''

rule merge_assemblies:
        input: expand(output_dir + "assembly/SPAdes/{sample}/transcripts.fasta", sample = samples)
        output: temporary(output_dir + "assembly/SPAdes/all_transcripts.fasta")
        shell:'''
         cat {input} >> {output}
        '''

rule cluster_assembly:
        input:
                output_dir + "assembly/SPAdes/all_transcripts.fasta"
        output:
                mmseqs_clu = scratch_dir + "mmseqs_clustering/clusterDB",
                tmp = scratch_dir + "mmseqs_clustering/scratch/",
                metatranscriptome = temporary(output_dir + "assembly/SPAdes/metatranscriptome.fasta"),
        params:
                percent_id = "0.95",
                  
        shell:'''
        unset OMP_PROC_BIND

        module load MMseqs2
        
        mmseqs easy-linclust {input} {output.mmseqs_clu} {output.tmp} --min-seq-id {params.percent_id}

        mv {output.mmseqs_clu}_rep_seq.fasta {output.metatranscriptome}
        '''

rule rename_transcripts:
        input:  output_dir + "assembly/SPAdes/metatranscriptome.fasta"
        output:  temporary(output_dir + "assembly/SPAdes/metatranscriptome_renamed.fasta")
        #conda: conda_dir + "anvio.yaml"
        # installation problems with conda, solved through module system
        shell:'''
        module load anvio

        anvi-script-reformat-fasta {input} -o {output} --simplify-names
        '''

rule identify_ERCCs:
        '''
        Remove the ERCC transcripts from the final assembly.
        
        Step 1
        ------
        Search for sequences in assembly that map to the ERCC reference file

        Step 2
        ------snak
        Gather the names from the contigs that identify as an ERCC sequence
        '''
        input: 
                ERCC = ERCC_folder + 'ERCC92.fa',
                metatranscriptome = output_dir + "assembly/SPAdes/metatranscriptome_renamed.fasta",
        output: 
                ERCC_matches = ERCC_folder + "ERCC_assembled_contigs.m8",
                ERCC_match_names = ERCC_folder + "ERCC_names_in_metaT.txt"
        params:
                tmp = scratch_dir + 'mmseqs'
        shell:"""
        unset OMP_PROC_BIND

        module load MMseqs2

        mmseqs easy-search {input.ERCC} {input.metatranscriptome} \
        {output.ERCC_matches} {params.tmp} --search-type 3

        cat {output.ERCC_matches} | cut -f2 | sed 's/^/>/' > {output.ERCC_match_names}
         """

COMPARE_HEADERS_CMD=r"""NR==FNR{a[$0];next} /^>/{f=($0 in a ? 1 : 0)} !f"""

rule remove_ERCCs:
        '''
        Remove the ERCC transcripts from the final assembly.

        Step 3
        ------
        Remove those sequences from the assembly before protein prediction, annotation, and read mapping
        '''
        input:
                metatranscriptome = output_dir + "assembly/SPAdes/metatranscriptome_renamed.fasta",
                ERCC_matches = ERCC_folder + "ERCC_assembled_contigs.m8",
                ERCC_match_names = ERCC_folder + "ERCC_names_in_metaT.txt"
        output: 
                output_dir + "assembly/SPAdes/final_metatranscriptome.fasta"
        shell:"""
        awk {COMPARE_HEADERS_CMD:q} {input.ERCC_match_names} {input.metatranscriptome} > {output}
        """

rule protein_prediction:
        input: 
                metaT = output_dir + "assembly/SPAdes/final_metatranscriptome.fasta",
                
        output:
                bed = output_dir + "assembly/protein/metatranscriptome.bed",
                cds = output_dir + "assembly/protein/metatranscriptome.cds",
                gff3 = output_dir + "assembly/protein/metatranscriptome.gff3",
                pep =  output_dir + "assembly/protein/metatranscriptome.pep",
        params:
                longORFs_dir = output_dir + "assembly/protein",
                dir = "final_metatranscriptome.fasta.transdecoder_dir",
                bed = "final_metatranscriptome.fasta.transdecoder.bed",
                cds = "final_metatranscriptome.fasta.transdecoder.cds",
                gff3 = "final_metatranscriptome.fasta.transdecoder.gff3",
                pep = "final_metatranscriptome.fasta.transdecoder.pep",
                checkpoints = "final_metatranscriptome.fasta.transdecoder_dir.__checkpoints"
        shell:'''
        cd {params.longORFs_dir}

        module load TransDecoder

        TransDecoder.LongOrfs -t {input.metaT} \
                --output_dir {params.longORFs_dir}

        TransDecoder.Predict -t {input.metaT} \
                --output_dir {params.longORFs_dir} --single_best_only

        mv {params.longORFs_dir}/{params.bed} {output.bed}
        mv {params.longORFs_dir}/{params.cds} {output.cds}
        mv {params.longORFs_dir}/{params.gff3} {output.gff3}
        mv {params.longORFs_dir}/{params.pep} {output.pep}
        # rm -rf {params.longORFs_dir}/{params.dir} {params.longORFs_dir}/{params.checkpoints}
        '''        

rule kallisto_index:
    input:
       fasta = output_dir + "assembly/SPAdes/final_metatranscriptome.fasta",
    output:
        index = output_dir + "assembly/SPAdes/final_metatranscriptome.idx",
    params:
        extra = "",  # optional parameters
    threads: 1
    wrapper:
        "v1.15.0/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq = [scratch_dir + 'cleanup/{sample}_R1.rRNA_removed.fastq.gz',
                scratch_dir + 'cleanup/{sample}_R2.rRNA_removed.fastq.gz'],
        index = output_dir + "assembly/SPAdes/final_metatranscriptome.idx",
    output:
        directory(output_dir + "kallisto/{sample}"),
    params:
        extra="",
    log:
        "data/logs/kallisto_quant_{sample}.log",
    threads: 1
    wrapper:
        "v1.15.0/bio/kallisto/quant"

rule merge_kallisto:
        input: output_dir + "kallisto/"
        output: directory(output_dir + "kallisto/merged/")
        shell:'''
        python3 {scripts_dir}merge_kallisto.py -i {input} -o {output}
        '''

rule kallisto_index_ERCC:
    input:
       fasta =  ERCC_folder + 'ERCC92.fa',
    output:
        index = ERCC_folder + "ERCC92.idx",
    params:
        extra = "",  # optional parameters
    threads: 1
    wrapper:
        "v1.15.0/bio/kallisto/index"

rule kallisto_quant_ERCC:
    input:
        fastq = [scratch_dir + 'cleanup/{sample}_R1.rRNA_removed.fastq.gz',
                scratch_dir + 'cleanup/{sample}_R2.rRNA_removed.fastq.gz'],
        index = ERCC_folder + "ERCC92.idx",
    output:
        directory(output_dir + "ERCC92/kallisto/{sample}"),
    params:
        extra="",
    log:
        "data/logs/kallisto_quant_ERCC_{sample}.log",
    threads: 1
    wrapper:
        "v1.15.0/bio/kallisto/quant"

rule merge_kallisto_ERCC:
        input: output_dir + "ERCC92/kallisto/"
        output: directory(output_dir + "ERCC92/kallisto/merged")
        shell:'''
        python3 {scripts_dir}merge_kallisto.py -i {input} -o {output}
        '''