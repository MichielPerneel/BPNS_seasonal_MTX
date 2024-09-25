<p align="center">
  <h1>
  A Year-long Metatranscriptomic Study in the North Sea
  </h1>
  A bioinformatic pipeline to process monthly raw metatranscriptomic samples from microeukaryotic plankton into interpretable data. This is the code used to analyse and visualize the data for [this publication](https://journals.asm.org/doi/10.1128/mbio.00383-24)
  </p>
<br/>


<p align="center">
  <a href="#about">About</a> •
  <a href="#data">Data</a> •
  <a href="#installation">Installation</a> •
  <a href="#structure">Project Structure</a> •
  <a href="#deploy snakemake">Deploy Snakemake</a> •
  <a href="#annotate-the-assembled-transcripts">Transcript Annotation</a> •
  <a href="#data-analysis">Data Analysis</a> •
</p>  

## About
Marine microeukaryotic plankton play key roles in global carbon and nutrient cycling, and knowledge of their diverse taxonomy and functional repertoire has increased substantially in the last decade largely due to global -omics surveys. However, few studies have been done to track seasonal turnover in microeukaryotic functional diversity at fixed locations. With this study we present a seasonal metatranscriptomic dataset based on monthly samples from July 2020 to July 2020, sampled from 6 fixed stations in the Southern North Sea. We observe the phenological and metabolic dynamics of a coastal microeukaryotic plankton assemblage in the region. Our analysis identifies key seasonal shifts in biodiversity and functional richness, underscoring the metabolic traits of dominant taxonomic groups. The data also elucidate spatial variations resulting from the convergence of oceanic and estuarine waters and highlight variations in trophic feeding modes of specific microeukaryotes. 

## Data
The raw sequence data needs to be stored in the `data/raw` folder to be able to run the Snakemake workflow. The workflow is designed to handle PE RNA-seq reads from mixed environmental communities.

Raw sequence data can be found on [SRA](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1021244) under BioProject PRJNA1021244.

## Installation
Fork this repository from github to an HPC environment. Make sure there is sufficient storage (A total of 1.7 TB was needed to analyse 62 samples of 25M reads - 281 GB of .fastq.gz files).  First, modify the [config file](./config.yaml) if necessary. Second, install the bioinformatic tools used in the separate rules through conda, and export their environments to .yaml files that are placed in the conda environment folder. Or check whether the modules can be loaded on your HPC (by including a load command in the rules' shell script). Third, have a look at the [sample csv](./samples.csv) with the following structure:

| sample    | date       | time  | month        | station |
| :-------: | :--------: | :---: | :----------: | :-----: |
| 1_130_S31 | 18/01/2021 | 9:05  | January_2021 | 130     |
| 1_330_S34 | 18/01/2021 | 14:05 | January_2021 | 330     |
| 1_700_S32 | 18/01/2021 | 11:45 | January_2021 | 700     |
| 1_780_S33 | 18/01/2021 | 13:00 | January_2021 | 780     |
| ... | ... | ... | ... | ... |

The data and figures folder is added to the .gitignore file to avoid pushing big result files to github.

## Structure
Before running the pipeline, the repository contains:

```bash
├── config.yaml
├── config.sh
├── data
│   ├── ERCC92
│   ├── illumina-adapters.fa
│   └── raw
│── figures
├── hpc_config
│   ├── cluster.yaml
│   └── config.yaml
├── README.md
├── samples.csv
├── scripts
│   ├── analysis
│   ├── trophic-mode-ml
│   ├── eggNOG_functional_annotation.pbs
│   ├── get_read_counts.pbs
│   ├── launch_snakemake.pbs
│   ├── merge_kallisto.py
│   ├── MMETSP_taxonomic_annotation.pbs
│   ├── eukprot_taxonomic_annotation.pbs
│   ├── phyloDB_taxonomic_annotation.pbs
│   ├── create_phyloDB_extended_annotation_file.sh
│   ├── phyloDB_expansion.ipynb
│   ├── phyloDB_extended_taxonomic_annotation.pbs
│   ├── transdecoder_extra_phyloDB_transcriptomes.pbs
│   ├── variance_stabilizing_transformation.R
│   ├── simka.pbs
│   ├── rename_fasta_headers.sh
│   ├── map_bpns.R
│   ├── get_sequence_lengths.rs
│   ├── Tara_MATOU_mapping.pbs
│   └── Tara_MAG_mapping.pbs
└── snakefile
```

## Deploy snakemake
To run the Snakemake workflow on the [VSC HPC](https://www.vscentrum.be), submit the launch snakemake script:

```sh
qsub scripts/launch_snakemake.pbs
```

The script deploys snakemake and submits jobs to the correct nodes specified in [the cluster configuration](./hpc_config/cluster.yaml) and according to the rules in the [snakefile](./snakefile). The cluster deployment is configured in this [config.yaml](./hpc_config/config.yaml) file. After running Snakemake the following folders and files should have generated:

```bash
├── config.yaml
├── data
│   ├── assembly
│   ├── ERCC92
│   ├── illumina-adapters.fa
│   ├── kallisto
│   ├── logs
│   ├── qc
│   ├── raw
│   └── scratch
├── hpc_config
│   ├── cluster.yaml
│   └── config.yaml
├── README.md
├── samples.csv
├── scripts
└── snakefile
```

## Annotate the assembled transcripts
First, generate the annotation folder to place results in:

```sh
mkdir data/annotation
```
There should exist a file called config.sh in the folder that specifies the paths to the different reference databases, an example can be found in [config_template.sh](config_template.sh) file. 

To functionally annotate the proteins predicted from the assembled transcripts, run:

```sh
mkdir data/annotation/functional_eggnog
qsub -m abe scripts/eggNOG_functional_annotation.pbs
```

To taxonomically annotate the proteins predicted from the assembled transcripts, run on a heavy enough cluster:

```sh
# MMETSP
mkdir data/annotation/taxonomy_MMETSP
qsub -m abe scripts/MMETSP_taxonomic_annotation.pbs

# PhyloDB
mkdir data/annotation/taxonomy_phyloDB
qsub -m abe scripts/phyloDB_taxonomic_annotation.pbs

# Extended version of PhyloDB
## Check out and run scripts/phyloDB_expansion.ipynb first
mkdir data/annotation/taxonomy_phylodb_extended
qsub -m abe scripts/phyloDB_extended_taxonomic_annotation.pbs

# EukProt
mkdir data/annotation/taxonomy_eukprot
qsub -m abe scripts/eukprot_taxonomic_annotation.pbs
```

## Data analysis
Now one can run these jupyter notebooks (locally) to generate the results discussed below.

First, an analysis of general expression patterns can be performed in [this notebook](./scripts/analysis/Expression_analysis.ipynb). The notebook plots the annotation rates of the metatranscriptome, the distribution of transcript length, and the level of expression vs the variance in expression. 

K-mer analysis of the assemblies is performed using [simka](./scripts/simka.pbs). The results can be visualised with [this notebook](./scripts/analysis/Simka_analysis.ipynb).

Representation of the reads in the Tara Oceans SMAGs can be checked using [this script](scripts/Tara_MAG_mapping.pbs).
Beware that this needs the SMAG data installed, and the location of that database specified in the config.sh file.

Analysis of ERCC92 spike-in reads and how to normalise samples can be found in [this notebook](./scripts/analysis/SpikeIn_analysis_Normalisation.ipynb).

Environmental data is visualized in [this notebook](./scripts/analysis/Environmental_analysis.ipynb).

Then, the calculated transcripts per L are visualised in [here](./scripts/analysis/Transcripts_per_L.ipynb).

When above notebooks ran, an overall [taxonomic analysis](./scripts/analysis/Taxonomic_analysis.ipynb) was done followed by [a functional analysis](./scripts/analysis/Functional_analysis.ipynb). The taxonomic results were verified using [FlowCam data](./scripts/analysis/Flowcam_Analysis.ipynb). We then related [functional to taxonomic data](./scripts/analysis/Taxonomic_vs_Funcional_Diversity.ipynb).

Based on the functional annotations, we mined the data for co-occurring functional identifiers using [WGCNA](./scripts/analysis/wgcna.R). Co-occuring functions are separated into distinct modules. The resulting eigengenes of the found modules of co-occurring functional identifiers were [correlated with taxonomic occurrence](./scripts/analysis/module_taxonomy_correlation.r). Using WGCNA, we identify important gene functions that closely follow a module's overall expression. To get an idea of metabolic functions that are typical of the obtained modules, we implemented [ranked Mann-Whitney U tests](./scripts/analysis/Functional_Module_Description.ipynb) to statistically test for enriched pathways in a given modules.

Pathways were analysed [here](./scripts/analysis/Functional_Pathway_Description.ipynb) and statistically described using [GAM regressions](./scripts/analysis/Pathway_GAM.r).

We zoomed in on the niches of diatoms, dinoflagellates, and *Phaeocystis* in [this notebook](./scripts/analysis/Diatoms_Dinoflagellates_and_Prymnesiophytes.ipynb).

Finally, we [implemented](./scripts/analysis/trophic_mode_prediction.ipynb) the machine learning prediction model from [Lambert et al.](https://www.pnas.org/doi/full/10.1073/pnas.2100916119) to classify our taxonomic species bins into phototrophs, heterotrophs, and mixotrophs.
