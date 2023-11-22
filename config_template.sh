# BPNS_seasonal_MTX/config.sh

# General directories
work_dir='BPNS_seasonal_MTX'
tmp='tmp'
conda_env='envs/'

# Final assembled metatranscriptome and protein translations
assembly=${work_dir}'/data/assembly/SPAdes/final_metatranscriptome.fasta'
assembly_dir=${work_dir}'data/assembly/SPAdes/'
metatranscriptome_mmseqsDB=${work_dir}'/data/assembly/SPAdes/final_metatranscriptome_mmseqsDB'
assembly_prot_dir=${work_dir}'/data/assembly/protein/'
assembly_pep='metatranscriptome.pep'
assembly_mmseqsDB='SPAdes_mmseqsDB'

# EggNOG functional reference
eggnog_ref_dir='ref_dbs/eggnog/'
eggnog_out_dir=${work_dir}'/data/annotation/functional_eggnog/'

# EukProt reference
eukprot_ref_dir='ref_dbs/eukprot/'
eukprof_ref_fasta='eukprot.reference.fasta.gz'
eukprot_reference='eukprot'
eukprot_out_dir=${work_dir}'/data/annotation/taxonomy_eukprot/'
eukprot_mmseqsDB='eukprot_DB'

# MMETSP reference
mmetsp_ref_dir='references/mmetsp/marmmetsp/'
mmetsp_ref_fasta='reference.pep.fa'
mmetsp_reference='marmmetsp'
mmetsp_out_dir=${work_dir}'/data/annotation/taxonomy_MMETSP/'
mmetsp_mmseqsDB='marmmetsp'

# PhyloDB and extension references
phylodb_ref_dir='references/phylodb/'
phylodb_ref_fasta='phylodb_1.076.pep.fa.gz'
phylodb_extended_ref_fasta='phylodb_1.076_extended.pep.fa.gz'
phylodb_reference='phylodb_1.076'
phylodb_extended_reference='phylodb_1.076_extended'
phylodb_out_dir=${work_dir}'/data/annotation/taxonomy_phyloDB/'
phylodb_mmseqsDB='phylodb_nobacktrace'
phylodb_extended_out_dir=${work_dir}'/data/annotation/taxonomy_phyloDB_extended/'
phylodb_extended_mmseqsDB='phylodb_extended'
phylodb_extra_transcriptome_dir=${work_dir}'/data/annotation/taxonomy_phyloDB/transcriptomes/'

# TARA MATOU and MAG references
MATOU='ref_dbs/MATOU/MATOU-v1.fa'
MATOU_mmseqsDB='ref_dbs/MATOU/MATOU_mmseqsDB'
MATOU_output=${work_dir}'data/MATOU_alignment/MATOU'
TARA_MAG_ref='ref_dbs/delmont_mags/MAG_contigs/MAGS_combined.fa'
TARA_MAG_out_dir=${work_dir}'data/SMAG_alignment/minimap/'
MAG_bamfilters_dir=${work_dir}'data/SMAG_alignment/bamfilters/'