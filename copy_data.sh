# download sample info
#aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/phenotype/vidus_hiv_acquisition_phenotype_variables.txt.gz
 aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/phenotype/vidus_hiv_acquisition_phenotype_variables.txt.gz \
    /shared/vidus/


#aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/phenotype/vidus_rna_seq_technical_variables.txt.gz
 aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/phenotype/vidus_rna_seq_technical_variables.txt.gz \
    /shared/vidus/

aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/phenotype/eigenstrat/vidus_hiv_acquisition_vl_suppressed_samples_top10_genotype_pcs.txt \
    /shared/vidus/
#aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/phenotype/eigenstrat/vidus_hiv_acquisition_all_samples_top10_genotype_pcs.txt \
    /shared/vidus/
gunzip /shared/vidus/vidus*variables.txt.gz


aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/gene_expression/deconvolution/cibersortx/hiv_status_candidate_gene_dge_cell_type_proportions_lm22.txt \
    /shared/vidus/

# gene expression data
aws s3 cp s3://rti-hiv/differential_gene_expression/vidus/data/gene_expression/whole_blood/eur/vidus_salmon_gene_data_gencode28.rds \
    /shared/vidus/
