#!/user/bin/env Rscript

## rreggiar@ucsc.edu

## testing job functionaliyt in Rstudio for running big chunks in background

if (!file.exists(file.path(output.data.dir, 'tcga.data', 'luad.data.frame.tsv'))){

	print('hello')
    
  # toil.survival <- read_tsv(file.path(output.data.dir, 'tcga.data/TCGA-LUAD.survival.tsv.gz')) 
  
  # luad.phenotypes <- read_tsv(file.path(output.data.dir, 'tcga.data/TCGA-LUAD.GDC_phenotype.tsv.gz')) %>% 
  #   select(sample = submitter_id.samples, 
  #          sample_type = sample_type.samples,
  #          age = age_at_index.demographic,
  #          sex = gender.demographic,
  #          race = race.demographic)
  
  # luad.mutect <- read_tsv(file.path(output.data.dir, 'tcga.data/TCGA-LUAD.mutect2_snv.tsv.gz')) %>% 
  #   filter(gene %in% c('KRAS')) %>% 
  #   select(sample = Sample_ID, gene, aa_change = Amino_Acid_Change) %>% 
  #   pivot_wider(names_from = gene, values_from = aa_change) %>% 
  #   mutate(across(.cols = where(is.list), paste)) 
  
  # luad.phenotypes.merge <- luad.phenotypes %>% 
  #   merge(luad.mutect, by = 'sample', all.x = T) %>% 
  #   filter(sample_type %in% c('Primary Tumor', 'Solid Tissue Normal')) %>% 
  #   merge(toil.survival, by = 'sample') %>% 
  #   mutate(sample = str_sub(sample, end = str_length(sample) - 1)) %>% 
  #   replace_na(list(KRAS = 'WT')) %>% 
  #   mutate(KRAS = ifelse(sample_type == 'Solid Tissue Normal', 'WT_normal', KRAS))

  # toil.norm.counts <- read_csv(file.path(output.data.dir, 'tcga.data/toil.counts.csv'))
  
  # norm.counts.tidy <- 
  #   toil.norm.counts %>% 
  #   select('ensg' = gene, contains('TCGA')) %>% 
  #   discard(~sum(.x == 0)/length(.x) >= 0.8) %>% 
  #   gather(sample, count, -ensg) 
  
  # rm(toil.norm.counts)
  
  # norm.counts.tidy.filt <- 
  #   norm.counts.tidy %>% 
  #   filter(sample %in% luad.phenotypes.merge$sample)
  
  # norm.counts.tidy.filt.merge <- 
  #   norm.counts.tidy.filt %>% 
  #   merge(gen.24.names, by = 'ensg') %>%
  #   select(-ensg)
  
  # rm(norm.counts.tidy.filt)
  
  # norm.counts.final <- 
  #   norm.counts.tidy.filt.merge %>% 
  #   merge(luad.phenotypes.merge, by = 'sample') %>% 
  #   distinct()
  
  # write_tsv(norm.counts.final, file = file.path(output.data.dir, 'tcga.data', 'luad.data.frame.tsv'))
  
} else if (file.exists(file.path(output.data.dir, 'tcga.data', 'luad.data.frame.tsv'))) {
  # norm.counts.final <- read_tsv(file.path(output.data.dir, 'tcga.data', 'luad.data.frame.tsv'))
		print(head(gen.34.lnc.names))

}