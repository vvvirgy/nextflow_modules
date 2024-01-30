#setwd("Z:/nextflow_module/results") # see if this stands anyway

require(CNAqc)
require(tidyverse)

# get the CNAqc best solution for each sample and patient from the table and the nextflow script 
#get the CNAqc output file names

#patients = $patient
#samples = $sample

patients = list.files() #$patient
#cohort = $datasetID
#sample = $sampleID

dir.create("join_table")

bestsol_path = sapply(patients, function(x) {
  files = list.files(x, full.names = T, recursive = T)
  grep(pattern = "*SEQUENZA_CNAqc/best_solution.rds", files, value = T)
  }
) 

## in the future: the sample id will be taken directly from the input csv

#read the results and store them in a list

results = list() # includes all results from both sequenza and cnaqc

cnaqc_segments = list() ## includes only cnaqc segments
cnaqc_mutations = list()

for (i in names(bestsol_path)){
  x = lapply(bestsol_path[[i]],readRDS) ## list with all the rds from one patient 
  names(x) = sample_id[[i]]
  results[[i]] = x
  
  patient_cnaqc = lapply(x, function(y) {
    y$cnaqc[[1]]
  })
  
  cnaqc_segments[[i]] = lapply(patient_cnaqc, function(x) {
    x_purity = x$purity
    x_ploidy = x$ploidy
    
    x$cna %>% 
      mutate(purity = x_purity) %>% 
      mutate(ploidy = x_ploidy)
  })

  cnaqc_mutations[[i]] = lapply(patient_cnaqc, function(x) {
    x_purity = x$purity
    x_ploidy = x$ploidy

    x$mutations %>% 
      mutate(purity = x_purity) %>% 
      mutate(ploidy = x_ploidy) %>% 
      mutate(normal_cn = 2)
  })
}


###################################

## add the sample_id column to the cna object

for (i in seq(names(cnaqc_segments))) {
  samples = names(cnaqc_segments)[[i]]
  
  id = sample_id[[samples]]
  cnaqc_segments[[samples]] = lapply(id, function(x) {
    cnaqc_segments[[samples]][[x]] %>%
      mutate(sample_id = x)
  })
}

## merge all samples together in a unique table of cn

tmp_cna = lapply(cnaqc_segments, function(x){do.call(rbind, x)})
piled_up_table_cna = do.call(rbind, tmp_cna)

# this is gonna take a while 

### step 1 = modifications on the mutations table

# add the sample_id column to it

for(i in seq(names(cnaqc_mutations))) {
  samples = names(cnaqc_mutations)[[i]]

  id = sample_id[[samples]]
  cnaqc_mutations[[samples]] = lapply(id, function(x) {
    cnaqc_mutations[[samples]][[x]] %>% 
      mutate(sample_id = x)
  }) 
}

# merge the mutation tables in a single one 

tmp_mut = lapply(cnaqc_mutations, function(x){do.call(rbind, x)})
piled_up_mutations = do.call(rbind, tmp_mut)

### step 2 = get the new segmentation and map the mutations on the new segments

source("joint_table_functions.R") ### change the name of the file! 

joint_new_segments = new_segmentation(piled_up_table_cna, piled_up_mutations)

# NB: the table produced in this way, considering the starting samples, is not including the information on the gene id atm

## step 3 = fix the format of the mutation table

joint_new_segments$mutation = get_right_mutTable(joint_new_segments$mutation)


### save new segmentation table and new mutation table 

saveRDS(joint_new_segments, file = "results/join_table/piled_up_table_updated.rds") # object with both cna and mutations

saveRDS(joint_new_segments$mutation, file = "results/join_table/joint_mutation_table.rds")
saveRDS(joint_new_segments$cna, file = "results/join_table/joint_cna_table.rds")

write_tsv(joint_new_segments$mutation, file = "results/join_table/join_mutation_table.tsv")
write_tsv(joint_new_segments$cna, path = "results/join_table/join_cna_table.tsv")