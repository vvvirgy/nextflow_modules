new_segmentation = function(piled_up_cna, mut_table) {
  
  #divide the table by chromosome
  
  chromosomes = unique(piled_up_cna$chr)
  
  piledbyChr = lapply(chromosomes, function(x) {
    tmp = piled_up_cna %>% 
      filter(chr == x)
  }) 
  
  names(piledbyChr) = chromosomes

  piledbyChr_updated <- list()
  piledbyChr_updated_mut <- list()

  for (chromo in chromosomes) {
    
    # get the table with the segment of the chromosome we are iterating over
    chromo_df = piledbyChr[[chromo]]

    samples = unique(chromo_df$sample_id)

    # get the indicator to end the segmentation process, which is the maximun value of the to column (end of the last original segment)

    max_global = max(sapply(samples, function(x) {
      chromo_df %>% 
        filter(sample_id == x) %>%
        pull(to) %>% 
        max()
      }))

    # get the min of from, which is the lower bound of the new segmentation, and the from of the new first segment

    min_from_list = sapply(samples, function(x) {
    chromo_df %>% 
      filter(sample_id == x) %>%
      pull(from) %>% 
      min()
    })

    all_min_total = min(min_from_list)
    all_min = all_min_total 
    
    # generate the vectors in which store the breakpoints

    from_breakpoints = all_min
    to_breakpoints = c()
    new_segments = c()

    ### start the iteration

    repeat {

      new_min_from_list = sapply(samples, function(x) {
        chromo_df %>% 
          filter(from > all_min) %>% 
          filter(sample_id == x) %>%
          select(from,to) %>% 
          min()
        })

      new_min = min(new_min_from_list)

      if(new_min == Inf) {
        new_min_from_list = sapply(samples, function(x) {
        chromo_df %>% 
          filter(to > all_min) %>% 
          filter(sample_id == x) %>%
          select(to) %>% 
          min()
        })

        new_min = min(new_min_from_list)
      }

      from_breakpoints = append(from_breakpoints, new_min)
      to_breakpoints = append(to_breakpoints, new_min)
      new_segments = append(new_segments, paste(chromo, all_min, new_min, sep = ":"))

      all_min = new_min

      print("going to next breakpoint")

      #check_to = sapply(samples, function(x) {
      #  chromo_df %>% 
      #    filter(to >= all_min) %>% 
      #    filter(sample_id == x) %>%
      #    select(to) %>% 
      #    min()
      #  }) %>% min()

      if(max_global == all_min) {

        print("reached the upper limit of the last breakpoint")

        #to_breakpoints = append(to_breakpoints, max_global)
        #new_segments = append(new_segments, paste(chromo, all_min, max_global, sep = ":"))

        break
      }
    }

    join_seg_table <- list()
    mut_table_new = list()

    for (i in seq(new_segments)) {

      from_n = from_breakpoints[i]
      to_n = to_breakpoints[i]

      tmp = chromo_df %>% 
        filter(from <= from_n & to >= to_n) %>% 
        mutate(from_new = from_n) %>% 
        mutate(to_new = to_n) %>% 
        mutate(updated_segment_id = #case_when(
      #(from <= from_n & to >= to_n) ~ 
            paste(chr, from_n, to_n, Major, minor, sep = ":")
      #)
          ) #%>% 
        #mutate(karyotype = #case_when(
      #(from <= from_n & to >= to_n) ~ 
            #paste(Major, minor, sep = ":"))#) 

      join_seg_table[[i]] <- tmp 
#
      tmp_mut = mut_table %>% 
        filter(from >= from_n & to <= to_n) %>% 
        mutate(new_segment_id = 
          paste(chr, from_n, to_n, karyotype, sep = ":")) %>% 
        mutate(mutation_id = paste(chr, from, to, ref, alt, sep = ":"))

      mut_table_new[[i]] <- tmp_mut

    }

    piled_up_table_updated = do.call(rbind, join_seg_table)
    piledbyChr_updated[[chromo]] <- piled_up_table_updated

    piled_up_mutations_updated = do.call(rbind, mut_table_new)
    piledbyChr_updated_mut[[chromo]] <- piled_up_mutations_updated
  }

  piled_up_cna_updated = do.call(rbind, piledbyChr_updated)
  piled_up_mut_updated = do.call(rbind, piledbyChr_updated_mut)

  Piled_up = list(mutation = piled_up_mut_updated, cna = piled_up_cna_updated)

  return(Piled_up)
}


############ get the mutation table as we want

get_right_mutTable <- function(mutations) { # takes mutations table and modifies the karyotype to get just the major and minor columns

  all_karyo = strsplit(mutations$karyotype, split = ":")

  Major = unlist(lapply(all_karyo, function(x) {as.numeric(x[[1]])}))
  minor = unlist(lapply(all_karyo, function(x) {as.numeric(x[[2]])}))

  mutation_new = mutations %>%
    #filter(from >= from & to <= to_seg) %>% 
    mutate(Major = Major,
         minor = minor) %>% 
    relocate(Major, .after = DP) %>% 
    relocate(minor, .after = Major) %>% 
    relocate(mutation_id, .before = chr) %>% 
    relocate(new_segment_id, .after = mutation_id) %>% 
    relocate(normal_cn, .after = minor) %>% 
    mutate(segment_id = NULL)
  
  return(mutation_new)
}