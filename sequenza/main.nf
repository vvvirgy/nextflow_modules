// Extract
 
process SEQUENZA_EXTRACT {

    input:
    path seqzFile

    output:
    var

    script:
    """
    #!/usr/bin/env Rscript
    
    library(sequenza)
   
    seqzExt <-
    sequenza.extract(
      file = $seqzFile,
      chromosome.list = chromosomes,
      normalization.method = 'median',
      window = 1e5,
      gamma = 280,
      kmin = 300,
      min.reads.baf = 50,
      min.reads = 50,
      min.reads.normal = 15,
      max.mut.types = 1
  )
    """
}

// Fit

process SEQUENZA_FIT {

    input:
    path seqzExt

    script:
    """
    #!/usr/bin/env Rscript

    library(sequenza)
    
     PARSE SEQZEXT FILE
    
    paraSpace <-
    sequenza.fit(
      sequenza.extract = $seqzExt,
      cellularity = seq(low_cell, up_cell, 0.01),
      ploidy = seq(low_ploidy, up_ploidy, 0.1),
      chromosome.list = chr.fit,
      female = as.logical(is_female)
    )

    """
}

//Results

process SEQUENZA_RESULTS {

    script:
    """
    #!/usr/bin/env Rscript

    library(sequenza)
    """
}
