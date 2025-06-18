pms2_realignment <- function(args, tools, vardict, output.samples.dir){
  bams <- read.table(args$bamTxt)
  bam_header <- system_safe(paste(tools$samtools, "view -H", bams[1,1]), intern=TRUE)
  use_chr <- any(grepl("^@SQ.*SN:chr", bam_header))
  
  dir.create(output.samples.dir, showWarnings = FALSE)
  for (b in seq_len(nrow(bams))){
    #Get sample name
    sample.name <-bams[b,1] %>% as.character %>%
      basename() %>%
      stringr::str_replace(".bam", "")
    full <- bams[b,1] %>% as.character()
    
    #Create results dir per sample
    results.dir.sample <- file.path(output.samples.dir, sample.name)
    dir.create (results.dir.sample, showWarnings = FALSE)
    
    #Create folder to store first approach
    results.dir.1 <- file.path(results.dir.sample, "Approach1")
    results.dir.2 <- file.path(results.dir.sample, "Approach2")
    dir.create(results.dir.1)
    dir.create(results.dir.2)
    
    #Location for different files
    ## APPROACH 1
    vcfvardict <- file.path (results.dir.1, paste0(sample.name, "freq", vardict$freq, "_VARDICT.vcf"))
    filtered.bam <- paste0(results.dir.1, "/filtered", sample.name, ".bam" )
    r1 <- paste0(results.dir.1,"/r1_", sample.name, ".fastq")
    r2 <- paste0(results.dir.1,"/r2_", sample.name, ".fastq")
    sam <- paste0(results.dir.1, "/segon_alineament", sample.name, ".sam")
    second.bam<- paste0(results.dir.1, "/second_alignment_", sample.name, ".bam")
    sorted.bam <- file.path(results.dir.1, paste0("second_alignment_", sample.name, "_sorted"))
    
    
    ## APPROACH 2
    vcfvardict.E11 <- file.path (results.dir.2, paste0(sample.name, "freq", vardict$freq, "_VARDICT.vcf"))
    filtered.bam.E11 <- paste0(results.dir.2, "/filtrar", sample.name, ".bam" )
    filtered.bam.E11.1 <- paste0(results.dir.2, "/filtrar_resta_", sample.name, ".bam" )
    filtered.bam.E11.2 <- paste0(results.dir.2, "/filtrar_E11_", sample.name, ".bam" )
    filtered.bam.E11.3 <- paste0(results.dir.2, "/filtrar_E11_picard", sample.name, ".bam" )
    read.names <- paste0(results.dir.2, "/read_names", sample.name, ".txt" )
    r1.E11 <- paste0(results.dir.2,"/r1_", sample.name, ".fastq")
    r2.E11 <- paste0(results.dir.2,"/r2_", sample.name, ".fastq")
    sam.E11 <- paste0(results.dir.2, "/segon_alineament", sample.name, ".sam")
    second.bam.E11 <- paste0(results.dir.2, "/second_aligment_", sample.name, ".bam")
    third.bam.E11 <- paste0(results.dir.2, "/third_alignment", sample.name, ".bam")
    sorted.bam.E11 <- file.path(results.dir.2, paste0("third_alignment_", sample.name, "_sorted"))
    
    
    #PIPELINE
    
    #APROACH 1 ----
    
    ### 1. Filter gen + pseudogen regions with samtools
    chr_prefix <- if(use_chr) "chr" else ""
    position.filter <- ifelse(args$genome=="hg19",
                              paste0(chr_prefix, "7:6012830-6049024 ", chr_prefix, "7:6774736-6791432"),
                              paste0(chr_prefix, "7:5973199-5983191 ", chr_prefix, "7:5989600-6009393"))
    
    cmd1 <- paste(tools$samtools, "view -b -h", full, position.filter, ">", filtered.bam)
    print(cmd1); system_safe(cmd1)
    
    ### 2. Go back to fastq with picard
    cmd2 <- paste("java -jar", tools$picard, "SamToFastq",
                  "-I", filtered.bam,
                  "-F", r1,
                  "-F2", r2,
                  "-VALIDATION_STRINGENCY LENIENT")
    print(cmd2); system_safe(cmd2)
    
    ### 3. Realign with bwa mem
    cmd3 <- paste(tools$bwa,  "mem -t 1", args$reference, r1, r2, ">", sam)
    print(cmd3); system_safe(cmd3)
    
    ### 4. Convert sam to bam
    cmd4 <- paste(tools$samtools, "view -b -S", sam, ">", second.bam)
    print(cmd4); system_safe(cmd4)
    
    ### 5. Sort the bam file amb index the bam file
    cmd5 <- paste(tools$samtools, "sort -o", paste0(sorted.bam, ".bam"), second.bam)
    print(cmd5); system_safe(cmd5)
    cmd6 <- paste (tools$samtools , "index", paste0(sorted.bam, ".bam"))
    print(cmd6); system_safe(cmd6)
    
    #delete intermediate files
    system_safe(paste("rm", filtered.bam))
    system_safe(paste("rm", r1))
    system_safe(paste("rm", r2))
    system_safe(paste("rm", sam))
    system_safe(paste("rm", second.bam))
    
    # 6. Variant calling with VardictJava
    range.position <- ifelse(args$genome=="hg19", 
                             "chr7:6012350-6049257:PMS2", 
                             "chr7:5972719-6009626:PMS2")
    
    cat("Running VarDict command:\n")
    print(tools$vardict)
    
    print("VarDict parameters:")
    print(str(vardict))
    
    cmd7 <- paste(tools$vardict,
                  "-G", args$reference,
                  "-X", vardict$X,
                  "-q", vardict$phred_score,
                  "-m", vardict$missmatches,
                  "-f", vardict$freq,
                  "-N", sample.name,
                  "-b", paste0(sorted.bam, ".bam"),
                  "-R", range.position,
                  "|", "teststrandbias.R",
                  "|", "var2vcf_valid.pl",
                  "-N", sample.name,
                  "-E -A",
                  "-d 8",
                  "-q", vardict$phred_score,
                  "-f", vardict$freq,
                  ">", vcfvardict)
    
    print(cmd7); system_safe(cmd7)
    
    # APPROACH 2 ----
    
    ## 1. Filter regions:
    
    ### 1.1 Outside E11 -> Filter with samtools the regions of interst
    position.filter.outside.E11 <- ifelse(args$genome == "hg19",
                                          paste0(chr_prefix, "7:6012830-6022822 ", chr_prefix, "7:6029231-6049024 ", chr_prefix, "7:6774736-6775220 ", chr_prefix, "7:6781116-6791432"),
                                          paste0(chr_prefix, "7:5973199-5983191 ", chr_prefix, "7:5989600-6009393 ", chr_prefix, "7:6735105-6735589 ", chr_prefix, "7:6741485-6751801"))
    
    cmd1 <-paste(tools$samtools, "view -b -h", full,  position.filter.outside.E11, ">", filtered.bam.E11.1)
    print(cmd1); system_safe(cmd1)
    
    ### 1.2 Inside E11
    #### 1.2.1 Filter considering the invariable variants detailed in Gould et al
    invariable.positions <- ifelse(args$genome == "hg19",
                                   "chr7:6026364-6026364 chr7:6026598-6026598 chr7:6026601-6026601 chr7:6026625-6026625 chr7:6026636-6026636 chr7:6026707-6026707 chr7:6027017-6027017 chr7:6027474-6027474 >",
                                   "chr7:5986733-5986733 chr7:5986967-5986967 chr7:5986970-5986970 chr7:5986994-5986994 chr7:5987005-5987005 chr7:5987076-5987076 chr7:5987386-5987386 chr7:5987843-5987843 >")
    cmd1.2 <- paste(tools$samtools, "view -b -h", full, invariable.positions, filtered.bam.E11.2)
    print(cmd1.2); system_safe(cmd1.2)
    #### 1.2.2 Get reads name that overlap with this 7 positions (without duplicates)
    cmd1.3 <- paste( tools$samtools, "view", filtered.bam.E11.2 , " | cut -f1 | sort | uniq >", read.names)
    print(cmd1.3); system_safe(cmd1.3)
    
    #### 1.2.3 Filter reads from the original bam file taking into consideration the name (picard).
    cmd1.4 <- paste("java -jar", tools$picard, "FilterSamReads",
                    "-I", full,
                    "-O", filtered.bam.E11.3,
                    "-FILTER includeReadList",
                    "-READ_LIST_FILE", read.names,
                    "-VALIDATION_STRINGENCY LENIENT")
    print(cmd1.4); system_safe(cmd1.4)
    
    
    ## 2.  Go back to fastq with picard:
    ### 2.1 Outside E11
    cmd2 <- paste("java -jar", tools$picard, "SamToFastq",
                  "-I", filtered.bam.E11.1,
                  "-F", r1.E11,
                  "-F2", r2.E11,
                  "-UNPAIRED_FASTQ /dev/null",
                  "-INCLUDE_NON_PF_READS true",
                  "-VALIDATION_STRINGENCY LENIENT")
    print(cmd2); system_safe(cmd2)
    
    ### 2.2 Inside E11 -> not necessary
    
    ## 3. Realign with bwa mem
    ### 3.1 Outside E11
    cmd3 <- paste(tools$bwa,  "mem -t 1", args$reference, r1.E11, r2.E11, ">", sam.E11)
    print(cmd3); system_safe(cmd3)
    
    ### 3.2 Inside E11 -> not necessary
    
    ## 4. Convert from sam to bam
    ### 4.1 Outside E11
    cmd4 <- paste(tools$samtools, "view -b -S", sam.E11, ">", second.bam.E11)
    print(cmd4); system_safe(cmd4)
    
    ### 4.2 Inside E11 -> not necessary
    
    ## 5. Merge the two bams
    cmd4.1 <- paste("java -jar", tools$picard, "MergeSamFiles",
                    "-I", second.bam.E11,
                    "-I", filtered.bam.E11.3,
                    "-O", third.bam.E11,
                    "-VALIDATION_STRINGENCY LENIENT")
    print(cmd4.1); system_safe(cmd4.1)
    
    ## 6. Sort + index
    cmd5 <- paste(tools$samtools, "sort -o", paste0(sorted.bam.E11, ".bam"), third.bam.E11)
    print(cmd5); system_safe(cmd5)
    cmd6 <- paste (tools$samtools , "index", paste0(sorted.bam.E11, ".bam"))
    print(cmd6); system_safe(cmd6)
    
    #remove intermediate files
    system_safe(paste("rm", filtered.bam.E11.1))
    system_safe(paste("rm", filtered.bam.E11.2))
    system_safe(paste("rm", filtered.bam.E11.3))
    system_safe(paste("rm", r1.E11))
    system_safe(paste("rm", r2.E11))
    system_safe(paste("rm", sam.E11))
    system_safe(paste("rm", second.bam.E11))
    system_safe(paste("rm", third.bam.E11))
    
    ## 7. Variant calling with VardictJava
    
    cat("Running VarDict command:\n")
    print(tools$vardict)
    
    print("VarDict parameters:")
    print(str(vardict))
    
    cmd7 <- paste(tools$vardict,
                  "-G", args$reference,
                  "-X", vardict$X,
                  "-q", vardict$phred_score,
                  "-m", vardict$missmatches,
                  "-f", vardict$freq,
                  "-N", sample.name,
                  "-b", paste0(sorted.bam.E11, ".bam"),
                  "-R", range.position,
                  "|", "teststrandbias.R",
                  "|", "var2vcf_valid.pl",
                  "-N", sample.name,
                  "-E -A",
                  "-d 8",
                  "-q", vardict$phred_score,
                  "-f", vardict$freq,
                  ">", vcfvardict.E11) 
    print(cmd7); system_safe(cmd7)
  }
}


zip_tabix <- function(resultsDir, patro){
  directoris <- list.dirs(resultsDir, recursive = TRUE)
  lapply(directoris, function(x){
    message(paste("Tabix index created for", x))
    file <- list.files(x, pattern=patro, recursive=FALSE)
    message(length(file))
    if (length(file)>0){
      f <- file.path(x, file)
      
      # First, check the chromosome format in the VCF
      vcf_content <- readLines(f, n=100)
      header_line <- grep("^#CHROM", vcf_content, value=TRUE)
      first_variant <- grep("^[^#]", vcf_content, value=TRUE)[1]
      
      message("Header line: ", header_line)
      message("First variant: ", first_variant)
      
      cmd9 <- paste("bgzip -f", f)
      print(cmd9); system_safe(cmd9)
      
      # Create tabix index with appropriate sequence naming
      cmd10 <- paste("tabix -f", paste0(f, ".gz"))
      print(cmd10); system_safe(cmd10)
    }
  })
}

convert_vcf_txt <- function(resultsDir, vcf.pattern, rng, args){
  files <- list.files(resultsDir, pattern=vcf.pattern, recursive=TRUE, full.names=TRUE)
  # bams <- read.table(args$bamTxt)$V1
  # samples.name <- bams[,1] %>%
  #   as.character() %>%
  #   basename() %>%
  #   stringr::str_replace(".bam", "") %>%
  #   rep(each=2)
  samples.name <- basename(files) %>% stringr::str_replace("freq.+", "")
  seqlevels(rng) <- "chr7"
  seqlevels(bedGR) <- "chr7"
  for (f in seq_len(length(files))){ 
    sample.name <- samples.name[f]
    print(sample.name)
    print(f)
    tab <- Rsamtools::TabixFile(files[f])
    # Read vcf
    vcf <- vcfR::read.vcfR(tab$path)
    # Ensure chromosome format matches
    vcf.rng <- VariantAnnotation::readVcf(tab, args$genome, param=rng)
    vcf.rng.exons <- VariantAnnotation::readVcf(tab, args$genome, param=bedGR)
    exons.df <- data.frame(ID= names(vcf.rng.exons@rowRanges),
                           exon = vcf.rng.exons@rowRanges$paramRangeID %>% as.character()) %>%
      dplyr::group_by(ID) %>% dplyr::slice (1)
    
    VariantAnnotation::writeVcf(vcf.rng, file.path(dirname(files[f]), "rng_vardict.vcf"))
    vcf.rng <- vcfR::read.vcfR(file.path(dirname(files[f]), "rng_vardict.vcf"))
    myID <- getID(vcf.rng)
    length(unique(myID, incomparables = NA)) == length(myID)
    vcf.rng <- vcf.rng[!duplicated(myID, incomparables = NA), ]
    #tidy vcfs and make a dataframe
    #a) all vcf
    vcf.tidy <- vcfR::vcfR2tidy(vcf, single_frame = TRUE)
    vcf.df <- as.data.frame(vcf.tidy$dat) %>%
      dplyr::mutate(ID= paste0(CHROM, ":", POS, "_", REF, "/", ALT))
    #b) in roi vcf
    vcf.tidy.rng <- vcfR::vcfR2tidy(vcf.rng, single_frame = TRUE)
    vcf.df.rng <- as.data.frame(vcf.tidy.rng$dat)
    
    
    
    #Merge two dataframes to mark in ROIS
    sel <- vcf.df$POS %in% vcf.df.rng$POS
    all.variants <- merge(vcf.df, exons.df, by="ID", all=TRUE) %>%
      mutate(INROI= sel,
             ID_sample=all_of(sample.name)) %>%
      relocate(CHROM, POS, ID, REF, ALT,INROI, exon)
    
    write.table(all.variants, file.path(dirname(files[f]), "all.variants.txt") , sep="\t", row.names = F)
  }
}

merge_pipelines <- function(resultsDir, vars.paralogous, classification, args) {
  samples <- list.dirs(resultsDir, recursive = FALSE)
  print(paste("Processing directory:", resultsDir))
  print(paste("Found", length(samples), "samples"))
  
  final.results <- file.path(resultsDir, "final_results")
  dir.create(final.results, showWarnings = FALSE)
  
  for(sample_dir in samples) {
    message(paste("Processing sample:", basename(sample_dir)))
    
    all.vars <- list.files(sample_dir, 
                           pattern = "all.variants.txt", 
                           recursive = TRUE, 
                           full.names = TRUE)
    
    if(length(all.vars) == 2) {
      tryCatch({
        # Read files and ensure column names are correct
        file1 <- read.delim(all.vars[1], sep = "\t", stringsAsFactors = FALSE)
        file2 <- read.delim(all.vars[2], sep = "\t", stringsAsFactors = FALSE)
        
        # Rename AF columns before merge to avoid confusion
        if("AF" %in% colnames(file1)) file1$AF.x <- file1$AF
        if("AF" %in% colnames(file2)) file2$AF.y <- file2$AF
        
        # Select columns for merge
        merge_cols <- c("ID", "CHROM", "POS", "REF", "ALT", "INROI", 
                        "exon", "FILTER", "SAMPLE", "TYPE", "AF.y")
        file2_subset <- file2[, merge_cols[merge_cols %in% colnames(file2)]]
        
        # Merge files
        both <- merge(file1, file2_subset, by = "ID", all = TRUE, 
                      suffixes = c(".x", ".y"))
        
        # Clean up merged data
        both <- both %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            CHROM = ifelse(is.na(CHROM.x), CHROM.y, CHROM.x),
            POS = ifelse(is.na(POS.x), POS.y, POS.x),
            REF = ifelse(is.na(REF.x), REF.y, REF.x),
            ALT = ifelse(is.na(ALT.x), ALT.y, ALT.x),
            TYPE = ifelse(is.na(TYPE.x), TYPE.y, TYPE.x),
            FILTER = ifelse(is.na(FILTER.x), FILTER.y, FILTER.x),
            exon = ifelse(is.na(exon.x), exon.y, exon.x),
            INROI = ifelse(is.na(INROI.x), INROI.y, INROI.x),
            # Handle AF more comprehensively
            AF = case_when(
              !is.na(AF.x) ~ as.numeric(AF.x),
              !is.na(AF.y) ~ as.numeric(AF.y),
              TRUE ~ NA_real_
            ),
            AF.x = as.numeric(AF.x),
            AF.y = as.numeric(AF.y)
          )
        
        # Process variants
        message("Running variant annotation...")
        all.variants.final <- cdnav2(dataset = both, 
                                     vars.paralogous = vars.paralogous, 
                                     classification = classification, 
                                     args = args)
        
        # NEW CODE HERE - Add required columns and reorder
        all.variants.final <- all.variants.final %>%
          dplyr::mutate(
            NM = "NM_000535.7",
            NP = "NP_000526.2"
          ) %>%
          dplyr::select(cdna, prot, NM, NP, class_paralogous, 
                        paralogous_above60, present_pipelines, 
                        everything())
        
        # Add QC metrics
        all.variants.final <- all.variants.final %>%
          dplyr::mutate(
            quality_score = case_when(
              as.numeric(MQ) >= 60 & as.numeric(DP) >= 20 ~ "HIGH",
              as.numeric(MQ) >= 40 & as.numeric(DP) >= 10 ~ "MEDIUM",
              TRUE ~ "LOW"
            ),
            strand_bias = ifelse(!is.na(SBF) & as.numeric(SBF) < 0.01, 
                                 "HIGH", "LOW")
          )
        
        # Apply filtering
        all.variants.final <- filter_variants(all.variants.final)
        
        # Write output
        output_file <- file.path(final.results, 
                                 paste0(basename(sample_dir), 
                                        "_PMS2_variants.txt"))
        write.table(all.variants.final, 
                    output_file, 
                    sep = "\t", 
                    row.names = FALSE, 
                    col.names = TRUE, 
                    quote = FALSE)
        
      }, error = function(e) {
        message(paste("Error processing sample:", basename(sample_dir)))
        message(e)
        # Print additional debugging info
        message("\nData structure:")
        if(exists("file1")) print(str(file1))
        if(exists("file2")) print(str(file2))
      })
    } else {
      warning(paste("Expected 2 variant files for sample", 
                    basename(sample_dir), 
                    "but found", 
                    length(all.vars)))
    }
  }
}

# Helper function to generate QC report
generate_qc_report <- function(variants) {
  qc_stats <- data.frame(
    metric = c(
      "Total variants",
      "Variants in ROI",
      "PASS filter variants",
      "High quality variants",
      "Variants with strand bias",
      "Paralogous variants",
      "Variants needing LR-PCR",
      "Mean depth",
      "Mean mapping quality"
    ),
    value = c(
      nrow(variants),
      sum(variants$INROI == TRUE, na.rm = TRUE),
      sum(variants$FILTER == "PASS", na.rm = TRUE),
      sum(variants$quality_score == "HIGH", na.rm = TRUE),
      sum(variants$strand_bias == "HIGH", na.rm = TRUE),
      sum(variants$class_paralogous == "paralogous", na.rm = TRUE),
      sum(variants$what_to_do == "Perform LR-PCR", na.rm = TRUE),
      mean(as.numeric(variants$DP), na.rm = TRUE),
      mean(as.numeric(variants$MQ), na.rm = TRUE)
    )
  )
  return(qc_stats)
}

cdnav2 <- function(dataset, vars.paralogous, classification, args) {
  NC <- ifelse(args$genome == "hg19",
               "NC_000007.13",
               "NC_000007.14")
  
  # Convert AF values properly
  result <- dataset %>%
    # First remove chr prefix from ID and CHROM
    dplyr::mutate(across(starts_with("CHROM"), ~gsub("chr", "", .))) %>%
    # Handle AF values
    dplyr::mutate(
      AF.x = as.numeric(gsub("\"", "", AF.x)),
      AF.y = as.numeric(gsub("\"", "", AF.y)),
      AF = as.numeric(gsub("\"", "", AF))
    ) %>%
    # Clean variant format
    dplyr::mutate(
      variant = case_when(
        str_length(REF) == 1 ~ paste0(NC, "%3Ag.", POS, REF, ">", ALT),
        TRUE ~ paste0(NC, "%3Ag.", as.numeric(POS+1), "_", 
                      as.numeric(POS+str_length(REF)-1), "del", 
                      str_sub(REF, 2, -1))
      )
    ) %>%
    # Handle insertions and complex variants
    dplyr::mutate(
      variant = case_when(
        str_length(ALT) > 1 & str_length(REF) == 1 ~ 
          paste0(NC, "%3Ag.", as.numeric(POS), "_", 
                 as.numeric(POS+1), "ins", str_sub(ALT, 2, -1)),
        str_length(ALT) > 1 & str_length(REF) > 1 ~
          paste0(NC, "%3Ag.", as.numeric(POS), "_",
                 as.numeric(POS+str_length(REF)-1), "delins", ALT),
        TRUE ~ variant
      )
    )
  
  # Initialize cdna and prot columns
  result$cdna <- ":"
  result$prot <- ":"
  
  # Process annotations through Mutalyzer
  server_mutalyzerv3 <- "https://v3.mutalyzer.nl/api/normalize/"
  
  for(i in seq_len(nrow(result))) {
    if(result$CHROM[i] == "7") {
      if(!(result$TYPE[i] %in% c("INV", "<DEL>", "<DUP>", "DEL", "DUP")) && 
         str_length(result$variant[i]) < 100) {
        
        mut <- NULL
        attempts <- 0
        while(is.null(mut) && attempts < 3) {
          mut <- tryCatch({
            jsonlite::fromJSON(paste0(server_mutalyzerv3, 
                                      URLencode(result$variant[i])))
          }, error = function(e) NULL)
          attempts <- attempts + 1
          if(is.null(mut)) Sys.sleep(1)  # Wait before retry
        }
        
        if(!is.null(mut) && !is.null(mut$equivalent_descriptions)) {
          nm_entries <- grep("NM_000535.7", 
                             unlist(mut$equivalent_descriptions), 
                             value = TRUE)
          if(length(nm_entries) > 0) {
            result$cdna[i] <- nm_entries[1]
            result$prot[i] <- if(length(nm_entries) > 1) nm_entries[2] else ":"
          }
        }
      }
    }
  }
  
  # Get variant classification lists based on genome version
  if(args$genome == "hg19") {
    vars.par <- vars.paralogous$ID.hg19
    vars.ben <- classification$ben$ID.vars
    vars.vus <- classification$vus$ID.vars
    vars.pat <- classification$pat$ID.vars
    vars.lpat <- classification$lpat$ID.vars
    vars.lben <- classification$lben$ID.vars
  } else {
    vars.par <- vars.paralogous$ID.hg38
    vars.ben <- classification$ben$ID.vars.hg38
    vars.vus <- classification$vus$ID.vars.hg38
    vars.pat <- classification$pat$ID.vars.hg38
    vars.lpat <- classification$lpat$ID.vars.hg38
    vars.lben <- classification$lben$ID.vars.hg38
  }
  
  # Final classification and organization
  result <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Check against both hg19 and hg38 IDs for paralogous variants
      class_paralogous = ifelse(ID %in% c(vars.paralogous$ID.hg19, vars.paralogous$ID.hg38), 
                                "paralogous", ""),
      class = case_when(
        ID %in% vars.ben ~ "BEN",
        ID %in% vars.vus ~ "VUS",
        ID %in% vars.pat ~ "PAT",
        ID %in% vars.lpat ~ "lPAT",
        ID %in% vars.lben ~ "lBEN",
        TRUE ~ "Not found"
      ),
      present_pipelines = case_when(
        is.na(AF.y) ~ "general",
        is.na(AF.x) ~ "E11",
        TRUE ~ "2 pipelines"
      ),
      paralogous_above60 = !is.na(AF.x) && class_paralogous == "paralogous" && AF.x >= 0.6
    )
  
  return(result)
}


filter_variants <- function(all.variants.final){
  cat("Starting filter_variants function\n")
  cat("Initial number of variants:", nrow(all.variants.final), "\n")
  
  # Print quality metrics before filtering
  cat("\nQuality metrics distribution:\n")
  cat("Mean mapping quality:", mean(as.numeric(all.variants.final$MQ), na.rm=TRUE), "\n")
  cat("Mean depth:", mean(as.numeric(all.variants.final$DP), na.rm=TRUE), "\n")
  cat("Mean allele frequency:", mean(as.numeric(all.variants.final$AF), na.rm=TRUE), "\n")
  
  cat("\nFilter types:\n")
  print(table(all.variants.final$FILTER.x, useNA="ifany"))
  
  # Check for any potential pathogenic variants that might be getting filtered
  potential_pathogenic <- all.variants.final %>%
    dplyr::filter(class %in% c("PAT", "lPAT") | grepl("splice|nonsense|frameshift", cdna, ignore.case=TRUE))
  
  if(nrow(potential_pathogenic) > 0) {
    cat("\n!!! POTENTIAL PATHOGENIC VARIANTS DETECTED !!!\n")
    # Fix: Ungroup before select to avoid the rowwise_df error
    print(potential_pathogenic %>% dplyr::ungroup() %>% dplyr::select(ID, cdna, prot, class, FILTER.x, AF))
  }
  
  all.variants.final2 <- all.variants.final %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # First determine basic filtering status
      filter_status = case_when(
        is.na(FILTER.x) ~ "FAILED",
        !FILTER.x %in% c("PASS", "NM5.25", "InIns", "NM5.25;InIns", "InIns;NM5.25") ~ "FAILED",
        INROI == FALSE ~ "OUT_OF_ROI",
        TRUE ~ "PASSED"
      ),
      
      # Determine variant classification status
      classification_status = case_when(
        filter_status != "PASSED" ~ filter_status,
        class %in% c("BEN", "lBEN") ~ "BENIGN",
        class %in% c("PAT", "lPAT") ~ "PATHOGENIC",
        class == "VUS" ~ "VUS",
        TRUE ~ "UNKNOWN"
      ),
      
      # Handle paralogous variants
      paralogous_status = case_when(
        filter_status != "PASSED" ~ filter_status,
        class_paralogous == "paralogous" & !paralogous_above60 ~ "PARALOGOUS_LOW_AF",
        class_paralogous == "paralogous" & paralogous_above60 ~ "PARALOGOUS_HIGH_AF",
        TRUE ~ "NON_PARALOGOUS"
      ),
      
      # Final decision logic with lower threshold for PATHOGENIC variants
      what_to_do = case_when(
        # Always recommend LR-PCR for pathogenic variants, even if quality filters failed
        class %in% c("PAT", "lPAT") & present_pipelines == "2 pipelines" ~ "Perform LR-PCR",
        class %in% c("PAT", "lPAT") & present_pipelines != "2 pipelines" ~ "Only perform LR-PCR if IHC PMS2-",
        # Check for potential pathogenic variants based on transcript consequences
        grepl("splice|nonsense|frameshift", cdna, ignore.case=TRUE) & INROI == TRUE ~ "Perform LR-PCR",
        # Regular filtering logic
        filter_status == "FAILED" ~ "QUALITY FILTER NOT PASSED",
        filter_status == "OUT_OF_ROI" ~ "Out of ROI",
        paralogous_status == "PARALOGOUS_LOW_AF" ~ "Do not classify, paralogous variant <60",
        classification_status == "BENIGN" ~ "Do not report, BEN/lBEN variant",
        classification_status == "VUS" & present_pipelines == "2 pipelines" ~ "Report VUS",
        classification_status == "VUS" & present_pipelines != "2 pipelines" ~ "Do not report VUS",
        TRUE ~ "Classify variant"
      )
    )
  
  # Add debug counts
  cat("\nFilter status distribution:\n")
  print(table(all.variants.final2$filter_status))
  
  cat("\nClassification status distribution:\n")
  print(table(all.variants.final2$classification_status))
  
  cat("\nParalogous status distribution:\n")
  print(table(all.variants.final2$paralogous_status))
  
  cat("\nFinal what_to_do distribution:\n")
  print(table(all.variants.final2$what_to_do))
  
  return(all.variants.final2)
}

system_safe <- function(cmd, intern = FALSE) {
  res <- system(cmd, intern = intern)
  if (!intern && res != 0) {
    stop(sprintf("Command failed with exit code %d: %s", res, cmd))
  }
  return(res)
}