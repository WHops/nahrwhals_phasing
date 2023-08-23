# Convert actions dataframe into a more efficient lookup structure
process_actions <- function(actions_df) {
  action_list <- list()
  
  for(i in seq(1, nrow(actions_df), by=2)) {
    region <- actions_df$V1[i]
    sample <- actions_df$V2[i]
    h1_action <- actions_df$V4[i]
    h2_action <- actions_df$V4[i+1]
    
    final_action <- ifelse(h1_action == h2_action, h1_action, 'unphase')
    action_list[[region]][[sample]] <- final_action
  }
  
  return(action_list)
}

# Extract the END from the INFO field
extract_end <- function(info_str) {
  as.numeric(stringr::str_extract(info_str, "(?<=END=)\\d+"))
}

process_vcf <- function(vcf_filename, action_list, output_filename) {
  vcf <- vcfR::read.vcfR(vcf_filename)
  fix_df <- as.data.frame(vcf@fix)
  
  stats <- list(keep=0, flip=0, unphase=0, missing_unphase=0)
  
  for(i in seq_along(fix_df$ID)) {
    region_id <- with(fix_df, paste(CHROM[i], as.numeric(POS[i]), extract_end(INFO[i]), sep="-"))
    genotypes <- setdiff(colnames(vcf@gt), "FORMAT")
    missing_idx <- which(is.na(vcf@gt[i, genotypes]))
    vcf@gt[i, missing_idx] <- './.'
    
    for(sample in genotypes) {
      genotype <- vcf@gt[i, sample]
      
      # Check if region and sample are in the action_list
      region_exists <- region_id %in% names(action_list)
      sample_exists <- region_exists && (sample %in% names(action_list[[region_id]]))
      
      # Assign action based on existence checks
      if(region_exists && sample_exists) {
        action <- action_list[[region_id]][[sample]]
      } else {
        action <- "missing_unphase"
      }
      
      print(paste(region_id, sample, action, genotype))
      
      switch(action,
             flip = {
               stats$flip <- stats$flip + 1
               vcf@gt[i, sample] <- ifelse(genotype == "0|1", "1|0", ifelse(genotype == "1|0", "0|1", genotype))
             },
             unphase = {
               stats$unphase <- stats$unphase + 1
               vcf@gt[i, sample] <- ifelse(genotype %in% c("0|1", "1|0"), "0/1", genotype)
             },
             missing_unphase = {
               stats$missing_unphase <- stats$missing_unphase + 1
               vcf@gt[i, sample] <- ifelse(genotype %in% c("0|1", "1|0"), "0/1", genotype)
             },
             keep = {
               stats$keep <- stats$keep + 1
             }
      )
    }
  }
  
  vcfR::write.vcf(vcf, file=output_filename)
  return(stats)
}

flip_actions <- '~/cluster10/g/korbel/hoeps/projects/nahr/phaselab/snake_approach/res/actions_all.tsv'
unphased_vcf <- '~/Desktop/all2.vcf'
output_filename <- "~/Desktop/test.vcf.gz"

actions_df <- read.table(flip_actions, header=FALSE, stringsAsFactors=FALSE)
actions <- process_actions(actions_df)
stats <- process_vcf(unphased_vcf, actions, output_filename)
print(stats)
