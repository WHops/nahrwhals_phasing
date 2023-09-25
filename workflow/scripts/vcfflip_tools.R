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
    region_id <- with(fix_df, paste(CHROM[i], gsub(" ", "", POS[i]), extract_end(INFO[i]), sep="-"))
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
      
      # Pre-compute common elements
      parts <- strsplit(genotype, ":", fixed=TRUE)[[1]]
      gt <- parts[1]
      custom_phase_block <- ifelse(length(parts) > 1, paste0(":", parts[2]), "")
      
      modified_gt = modify_genotype(gt, action, stats)
      stats[[action]] <- stats[[action]] + 1
      print(stats)
      # Update the genotype with any modifications, keeping the custom phase block
      vcf@gt[i, sample] <- paste0(modified_gt, custom_phase_block)
      
    }
  }
  
  vcfR::write.vcf(vcf, file=output_filename)
  return(stats)
}


#
# modify_genotype - Modify a genotype based on a specified action and update statistics.
#
# Args:
#   gt: A character vector representing the genotype to modify.
#   action: A character string representing the action to perform on the genotype.
#   stats: A named list representing the statistics to update.
#
# Returns:
#   A character vector representing the modified genotype.
#
# Example:
#   modify_genotype(c("0|1", "1|0"), "flip", list("flipped" = 0))
modify_genotype <- function(gt, action, stats) {
  
  switch(action,
         flip = {
           gt_parts <- strsplit(gt, "\\|", fixed=FALSE)[[1]]
           modified_gt <- paste(rev(gt_parts), collapse="|")
         },
         unphase = {
           gt_parts <- strsplit(gt, "\\|", fixed=FALSE)[[1]]
           sorted_gt_parts <- sort(gt_parts)
           modified_gt <- paste(sorted_gt_parts, collapse="/")
         },
         missing_unphase = {
           gt_parts <- strsplit(gt, "\\|", fixed=FALSE)[[1]]
           sorted_gt_parts <- sort(gt_parts)
           modified_gt <- paste(sorted_gt_parts, collapse="/")
         },
         keep = {
           modified_gt <- gt  # No changes to the genotype
         },
         modified_gt <- gt  # Default case, no changes to the genotype
  )
  
  return(modified_gt)
}


library(argparse)

parser <- ArgumentParser()

parser$add_argument("--vcf", help="VCF file to flip")
parser$add_argument("--actions", help="Actions file")
parser$add_argument("--vcfout", help="Output VCF file")
parser$add_argument("--logout", help="Output log file")

args <- parser$parse_args()

actions_df <- read.table(args$actions, header=FALSE, stringsAsFactors=FALSE, sep='\t')
actions <- process_actions(actions_df)
stats <- process_vcf(args$vcf, actions, args$vcfout)
print('hi')
print(stats)
writeLines(as.character(stats), args$logout)


# if (interactive()) {
#   flip_actions <- '~/cluster11/g/korbel/hoeps/projects/nahr/phaselab/snake_approach/res/actions_all_50000.tsv'
#   unphased_vcf <- '~/PhD/projects/nahrcall/phaselab/nahrwhals_phasing/data/all.vcf'
#   output_filename <- "~/Desktop/test.vcf.gz"
# }
