#' Extracts the name and length of each sequence in a fasta file.
#'
#' This function takes an input fasta file and an output file as arguments,
#' and runs a 'gawk' command to extract the name and length of each sequence
#' in the fasta file. The extracted information is saved to the output file.
#'
#' @param inputfa A character string specifying the path to the input fasta file.
#' @param outputinfo A character string specifying the path to the output file.
#'
#' @return Nothing is returned by the function, but the sequence names and lengths are saved to the output file.
#' @author Wolfram Hoeps
#' @export
get_fasta_info <- function(inputfa, outputinfo) {
  # Run the gawk command to get the fasta file info
  cmd <- paste0(
    "gawk '/^>/{if (l!=\"\") print l; print; l=0; next}{l+=length($0)}END{print l}' ",
    inputfa, " > ", outputinfo
  )
  system(cmd)
}

#' Chunkify query fasta
#'
#' @description This is a helperfunction calling bedtools to
#' chop a query sequence into chunks.
#'
#' @param infasta A link to a single-seq fasta to be chopped
#' @param outfasta_chunk A link to the output path for a chopped multi-seq fasta.
#' @param chunklen length of sequence chunks in bp
#' @param params a list with all NAHRwhals parameters
#' @return nothing. But output files written.
#'
#' @author Wolfram Hoeps
#' @export
shred_seq_bedtools <- function(infasta,
                               outfasta_chunk,
                               chunklen,
                               bedtools_bin) {

  chunklen = as.numeric(chunklen)
  # Write a temporary bedfile that will be removed at the end of the function
  bed_tmp_file <- paste0("tmpbed_deleteme_", sprintf("%.0f", runif(1, 1e13, 1e14)), ".bed")
  fasta_awk_tmp_file <- paste0("tmpinfo_deleteme_", sprintf("%.0f", runif(1, 1e13, 1e14)), ".bed")

  get_fasta_info(infasta, fasta_awk_tmp_file)
  # Collect information about the fasta we want to chop. This info is needed for bedtools getfasta
  # To work well.
  name_len_df <- read.table(fasta_awk_tmp_file)
  contigname <- sub(">", "", (name_len_df)[1, ])
  contiglen <- as.numeric((name_len_df)[2, ])

  bed_df <- data.frame(
    seqnames = contigname,
    start = sprintf("%d", seq(0, contiglen - (contiglen %% chunklen), by = chunklen)),
    end = sprintf("%d", pmin(seq(0, contiglen - (contiglen %% chunklen), by = chunklen) + (chunklen - 1), contiglen))
  )

  write.table(bed_df, file = bed_tmp_file, sep = "\t", quote = F, row.names = F, col.names = F)


  system(paste0("rm ", infasta, ".fai"))

  sedcmd <- "sed -r \'s/(.*):/\\1_/'"
  system(paste0(bedtools_bin, " getfasta -fi ", infasta, " -bed ", bed_tmp_file, " | ", sedcmd, " > ", outfasta_chunk))

  system(paste0("rm ", bed_tmp_file))
  system(paste0("rm ", fasta_awk_tmp_file))

  # Make a print statement informing the user about the files that have been written, and the chunklen used.
  print(paste0("Wrote ", outfasta_chunk, " with chunklen ", chunklen, " bp."))
}



#' @title: aln_chunks_to_minimap
#' @export
get_fasta_and_shred <- function(sample, hap, region, chunklen, res_path, out_fasta, bedtools_bin){
    # Find the correct chunked reads fasta
  dir_path = paste0(res_path, region, '/fasta')


  #Criteria: filename has to match sample, hapx or hx and end on y.fa
  regex <- paste0("^.*", sample, ".*(h|hap)", hap, ".*y\\.fa$")
  
  # And then we grep those out of the whole list. 
  asm_fasta <- grep(regex, list.files(path = dir_path, full.names = TRUE), value = TRUE)
  #asm_chunked_fasta = paste0(asm_fasta, '_chunked_phasing.fa')

  # Use shred_seq_bedtools to turn this into chunks.
  shred_seq_bedtools(asm_fasta, out_fasta, chunklen, bedtools_bin)
}



#' @title determine_phase_with_whatshap
#' @export
subset_vcf_to_singlesample <- function(input_vcf_allsamples, sample, 
                                       output_vcf_singlesample_vcf, output_vcf_singlesample_vcf_gz, 
                                       bgzip_bin, bcftools_bin){

  
  # Some not so pleasant renamings. The age-old problem with GM vs NA names...

  vcf_samplename = sample
  if (sample == 'HG002'){
    vcf_samplename = 'NA24385'
  } else if(startsWith(sample, "GM")){
    vcf_samplename <- sub("^GM", "NA", sample)
  }
  
  bcftools_cmd = paste0(bcftools_bin, ' view -s ', vcf_samplename, ' ', input_vcf_allsamples, ' > ', output_vcf_singlesample_vcf)

  # bgzip system command that compresses output_vcf_singlesample_vcf and writes the file to output_vcf_singlesample_vcf_gz
  bgzip_cmd_2 = paste0(bgzip_bin , ' -c ',  output_vcf_singlesample_vcf, ' > ', output_vcf_singlesample_vcf_gz)

  system(bcftools_cmd)
  system(bgzip_cmd_2)
  Sys.sleep(1)
}


# @title: collect_whatshap_res
# @export
collect_whatshap_res <- function(haptags, sample, region, hap, summarylist_link){
     awk_cmd = paste0('tail -n +2 ', haptags, ' | awk \'BEGIN {FS=OFS="\t"} {print "', sample, '\t', hap, '\t', region, '",$1,$2}\' > ', summarylist_link)
     print(awk_cmd)
     system(awk_cmd)
     print('done')
     Sys.sleep(1)
  
}

# @title: evaluate_åsummarylist
# @export
evaluate_summarylist <- function(summarylist, actionlist){

  summ = read.table(summarylist, sep='\t')

  asm_hap = summ[1,'V2']
  sample = summ[1,'V1']
  region = summ[1,'V3']

  xx = as.data.frame(t(as.data.frame.character(table(summ$V5))))

  print(xx)
  print('##################')
  if (is.null(xx$H1)){
    xx$H1 = 0
  }
  if (is.null(xx$H2)){
    xx$H2 = 0
  }
  if (is.null(xx$none)){
    xx$none = 0
  }

  # Dirty bug killing:
  if ((xx$H1 + xx$H2) == 0){
  xx$H1 = 1
  xx$H2 = 1
  }


  h1_to_h2_fract =  xx$H1 / (xx$H1 + xx$H2)
  if (is.na(h1_to_h2_fract)){
    h1_to_h2_fract = 0.5
  }
  none_fract = xx$none / (xx$H1 + xx$H2 + xx$none)

  # Determine mapped haplotype
  if (h1_to_h2_fract > 0.5){
    mapped_hap = '1'
  } else if (h1_to_h2_fract < 0.5){
    mapped_hap = '2'
  } else {
    mapped_hap = 'unclear'
  }

  # Determine action.
  if (asm_hap == mapped_hap){
    action = 'keep'
  } else if ((asm_hap == 1) & (mapped_hap == 2) | (asm_hap == 2) & (mapped_hap == 1)){
    action = 'flip'
  } else if (mapped_hap == 'unclear'){
    action = 'unclear'
  }

  # Columns of the final_actionlist: region, sample, asm_hap, action, reads_H1, reads_H2, reads_unphased, h1_to_h2_fract, none_fract
  final_actionlist = paste0(paste(region, 
                                  sample,
                                  asm_hap, 
                                  action, 
                                  as.integer(as.numeric(xx$H1)), 
                                  as.integer(as.numeric(xx$H2)), 
                                  as.integer(as.numeric(xx$none)), 
                                  round(h1_to_h2_fract, 3), 
                                  round(none_fract,3), 
                                  sep='\t'), 
                           '\n')

  cat(final_actionlist, file=actionlist, append = TRUE)
}





library(argparse)

parser <- ArgumentParser()

parser$add_argument("--function_name")
parser$add_argument("--res_path")
parser$add_argument("--region")
parser$add_argument("--sample")
parser$add_argument("--hap")
parser$add_argument("--chunklen")
parser$add_argument("--aln_bam")
parser$add_argument("--haptags")
parser$add_argument("--summarylist_link")
parser$add_argument("--summarylist")
parser$add_argument("--actionlist")
parser$add_argument("--hg38_mmi")
parser$add_argument("--hg38_fa")
parser$add_argument("--vcf_dir")
parser$add_argument("--mm2_bin")
parser$add_argument("--samtools_bin")
parser$add_argument("--whatshap_bin")
parser$add_argument("--bedtools_bin")
parser$add_argument("--subset_vcf_allsamples")
parser$add_argument("--tabix_bin")
parser$add_argument("--bgzip_bin")
parser$add_argument("--bcftools_bin")
parser$add_argument("--subset_vcf_singlesample_vcf")
parser$add_argument("--subset_vcf_singlesample_vcf_gz")
parser$add_argument("--subset_vcf_singlesample_vcf_tbi")
parser$add_argument("--subset_vcf_singlesample_vcf_gz_tbi")
parser$add_argument("--asm_chunked_fasta")




args <- parser$parse_args()


if (args$function_name == "collect_whatshap_res") {
  collect_whatshap_res(args$haptags, args$sample, args$region, args$hap, args$summarylist_link)
} else if (args$function_name == "evaluate_summarylist") {
  evaluate_summarylist(args$summarylist, args$actionlist)
} else if (args$function_name == "subset_vcf_to_singlesample") {
  subset_vcf_to_singlesample(args$subset_vcf_allsamples, args$sample, args$subset_vcf_singlesample_vcf, args$subset_vcf_singlesample_vcf_gz, args$bgzip_bin, args$bcftools_bin)
} else if (args$function_name == "get_fasta_and_shred") {
  get_fasta_and_shred(args$sample, args$hap, args$region, args$chunklen, args$res_path, args$asm_chunked_fasta, args$bedtools_bin)
} else {
  print("No function name given.")
}

