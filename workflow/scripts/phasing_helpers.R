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






#' @title aln_chunks_to_minimap
#' @description Aligns chunks of the assembly to the reference genome using minimap2
#' @param res_path Path to the results folder
#' @param region Region to be aligned
#' @param sample Sample to be aligned
#' @param hap Haplotype to be aligned
#' @param hg38_mmi Path to the minimap2 index of the reference genome
#' @param minimap2_bin Path to the minimap2 binary
#' @param samtools_bin Path to the samtools binary
#' @return Path to the bam file
#' @export
#' @author Wolfram Hoeps
aln_chunks_to_minimap <- function(res_path, region, sample, hap,
                                  hg38_mmi, minimap2_bin, samtools_bin,
                                  bedtools_bin, chunklen = 10000) {
  # Find the correct chunked reads fasta
  dir_path = paste0(res_path, region, '/fasta')


  #Criteria: filename has to match sample, hapx or hx and end on y.fa
  regex <- paste0("^.*", sample, ".*(h|hap)", hap, ".*y\\.fa$")
  
  # And then we grep those out of the whole list. 
  asm_fasta <- grep(regex, list.files(path = dir_path, full.names = TRUE), value = TRUE)
  asm_chunked_fasta = paste0(asm_fasta, '_chunked_phasing.fa')

  # Use shred_seq_bedtools to turn this into chunks.
  shred_seq_bedtools(asm_fasta, asm_chunked_fasta, chunklen, bedtools_bin)

  # Make some assert statements to make sure that asm_chunked_fasta has (1) been written within the last second and (2) is roughly the same size as asm_fasta 
  assert_that(file.exists(asm_chunked_fasta))
  assert_that(file.info(asm_chunked_fasta)$mtime > Sys.time() - 1)
  assert_that(file.info(asm_chunked_fasta)$size > file.info(asm_fasta)$size * 0.9)
  assert_that(file.info(asm_chunked_fasta)$size < file.info(asm_fasta)$size * 1.1)
  
  outfile_sam = paste0('/scratch/hoeps/bamsam/', sample, '_h', hap, '_', region, '.sam')
  outfile_bam_unsrt = paste0('/scratch/hoeps/bamsam/', sample, '_h', hap, '_', region, '_unsrt.bam')
  outfile_bam = paste0('/scratch/hoeps/bamsam/', sample, '_h', hap, '_', region, '.bam')

  minimap2_cmd = paste0(minimap2_bin, " -a ", hg38_mmi, " ", 
                        asm_chunked_fasta, " > ",  outfile_sam)
  
  system(minimap2_cmd)
  
  system(paste0(samtools_bin, ' view -h -b ', outfile_sam, ' > ', outfile_bam_unsrt ))
  system(paste0(samtools_bin, ' sort ', outfile_bam_unsrt , ' > ', outfile_bam))

  system(paste0(samtools_bin, ' index ', outfile_bam ))

  
  return(outfile_bam)
}


determine_phase_with_whatshap <- function(aln_bam, region, sample, hap, hg38_fa, 
                                          vcf_dir, subset_vcf_allsamples, whatshap_bin){
  
  
  # Some not so pleasant renamings. The age-old problem with GM vs NA names...

  vcf_samplename = sample
  if (sample == 'HG002'){
    vcf_samplename = 'NA24385'
  } else if(startsWith(sample, "GM")){
    vcf_samplename <- sub("^GM", "NA", sample)
  }

  # Lets be brave
  tabix_bin = 'tabix'
  bcftools_bin = 'bcftools'

  chr = strsplit(region, '-')[[1]][1]
  start = strsplit(region, '-')[[1]][2]
  end = strsplit(region, '-')[[1]][3]
  phased_vcf_regex <- paste0("^.*", chr,"\\..*\\.vcf\\.gz$")
  phase_vcf_path = grep(phased_vcf_regex, list.files(path = vcf_dir, full.names = TRUE), value = TRUE)
  
  # Prepare vcf
  # subset_vcf_allsamples = paste0('/scratch/hoeps/nygc_subsets/', region, '_allsamples.vcf')
  subset_vcf_singlesample = paste0('/scratch/hoeps/nygc_subsets/', region, '_', sample, '_', hap, '.vcf')

  # tabix_cmd = paste0(tabix_bin, ' ', phase_vcf_path, ' ', chr, ':', start, '-', end, ' -h > ', subset_vcf_allsamples)
  # bgzip_cmd_1 = paste0('bgzip ', subset_vcf_allsamples)

  bcftools_cmd = paste0(bcftools_bin, ' view -s ', vcf_samplename, ' ', subset_vcf_allsamples, ' > ', subset_vcf_singlesample)
  bgzip_cmd_2 = paste0('bgzip ', subset_vcf_singlesample)
  index_cmd = paste0('tabix -p vcf ', subset_vcf_singlesample, '.gz')

  if (!file.exists(subset_vcf_allsamples)){
    system(tabix_cmd)
    system(bgzip_cmd_1)
  }

  if (!file.exists(paste0(subset_vcf_singlesample, '.gz.tbi'))){
    system(bcftools_cmd)
    system(bgzip_cmd_2)
    system(index_cmd)
  }
  whatshap_cmd = paste0(
    whatshap_bin, 
    ' haplotag ',
    '-o ',  aln_bam, '_tagged.bam ',
    '--reference ', hg38_fa, ' ',
    subset_vcf_singlesample,'.gz', ' --sample ', vcf_samplename, ' ',
    '--output-haplotag-list ', aln_bam, '_tags.tsv --ignore-read-groups ',
    aln_bam
  )
  
  print(whatshap_cmd)
  system(whatshap_cmd)
  
}

collect_whatshap_res <- function(haptags, sample, region, hap, summarylist_link){
     awk_cmd = paste0('tail -n +2 ', haptags, ' | awk \'BEGIN {FS=OFS="\t"} {print "', sample, '\t', hap, '\t', region, '",$1,$2}\' > ', summarylist_link)
     print(awk_cmd)
     system(awk_cmd)
     print('done')
  
}

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
  none_fract = xx$none / (xx$H1 + xx$H2 + xx$none)

  # Determine mapped haplotype
  if (h1_to_h2_fract > 0.5){
    mapped_hap = '1'
  } else if (h1_to_h2_fract < 0.5){
    mapped_hap = '2'
  } else {
    mapped_hap = 'unclear'
  }

  # Determine action 
  if (asm_hap == mapped_hap){
    action = 'keep'
  } else if ((asm_hap == 1) & (mapped_hap == 2) | (asm_hap == 2) & (mapped_hap == 1)){
    action = 'flip'
  } else if (mapped_hap == 'unclear'){
    action = 'unclear'
  }


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




args <- parser$parse_args()


if (args$function_name == "aln_chunks_to_minimap") {
  aln_chunks_to_minimap(args$res_path, args$region, args$sample, args$hap, args$hg38_mmi, args$mm2_bin, args$samtools_bin, args$bedtools_bin, args$chunklen)
} else if (args$function_name == "determine_phase_with_whatshap") {
  determine_phase_with_whatshap(args$aln_bam, args$region, args$sample, args$hap, args$hg38_fa, args$vcf_dir, args$subset_vcf_allsamples, args$whatshap_bin)
} else if (args$function_name == "collect_whatshap_res") {
  collect_whatshap_res(args$haptags, args$sample, args$region, args$hap, args$summarylist_link)
} else if (args$function_name == "evaluate_summarylist") {
  evaluate_summarylist(args$summarylist, args$actionlist)
}