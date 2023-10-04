# We try out the SNPRelate package here. 
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(grid)
library(vcfR)
library(dendextend)

extract_genotypes_for_variant <- function(vcf, variant_id) {
  # Extract the genotypes
  gt <- extract.gt(vcf)
  
  # Find the index of the desired variant
  
  # The old line was: variant_index <- which(row.names(gt) == variant_id)
  # Here instead find everything that has the part in the variant id
  
  variant_index <- which(grepl(variant_id, row.names(gt)))

  # If the variant is not found, return NULL
  if (length(variant_index) == 0) return(NULL)
  
  # Turn NA into './.'
  to_return = gt[variant_index,]
  to_return[is.na(to_return)] = './.'
  
  # if to_return is a character and not a data frame (because only one entry), convert it to a 1-dim dataframe anyway
  if (is.null(nrow(to_return))){
    to_return = t(as.data.frame(to_return))
    row.names(to_return) = variant_id
  }
  # Otherwise, return the genotypes for the variant
  return(to_return)
}

convert_genotypes_to_numeric <- function(genotypes) {
  sapply(genotypes, function(gt) {
    if (gt == "0|0" || gt == "0/0") return(0)
    if (gt == "0|1" || gt == "1|0" || gt == "0/1" || gt == "1/0") return(1)
    if (gt == "1|1" || gt == "1/1") return(2)
    if (grepl("1", gt) && grepl("\\.", gt)) return(1.5)
    if (grepl("0", gt) && grepl("\\.", gt)) return(0.5)
    if (gt == "./.") return(3)  # If genotype is missing, return 3
    return(NA)  # If none of the above match, return NA
  })
}

# Modified to accept multiple IDs and return a data frame
extract_convert_multi <- function(vcf_link, ids){
  vcf = read.vcfR(vcf_link)
  
  # Initialize an empty list to store results
  results_list <- list()
  gtstrings = extract_genotypes_for_variant(vcf, ids)
  #results_list[[id]] = numeric_genotypes
  for(gtrow in 1:nrow(gtstrings)) {
    numeric_genotypes = convert_genotypes_to_numeric(gtstrings[gtrow,])
    results_list[[gtrow]] = numeric_genotypes
  }
  # Convert list to data frame
  df <- do.call(data.frame, results_list)
  names(df) <- row.names(gtstrings)
  
  return(df)
}

# Map NW_gts values to palette indices. This assumes NW_gts can only have values: 0, 0.5, 1, 1.5, 2, and 3
value_to_index <- function(value) {
  mapping <- c("0" = 1, "0.5" = 2, "1" = 3, "1.5" = 4, "2" = 5, "3" = 6)
  return(mapping[as.character(value)])
}

plot_scatter <- function(x, pc.percent, superpopCol) {
  # Plotting
  txtFontSize=16
  axisFontSize=16
  axisTtlFontSize=22
  lgdTtlFontSize=22
  lgdFontSize=16

  # Fist, we plot as a scatterplot.
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,1)))
  xttl=paste("PC1 (",pc.percent[1],"%)",sep="")
  yttl=paste("PC2 (",pc.percent[2],"%)",sep="")
  p1=ggplot(x, aes(x=PC1, y=PC2, fill=ancestry)) + geom_point(shape=21, colour="black", size=2, alpha=0.5) + xlab(xttl) + ylab(yttl)
  p1=p1 + scale_fill_manual(values=as.character(superpopCol)) 
  xttl=paste("PC3 (",pc.percent[3],"%)",sep="")
  yttl=paste("PC4 (",pc.percent[4],"%)",sep="")
  p2=ggplot(x, aes(x=PC3, y=PC4, fill=ancestry)) + geom_point(shape=21, colour="black", size=2, alpha=0.5) + xlab(xttl) + ylab(yttl)
  p2=p2 + scale_fill_manual(values=as.character(superpopCol)) + labs(colour="Ancestry")
  print(p1, vp = viewport(layout.pos.row=1, layout.pos.col=1))
  print(p2, vp = viewport(layout.pos.row=2, layout.pos.col=1))
}

load_and_prepare <- function(vcf_gz, gds_file, ld.threshold){


  # Load genotypes
  snpgdsVCF2GDS(vcf_gz, gds_file)
  genofile <- snpgdsOpen(gds_file, readonly=T, allow.duplicate=T)

  return(genofile)
}
# Do PCA
do_pca <- function(genofile, pop_file, ld.threshold){

  # load the pop_file
  pop <- read.table(pop_file, sep="\t", header=TRUE)

  set_seed <- 1000

  snpset <- snpgdsLDpruning(genofile, ld.threshold=ld.threshold)#, maf=0.01, missing.rate=0.8, slide.max.bp = 500000)
  snpset.id <- unlist(unname(snpset))
  pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

  pc.percent <- pca$varprop*100

  sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  tab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1], 
                    PC2 = pca$eigenvect[,2], 
                    PC3 = pca$eigenvect[,3], 
                    PC4 = pca$eigenvect[,4], 
                    stringsAsFactors = FALSE)
  # Add populations
  tab = left_join(tab, pop, by='sample.id')
  tab$ancestry = as.factor(tab$pop)
  tab$pop = NULL
  x = tab
  x$ancestry=factor(x$ancestry, levels=c("AFR", "AMR", "EAS", "EUR", "SAS"))

  return(list(x, pc.percent))
}


extract_genotypes_by_variant <- function(vcf_gz, snp_id_targets, samp.order){
  NW_gts <- extract_convert_multi(vcf_gz, snp_id_targets)
  NW_gts = NW_gts[samp.order,,drop=F]
  return(NW_gts)
}

prepare_dendrogram_data <- function(genofile, x, superpopCol, snp_id_targets, my_palette) {
  
  # QC: Check if inputs are non-null and appropriately defined
  stopifnot(!is.null(genofile), !is.null(x), !is.null(superpopCol), !is.null(snp_id_targets))
  
  # Compute the IBS matrix
  ibs <- snpgdsIBS(genofile, num.thread=2)
  
  # Hierarchical clustering on the IBS matrix
  ibs.hc <- snpgdsHCluster(ibs)
  rv <- snpgdsCutTree(ibs.hc)
  
  # Assign colors
  dend <- rv$dendrogram
  labels_colors(dend) <- superpopCol[x$ancestry[match(labels(dend), x$sample.id)]]
  NW_gts <- extract_genotypes_by_variant(vcf_gz, snp_id_targets, rv$samp.order)
  # Color coding for the colored_bars 
  indices_matrix <- apply(NW_gts, 2, function(col) sapply(col, value_to_index))
  NW_gts_colors_matrix <- matrix(my_palette[as.vector(indices_matrix)], ncol = ncol(NW_gts))
  
  colnames(NW_gts_colors_matrix) = colnames(NW_gts)
  row.names(NW_gts_colors_matrix) = row.names(NW_gts)
  
  return(list(dend = dend, NW_gts_colors_matrix = NW_gts_colors_matrix))
}

plot_dendrogram <- function(dend, x, superpopCol, NW_gts_colors_matrix, my_palette, bottom, left, top, right, title) {
  
  # QC: Check if inputs are non-null and appropriately defined
  stopifnot(!is.null(dend), !is.null(x), !is.null(superpopCol), !is.null(NW_gts_colors_matrix))
  
  legend_colors <- my_palette[c(1,2,3,4,5,6)]
  legend_labels <- c("0.0", "0/.", "0/1", "1/.", "1/1", "./.")
  
  par(mar=c(left, bottom, top, right))
  plot(hang.dendrogram(dend, hang = -1), ylab = "Height", leaflab="perpendicular", main=title)
  colored_bars(dend = hang.dendrogram(dend, hang = -1), NW_gts_colors_matrix, rowLabels = colnames(NW_gts_colors_matrix), sort_by_labels_order = F)
  
  legend('topright', legend=levels(x$ancestry), fill=superpopCol, col=superpopCol, xpd=NA)
  legend(x=5, y=-0.08, legend=legend_labels, fill=legend_colors, title="Values", horiz=TRUE, bty="n", xpd=NA)
}



main <- function(vcf_gz, pop_file, gds_file, snp_id_targets, ld.threshold, title){

  gds_file=paste0('trash', as.numeric(runif(1)),'90230wddd3.gds')
  ld.threshold = 0.5
  # Load stuff
  genofile = load_and_prepare(vcf_gz, gds_file, ld.threshold)
  
  # Do PCA
  x_pcapct_list <- do_pca(genofile, pop_file, ld.threshold)
  x <- x_pcapct_list[[1]]
  pc.percent <- x_pcapct_list[[2]]

  superpopCol <- c(AFR="#FFCD33", AMR="#FF3D3D", EAS="#ADFF33", EUR="#64EBFF", SAS="#FF30FF")

  # The colors defined here are for the dendrogram:
  # 0/0, 0/., 0/1, 1/., 1/1, ./.
  my_palette <- colorRampPalette(c("grey", "white", "yellow", "yellow", "red", "white"))(6)
  
  # Plottings
  # plot_scatter(x, pc.percent, superpopCol)
  
  # Now towards the dendro
  dendro_data <- prepare_dendrogram_data(genofile, x, superpopCol, snp_id_targets, my_palette)
  plot_dendrogram(dendro_data$dend, x, superpopCol, dendro_data$NW_gts_colors_matrix, my_palette, 18,15,1,1, title)
}





#vcf_gz = "/Users/hoeps/PhD/projects/nahrcall/revisions/popgen/data/testcases/vcf_of_that/final_chr1.vcf.gz"
#vcf_gz = "/Users/hoeps/PhD/projects/nahrcall/revisions/popgen/data/testcases/vcf_of_that/snptest.vcf.gz"
vcf_gz = "/Users/hoeps/PhD/projects/nahrcall/revisions/popgen/data/testcases/vcf_of_that/out_fixfix.vcf.gz"
pop_file = "/Users/hoeps/PhD/projects/nahrcall/revisions/popgen/data/ancestries/anc_only.tsv"

snp_id_targets = c('<NWhal.inv+del_chr1-119747586-121609789>', 
                   #'<NWhal.inv+dup+del_chr1-119747586-121609789>',
                   #'<NWhal.inv+inv+inv_chr1-119747586-121609789>',
                   '<NWhal.inv_chr1-119747586-121609789>')

snp_id_targets = c(
  '<NWhal.del+inv_chr11-50055000-50424100>',
  '<NWhal.inv+del_chr11-50055000-50424100>',
  '<NWhal.inv+dup_chr11-50055000-50424100>',
  '<NWhal.inv+inv+del_chr11-50055000-50424100>',
  '<NWhal.inv_chr11-50055000-50424100>'     
)

snp_id_targets = c(
  '<NWhal.inv_chr2-138240963-138257877>',
  '<NWhal.inv+dup_chr11-50055000-50424100>'
)
main(vcf_gz, pop_file, 'whatever.gds', 'hr17-18499564-18948', 0.5, '')






