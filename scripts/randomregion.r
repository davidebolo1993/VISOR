#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos='http://cran.us.r-project.org')
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",repos='http://cran.us.r-project.org')
if (!requireNamespace("regioneR", quietly = TRUE))
  BiocManager::install("regioneR")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(regioneR))

option_list = list(
  make_option(c("-d", "--dimensions"), action="store", type='character', help="TSV file with chromosome dimensions (as result from 'cut -f1,2 genome.fa.fai') [required]"),
  make_option(c("-n", "--number"), action="store", default=100, type='numeric', help="number of intervals [100]"),
  make_option(c("-l", "--length"), action="store", default=200000, type='numeric', help="mean intervals length [200000]"),
  make_option(c("-s", "--standarddev"), action="store", default=0, type='numeric', help="standard deviation for intervals length [0]"),
  make_option(c("-x", "--exclude"), action="store", default=NULL, type='character', help="exclude regions in BED [optional]"),
  make_option(c("-v", "--variants"), action="store", type='character', help="variants types (-v 'deletion,tandem duplication,inversion') [required]"),
  make_option(c("-r", "--ratio"), action="store", type='character', help="variants proportions (-r '30:30:40')) [required]"),
  make_option(c("-i", "--idhaplo"), action="store", type='numeric', default=1, help="haplotype number [1]")
)

opt = parse_args(OptionParser(option_list=option_list))

possiblevariants<-c('deletion', 'insertion', 'inversion', 'tandem duplication', 'inverted tandem duplication', 'translocation copy-paste', 'translocation cut-paste' ,  'reciprocal translocation')

if (is.null(opt$dimensions)) {
  stop('-d/--dimensions TSV file is required')
} else {
  genome<-read.table(file.path(opt$dimensions), sep='\t', header = F)
}


if (is.null(opt$exclude)) {
  exclude<-NULL
} else {
  exclude<-read.table(file.path(opt$exclude), sep='\t', header=F)
}

if (is.null(opt$variants)) {
  stop('variants in -v/--variants are required')
}

if (is.null(opt$ratio)) {
  stop('variants ratios in -r/--ratio are required')
}

variants<-unlist(strsplit(opt$variants, ','))

if (!all(variants %in% possiblevariants)) {
  stop('accepted variants are: ', paste(possiblevariants, collapse=', '))
}

number<-as.numeric(unlist(strsplit(opt$ratio, ':')))

if (length(variants) != length(number)) {
  stop('for each variant a ratio must be specified')
}

varwithproportions<-rep(variants, (opt$number/100)*number)

regions<-createRandomRegions(opt$number, length.mean=opt$length, length.sd=opt$standarddev, genome=genome, mask=exclude, non.overlapping=TRUE)

# just in case we need another list of regions to exclude for translocation
newexclude<-rbind(exclude[,1:3], cbind(as.character(seqnames(regions)), as.numeric(start(regions)-1), as.numeric(end(regions))))

df <- data.frame(chromosome=seqnames(regions),start=start(regions)-1,end=end(regions), type=varwithproportions, stringsAsFactors = FALSE)

info<-rep('INFO',nrow(df))
breakseqlen<-rep(0, nrow(df))

for(i in (1:nrow(df))) {
  if (df$type[i] == 'deletion') {
    info[i]<-'None'
    breakseqlen[i]<-sample(0:10,1)
  } else if (df$type[i] == 'tandem duplication') {
    info[i] <- '2'
    breakseqlen[i]<-sample(0:10,1)
  } else if (df$type[i] == 'inversion'){
    info[i] <- 'None'
    breakseqlen[i]<-sample(0:10,1)
  } else if (df$type[i] == 'inverted tandem duplication') {
    info[i] <- '2'
  } else if (df$type[i] == 'insertion') {
    df$start[i]<-df$end[i]-1
    num<-sample(10:1000,1)
    alphabet<-c('A', 'T', 'C', 'G')
    motif<-paste(sample(alphabet, num, replace = T), collapse='')
    info[i]<-motif
    breakseqlen[i]<-sample(0:10,1)    
  } else if (df$type[i] == 'translocation cut-paste' | df$type[i] == 'translocation copy-paste') {
    newregion<-createRandomRegions(1, length.mean=opt$length, length.sd=opt$standarddev, genome=genome, mask=newexclude, non.overlapping=TRUE)
    chromosome<-as.character(seqnames(newregion))
    start<-as.numeric(start(newregion))-1
    end<-as.numeric(end(newregion))
    orientation<-sample(c('forward', 'reverse'),1)
    newexclude<-rbind(newexclude,c(chromosome, start, end))
    info[i]<-paste(paste0('h',opt$idhaplo),chromosome,start,orientation,sep=':')
    breakseqlen[i]<-sample(0:10,1) 
  } else {
    newregion<-createRandomRegions(1, length.mean=opt$length, length.sd=opt$standarddev, genome=genome, mask=newexclude, non.overlapping=TRUE)
    chromosome<-as.character(seqnames(newregion))
    start<-as.numeric(start(newregion))-1
    end<-as.numeric(end(newregion))
    orientation1<-sample(c('forward', 'reverse'),1)
    orientation2<-sample(c('forward', 'reverse'),1)
    newexclude<-rbind(newexclude,c(chromosome, start, end))
    info[i]<-paste(paste0('h',opt$idhaplo+1),chromosome,start,orientation1, orientation2, sep=':')
    breakseqlen[i]<-sample(0:10,1) 
  }
} 

final<-data.frame(df,info,breakseqlen, stringsAsFactors = FALSE)

write.table(final,file = '',quote = FALSE, col.names = FALSE, row.names = FALSE, sep='\t')