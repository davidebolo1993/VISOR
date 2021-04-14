if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos='http://cran.us.r-project.org')
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",repos='http://cran.us.r-project.org')
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table",repos='http://cran.us.r-project.org') 
if (!requireNamespace("regioneR", quietly = TRUE))
  BiocManager::install("regioneR")
if (!requireNamespace("Rsamtools", quietly = TRUE))
  BiocManager::install("Rsamtools")


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(regioneR))
suppressPackageStartupMessages(library(Rsamtools))

option_list = list(
  make_option(c("-d", "--dimensions"), action="store", type='character', help="TSV file with chromosome dimensions (as result from 'cut -f1,2 genome.fa.fai') [required]"),
  make_option(c("-n", "--number"), action="store", default=100, type='numeric', help="number of intervals [100]"),
  make_option(c("-l", "--length"), action="store", default=50000, type='numeric', help="mean intervals length [50000]"),
  make_option(c("-s", "--standarddev"), action="store", default=0, type='numeric', help="standard deviation for intervals length [0]"),
  make_option(c("-x", "--exclude"), action="store", default=NULL, type='character', help="exclude regions in BED [optional]"),
  make_option(c("-v", "--variants"), action="store", type='character', help="variants types (-v 'deletion,tandem duplication,inversion') [required]"),
  make_option(c("-r", "--ratio"), action="store", type='character', help="variants proportions (-r '30:30:40')) [required]"),
  make_option(c("-i", "--idhaplo"), action="store", type='numeric', default=1, help="haplotype number [1]"),
  make_option(c("-m", "--microsatellites"), action="store", type='character', default=NULL, help="BED file with microsatellites (ucsc format) [optional]"),
  make_option(c("-g", "--genome"), action="store", type="character", default=NULL, help="FASTA genome [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

possiblevariants<-c('deletion', 'insertion', 'inversion', 'tandem duplication', 'inverted tandem duplication', 'translocation copy-paste', 'translocation cut-paste' ,  'reciprocal translocation', 'SNP', 'MNP', 'tandem repeat contraction', 'tandem repeat expansion', 'perfect tandem repetition', 'approximate tandem repetition')
possiblemicro<-c('AAC','AAG','AAT','AC', 'ACA','ACC','ACT','AG','AGA','AGC','AGG','AT','ATA','ATC','ATG','ATT','CA','CAA','CAC','CAG','CAT','CCA','CCG','CCT','CT','CTA','CTC','CTG','CTT','GA','GAA','GAT','GCA','GCC','GGA','GGC','GT','GTG','GTT','TA','TAA','TAC','TAG','TAT','TC','TCA','TCC','TCT','TG','TGA','TGC','TGG','TGT','TTA','TTC','TTG')

if (is.null(opt$dimensions)) {
  stop('-d/--dimensions TSV file is required')
} else {
  genome<-fread(file.path(opt$dimensions), sep='\t', header = F)
}

if (is.null(opt$exclude)) {
  exclude<-NULL
} else {
  exclude<-fread(file.path(opt$exclude), sep='\t', header=F)
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

testt<-c("tandem repeat contraction", "tandem repeat expansion")

if (any(testt %in% variants)) {
  if (is.null(opt$microsatellites)) {
    stop('when adding contractions/expansions of known microsatellites, a BED file specifying known microsatellites must be given')
  } else {
    microbed<-fread(file.path(opt$microsatellites), sep='\t', header=FALSE)
    if (!is.null(exclude)) {
      exclude<-data.frame(rbind(exclude[,1:3],data.frame(V1=microbed$V1, V2=microbed$V2-1, V3=microbed$V3+1)))
    } else {
      exclude<-data.frame(V1=microbed$V1, V2=microbed$V2-2, V3=microbed$V3+2)
    }
  }
}

testv<-c("SNP", "MNP")

if (any(testv %in% variants)) {
  if (is.null(opt$genome)) {
    stop('when adding SNPs and MNPs, a genome FASTA must be given')
  } else {
    if (!file.exists(file.path(paste0(opt$genome, '.fai')))) {
      indexFa(file.path(opt$genome))
    }
  }
}

number<-as.numeric(unlist(strsplit(opt$ratio, ':')))

if (cumsum(number)[length(cumsum(number))] != 100) {
  stop('sum of variant ratios must be 100')
}

if (length(variants) != length(number)) {
  stop('for each variant a ratio must be specified')
}

varwithproportions<-rep(variants, (opt$number/100)*number)

#the number of regions here is all minus the number of repeat contractions/expansions for which we sample from known regions
keepdist<-varwithproportions[which(varwithproportions %in% testt)]
dft<-data.frame(matrix(NA, nrow = length(keepdist), ncol = 6), stringsAsFactors=FALSE)
colnames(dft)<-c("chromosome","start","end","type","info","breakseqlen")
if (length(keepdist) != 0) {
  varwithproportions<-varwithproportions[-which(varwithproportions %in% testt)] #remove
  randomicro<-microbed[sample(nrow(microbed), length(keepdist)),]
  for (i in 1:nrow(randomicro)) {
    dft$chromosome[i]<-randomicro$V1[i]
    dft$start[i]<-as.numeric(randomicro$V2[i])
    dft$end[i]<-as.numeric(randomicro$V3[i])
    dft$type[i]<-keepdist[i]
    motif<-unlist(strsplit(randomicro$V4[i],'x'))[2]
    number_<-as.numeric(unlist(strsplit(randomicro$V4[i],'x'))[1])
    minimum<-1
    if (keepdist[i] == 'tandem repeat contraction') {
      maximum<-number_-1
    } else {#is expansion
      maximum<-500
    }
    dft$info[i]<-paste(motif, sample(minimum:maximum,1),sep=':')
    dft$breakseqlen[i]<-0
  }
}

regions<-createRandomRegions(length(varwithproportions), length.mean=opt$length, length.sd=opt$standarddev, genome=genome, mask=exclude, non.overlapping=TRUE)

# just in case we need another list of regions to exclude for translocation
newexclude<-data.frame(rbind(exclude[,1:3], cbind(as.character(seqnames(regions)), as.numeric(start(regions)-2), as.numeric(end(regions)+2))))

df <- data.frame(chromosome=seqnames(regions),start=start(regions),end=end(regions), type=varwithproportions, stringsAsFactors = FALSE)

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
    num<-round(rnorm(1, mean=opt$length, sd=opt$standarddev)) #adapt mean and stdev for insertions to those specified by the user
    alphabet<-c('A', 'T', 'C', 'G')
    motif<-paste(sample(alphabet, num, replace = T), collapse='')
    info[i]<-motif
    breakseqlen[i]<-sample(0:10,1)    
  } else if (df$type[i] == 'perfect tandem repetition') {
    df$start[i]<-df$end[i]-1
    motif<-sample(possiblemicro,1)
    num<-sample(6:500,1)
    info[i]<-paste(motif,num,sep=':')
    breakseqlen[i]<-0
  } else if (df$type[i] == 'approximate tandem repetition') {
    df$start[i]<-df$end[i]-1
    motif<-sample(possiblemicro,1)
    num<-sample(6:500,1)
    err<-sample(1:round(num/3),1)
    info[i]<-paste(motif,num,err,sep=':')
    breakseqlen[i]<-0
  } else if (df$type[i] == 'SNP') {
    df$start[i]<-df$end[i]-1
    tog<-makeGRangesFromDataFrame(data.frame(chromosome=df$chromosome[i], start=df$end[i], end=df$end[i]))
    nuc<-as.character(scanFa(opt$genome,tog))[[1]]
    if (nuc == 'N') {
      alphabet<-c('A', 'T', 'C', 'G', 'N')
    } else {
      alphabet<-c('A', 'T', 'C', 'G')
    }
    alphabet<-alphabet[-which(alphabet == nuc)]
    info[i]<-sample(alphabet,1)
    breakseqlen[i]<-0
  } else if(df$type[i] == 'MNP') {
    smallseq<-sample(2:10,1)
    df$start[i]<-df$end[i]-smallseq
    tog<-makeGRangesFromDataFrame(data.frame(chromosome=df$chromosome[i], start=df$start[i], end=df$end[i]))
    nucseq<-as.character(scanFa(opt$genome,tog))[[1]]
    for (l in 1:nchar(nucseq)) {
      nuc<-substr(nucseq, l, l)
      if (nuc == 'N') {
        alphabet<-c('A', 'T', 'C', 'G', 'N')
      } else {
        alphabet<-c('A', 'T', 'C', 'G')
      }
      alphabet<-alphabet[-which(alphabet == nuc)]
      substr(nucseq, l, l)<-sample(alphabet,1)
    }
    info[i]<-nucseq
    breakseqlen[i]<-0
  } else if (df$type[i] == 'translocation cut-paste' | df$type[i] == 'translocation copy-paste') {
    newregion<-createRandomRegions(1, length.mean=opt$length, length.sd=opt$standarddev, genome=genome, mask=newexclude, non.overlapping=TRUE)
    chromosome<-as.character(seqnames(newregion))
    start<-as.numeric(start(newregion))
    end<-as.numeric(end(newregion))
    orientation<-sample(c('forward', 'reverse'),1)
    newexclude<-data.frame(rbind(newexclude,c(chromosome, start-2, end+2)))
    info[i]<-paste(paste0('h',opt$idhaplo),chromosome,start,orientation,sep=':')
    breakseqlen[i]<-sample(0:10,1) 
  } else {
    newregion<-createRandomRegions(1, length.mean=opt$length, length.sd=opt$standarddev, genome=genome, mask=newexclude, non.overlapping=TRUE)
    chromosome<-as.character(seqnames(newregion))
    start<-as.numeric(start(newregion))
    end<-as.numeric(end(newregion))
    orientation1<-sample(c('forward', 'reverse'),1)
    orientation2<-sample(c('forward', 'reverse'),1)
    newexclude<-data.frame(rbind(newexclude,c(chromosome, start-2, end+2)))
    info[i]<-paste(paste0('h',opt$idhaplo+1),chromosome,start,orientation1, orientation2, sep=':')
    breakseqlen[i]<-sample(0:10,1) 
  }
} 

final<-data.frame(df,info,breakseqlen, stringsAsFactors = FALSE)

fwrite(rbind(final,dft),file = '',quote = FALSE, col.names = FALSE, row.names = FALSE, sep='\t')

