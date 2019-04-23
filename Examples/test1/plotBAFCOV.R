vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


#### Input Parameters ######
FileIn <- split.vars[1]
fileout <- split.vars[2]


systemwcl<-paste("wc -l ",FileIn,sep="")
TotalLines<-system(systemwcl,intern = TRUE)
TotalLines<-as.numeric(unlist(strsplit(TotalLines," "))[1])

con <- file(FileIn, "r")

TableCountTotal<-matrix("",nrow=TotalLines,ncol=13)
Nlines<-50000
control<-1
StartRow<-1
while (control==1)
{
  TableIn<-readLines(con, n = Nlines)
  NRead<-length(TableIn)
  TableInSplit<-unlist(strsplit(TableIn,"\t"))
  EndRow<-StartRow+NRead-1
  
  indChr<-seq(1,length(TableInSplit),by=13)
  indPos<-seq(2,length(TableInSplit),by=13)
  indRef<-seq(3,length(TableInSplit),by=13)
  indA<-seq(4,length(TableInSplit),by=13)
  indT<-seq(5,length(TableInSplit),by=13)
  indC<-seq(6,length(TableInSplit),by=13)
  indG<-seq(7,length(TableInSplit),by=13)
  inda<-seq(8,length(TableInSplit),by=13)
  indt<-seq(9,length(TableInSplit),by=13)
  indc<-seq(10,length(TableInSplit),by=13)
  indg<-seq(11,length(TableInSplit),by=13)
  indIn<-seq(12,length(TableInSplit),by=13)
  indDel<-seq(13,length(TableInSplit),by=13)
  
  TableCountTotal[c(StartRow:EndRow),]<-cbind(TableInSplit[indChr],TableInSplit[indPos],TableInSplit[indRef],TableInSplit[indA],TableInSplit[indT],TableInSplit[indC],TableInSplit[indG],TableInSplit[inda],TableInSplit[indt],TableInSplit[indc],TableInSplit[indg],TableInSplit[indIn],TableInSplit[indDel])
  StartRow<-EndRow+1
  if (NRead<Nlines)
  {
    control<-0
  }
}

close(con)

#### Filtering Position with small Insertions or Deletions ####
InCol<-TableCountTotal[,12]
DelCol<-TableCountTotal[,13]

indInDelF<-unique(c(which(InCol!="NA"),which(DelCol!="NA")))

if (length(indInDelF)!=0)
{
TableCountTotal<-TableCountTotal[-indInDelF,]
}


#### Filtering Strand Bias ####
RefCol<-TableCountTotal[,3]
BaseVec<-c("A","T","C","G")
l<-1
TableCountFinal<-matrix("Z",nrow=nrow(TableCountTotal),ncol=4)
for (i in 1:nrow(TableCountTotal))
{
  indRef<-which(BaseVec==RefCol[i])
  forwardCount<-as.numeric(TableCountTotal[i,c(4:7)])
  reverseCount<-as.numeric(TableCountTotal[i,c(8:11)])
  forwardRef<-forwardCount[indRef]
  forwardAlt<-sum(forwardCount[-indRef])
  reverseRef<-reverseCount[indRef]
  reverseAlt<-sum(reverseCount[-indRef])
  TableCont<-rbind(c(forwardRef,forwardAlt),c(reverseRef,reverseAlt))
  pvalue<-fisher.test(TableCont,alternative = "two.sided")$p.value
  if (pvalue>0.1)
  {
    TableCountFinal[l,]<-c(as.character(TableCountTotal[i,1]),as.character(TableCountTotal[i,2]),sum(forwardCount)+sum(reverseCount),max(forwardCount[-indRef])+max(reverseCount[-indRef]))
    l<-l+1
  }
}

TableCountFinal<-TableCountFinal[c(1:(l-1)),]
TableCountFinal<-data.frame(TableCountFinal)


colnames(TableCountFinal)<-c("chr","pos","dp","altdp")
mycount.filtered<-TableCountFinal

##### PLOT BAF AND NORMCOV ###

pos<-as.numeric(mycount.filtered[,2])
BAF<-as.numeric(mycount.filtered[,4])/as.numeric(mycount.filtered[,3])
cov<-as.numeric(mycount.filtered[,3])
med<-median(cov)
cov2<-cov/median


colors<-rep('black', nrow(mycount.filtered))
indhetdel<-which(pos >= 20000000 & pos <= 22000000)
indhetdup<-which(pos >= 35000000 & pos <= 37000000)

colors[indhetdel]<-'red'
colors[indhetdup]<-'red'


pdf(fileout,height=10,width=15)
par(mfrow=c(2,1))
plot(pos, BAF, cex=0.1, ylim=c(-1,1), xlab="coordinates",ylab="BAF", col=colors)
plot(pos, cov2, ylim=c(-1,1),xlab="coordinates",ylab="coverage (norm)", col=colors)
dev.off()

