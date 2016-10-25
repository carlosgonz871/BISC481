######################################
# 
# Retrieval of ChIP-seq data
# Generation of bound and unbound data
# BISC 481
######################################

## Install packages
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("AnnotationHub")
biocLite("GenomicRanges")
biocLite("BSgenome.Mmusculus.UCSC.mm10")
##DnaShaper already intall into R studio
#had to install package 'shiny' manually for AnnotationHub to function using code:
install.packages("shiny", repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
##caret has been previously installed

## Initialization
library(DNAshapeR)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(caret)


## Data retreive
seqLength <- 50 #500
sampleSize <- 2000 #42045
workingPath <- "C:\\Users\\carlo\\Downloads\\BISC481-master\\BISC481-master\\CTCF\\"

# Bound (ChIP-seq)
ah <- AnnotationHub()
#unique(ah$dataprovider)
#query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "bound.fa"))

# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()

while ( length(randomGr) < sampleSize )  {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:tmpChrLength,1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}
randomGr

# Overlap checking
findOverlaps(ctcfPeaks, randomGr)

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, 
         filename = paste0(workingPath, "bound.fa"))

head(fileName)




# Bound (ChIP-seq)
ah <- AnnotationHub()
#unique(ah$dataprovider)
#query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "unbound.fa"))

# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()

while ( length(randomGr) < sampleSize )  {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:tmpChrLength,1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}
randomGr

# Overlap checking
findOverlaps(ctcfPeaks, randomGr)

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, 
         filename = paste0(workingPath, "unbound.fa"))











# Bound (ChIP-seq)
ah <- AnnotationHub()
#unique(ah$dataprovider)
#query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "ctcf.fa.HelT"))

# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()

while ( length(randomGr) < sampleSize )  {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:tmpChrLength,1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}
randomGr

# Overlap checking
findOverlaps(ctcfPeaks, randomGr)

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, 
         filename = paste0(workingPath, "ctcf.fa.HelT"))

head(fileName(getFasta))









# Bound (ChIP-seq)
ah <- AnnotationHub()
#unique(ah$dataprovider)
#query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta(( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "ctcf.fa.ProT"))
# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22],
chrLength <- seqlengths(Mmusculus)[1:22],
randomGr <- GRanges(),

while ( length(randomGr) < sampleSize )  {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:tmpChrLength,1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}
randomGr,

# Overlap checking
findOverlaps(ctcfPeaks, randomGr),

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, 
         filename = paste0(workingPath, "ctcf.fa.ProT"))


