#!/usr/bin/Rscript

##---------------------------------------------------------
# Compute he proportion of read pairs which map to more than one gene
# Input: 
# 	dir		directory contains the t2g.txt, .M and .k files
# 	files		file name of .M and .k, keep them same
# 			look at misc/transcript_pool_estimation/cmd.txt for 
# 			example to create t2g.txt, .M and .k files
# 	ncores		number of threads
# Example:
# 	Rscript paralog.R ~/rsss/human expr 8
##---------------------------------------------------------

# options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
	cat(paste("Less than three args.", "Usage: Rscript paralog.R [DIR] [FILENAME] [NCORES]", "e.g.: paralog.R ~/rsss/human expr 8", sep="\n"), "\n")
	quit()
}

dir <- args[1]
file <- args[2]
ncore <- as.integer(args[3])

require(Matrix)
require(parallel)

transgene.ids = read.table(paste(dir, "/t2g.txt", sep=""), header = F, col.names = c("transcript.id", "gene.id"))
cat(paste("Read ", dir, "/t2g.txt", sep=""), "\n")

k = read.delim(paste(dir, "/", file, ".k", sep=""), header = F, col.names="reads")
cat(paste("Read ", dir, "/", file, ".k", sep=""), "\n")

transcripts = system(paste("head -n 1 ", dir, "/", file, ".M | tr \"\\t\" \"\\n\"", sep=""), wait=T, intern = T)
transcripts = transcripts[grep("^#", transcripts, invert = TRUE)]
cat(paste("Read header ", dir, "/", file, ".M", sep=""), "\n")

m = read.delim(paste(dir, "/", file, ".M", sep=""), skip = 1, header = F, col.names = c("line.index", "transcript"))
m = m+1
cat(paste("Read table ", dir, "/", file, ".M", sep=""), "\n")

sm=sparseMatrix(m[,1],m[,2])
gene.count.test <- mclapply(seq(1, dim(sm)[1], by=1), function (i) {
  sm.idx = which(sm[i,]=="TRUE")
  t.names = transcripts[sm.idx]
  length(unique(transgene.ids$gene.id[which(transgene.ids$transcript.id %in% t.names)]))
}, mc.cores=ncore)

ul = unlist(gene.count.test)
p = sum(k[which(ul > 1),])/sum(k)
cat(paste("The proportion of read pairs which map to more than one gene: ", p, sep=""), "\n")