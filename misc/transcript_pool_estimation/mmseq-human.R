#!/usr/bin/Rscript

source("mmseq-latest/src/R/mmseq.R")

thres= -6
#exp(thres)

ms.gene=readmmseq("bodymapliver.gene.mmseq")
#hist(ms.gene$log_mu)
#sum(ms.gene$log_mu>thres)/length(ms.gene$log_mu)
#sum(ms.gene$observed)/length(ms.gene$observed)
#sum(ms.gene$log_mu>thres)/length(ms.gene$observed)

ms.trans=readmmseq("bodymapliver.mmseq")
#sum(ms.trans$observed)/length(ms.trans$observed)
#sum(ms.trans$log_mu>thres)/length(ms.trans$observed)

t2g <- read.table("t2g.txt")
# expressed/active genes: ms.gene$log_mu>thres
trans.from.expr.genes=t2g[t2g[,2] %in% rownames(ms.gene$log_mu)[ms.gene$log_mu>thres],1]

# total no. of genes 44635
length(ms.gene$log_mu)
# total no. of active genes 18374
sum(ms.gene$log_mu>thres)
# total no. of transcripts 170778
length(ms.trans$log_mu)
# total no. of expressed transcripts 46223
sum(ms.trans$log_mu>thres)
# total no. of transcripts in active genes 128604
length(ms.trans$log_mu[match(trans.from.expr.genes, rownames(ms.trans$log_mu))])
# total no. of expressed transcripts in active genes 46223
sum(ms.trans$log_mu[match(trans.from.expr.genes, rownames(ms.trans$log_mu))] > thres)

# f_g the proportion of active genes. est 0.4116501
f_g.thres=sum(ms.gene$log_mu>thres)/length(ms.gene$log_mu) 
# f_t the proportion of transcripts expressed from active genes. est 0.3594212
f_t.thres=sum(ms.trans$log_mu[match(trans.from.expr.genes, rownames(ms.trans$log_mu))] > thres)/length(ms.trans$log_mu[match(trans.from.expr.genes, rownames(ms.trans$log_mu))])

message("f_g.thres: ", f_g.thres)
message("f_t.thres: ", f_t.thres)

# or estimate by observed
# f_g 0.519077
f_g.observed=sum(ms.gene$observed)/length(ms.gene$observed)

# f_t 0.8954853
trans.from.expr.genes.observed=t2g[t2g[,2] %in% rownames(ms.gene$observed)[ms.gene$observed == 1],1]
f_t.observed=sum(ms.trans$observed[match(trans.from.expr.genes.observed, rownames(ms.trans$observed))])/length(ms.trans$observed[match(trans.from.expr.genes.observed, rownames(ms.trans$observed))])

message("f_g.observed: ", f_g.observed)
message("f_t.observed: ", f_t.observed)
