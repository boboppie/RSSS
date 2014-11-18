#!/usr/bin/env python

import sys
sys.path.append("core")
sys.path.append("sens_prec")
import ductape as tp
from steps import *
from sens_prec_pipeline import *
import rpy2
import rpy2.robjects
import itertools
import collections

if __name__ == '__main__':
    #unittest.main(exit=True, verbosity=2)

    # Compute annotations + RABT for varying (sens, prec)
    coverage = 1000
    #PBS -t 0-24
    lSensitivity = [0.2, 0.4, 0.6, 0.8, 1.0]
    lPrecision = [0.2, 0.25, 0.3, 0.35, 0.4]

    # Limit to a single (sens, prec) if a cluster thread
    jobid = jobid()
    if not (jobid is None):
        sensPrec = [(sens, prec) for (sens, prec) in itertools.product(lSensitivity, lPrecision)][jobid]
        lSensitivity = [ sensPrec[0] ]
        lPrecision = [ sensPrec[1] ]

    # Compute/load (sens, prec)-independant computations
    transcriptPool = TranscriptPool2(flSubset=True, seedVal=1)
    simulatedReads = SimulatedReads2(nTotalFrags=len(transcriptPool.lTruthSet)*coverage, transcriptPool=transcriptPool)#, recompute=True)
    tophatBowtieIndex = TophatBowtieIndex(transcriptPool=transcriptPool)#,recompute=True)
    tophatAlign = TophatAlign(simulatedReads=simulatedReads, tophatBowtieIndex=tophatBowtieIndex)#,recompute=True)

    # Compute/load (sens, prec)-dependant computations
    dRes = {}
    for sensitivity in lSensitivity:
        for precision in lPrecision:
            dRes[(sensitivity, precision)] = pipelineCufflinks(transcriptPool=transcriptPool, simulatedReads=simulatedReads, 
                                                      tophatBowtieIndex=tophatBowtieIndex, tophatAlign=tophatAlign,
                                                      coverage=coverage, sensitivity=sensitivity, precision=precision, RABTassembly=True)

    dRes = collections.OrderedDict(sorted(dRes.items(), key=lambda t: t[0]))

    # Compute/load Cufflinks non-RABT pipeline
    cuffRes = pipelineCufflinks(transcriptPool=transcriptPool, simulatedReads=simulatedReads, 
                       tophatBowtieIndex=tophatBowtieIndex, tophatAlign=tophatAlign,
                       coverage=coverage, sensitivity=lSensitivity[0], precision=lPrecision[0], RABTassembly=False)
    # Compute/load Oases pipeline
    oasesRes = pipelineOases(transcriptPool=transcriptPool, simulatedReads=simulatedReads)

    # Exit if running a single run on the cluster
    if not (jobid is None): exit()

    # ip.embed()
    # Compute plot vectors
    plotAnSensitivity = [sensitivity for (sensitivity, precision) in dRes.keys()]
    plotAnPrecision = [precision for (sensitivity, precision) in dRes.keys()]
    plotRbSensitivity = [r.clTranscriptome.getSensitivity() for r in dRes.values()]
    plotRbPrecision = [r.clTranscriptome.getPrecision() for r in dRes.values()]
    plotAnExpressionTP = ["%.2f" % (r.plotMmseqPrediction.logCorrelationTP,) for r in dRes.values()]
    plotRbExpressionTP = ["%.2f" % (r.plotClMmseqPrediction.logCorrelationTP,) for r in dRes.values()]
    plotClExpressionTP =  "%.2f" % (cuffRes.plotClMmseqPrediction.logCorrelationTP,)
    plotOaExpressionTP =  "%.2f" % (oasesRes.oaPlotMmseqPrediction.logCorrelationTP,) 
    plotAnExpressionFP = ["%.2f" % (r.plotMmseqPrediction.meanFPExprFrac,) for r in dRes.values()]
    plotRbExpressionFP = ["%.2f" % (r.plotClMmseqPrediction.meanFPExprFrac,) for r in dRes.values()]
    plotClExpressionFP =  "%.2f" % (cuffRes.plotClMmseqPrediction.meanFPExprFrac,) 
    plotOaExpressionFP =  "%.2f" % (oasesRes.oaPlotMmseqPrediction.meanFPExprFrac,)

    R = rpy2.robjects.r
    R.source(os.path.join(os.environ["RSSS_MMSEQ_DIR"], "mmseq.R"))

    # Plot histograms of true expression and Cufflinks reconstructed expression
    # Read mmseq expression estimates
    R.assign("fn.expr.mmseq", cuffRes.clMmseqExpr.fpExprMmseq)
    R("mmseq.out <- readmmseq(mmseq_files=fn.expr.mmseq, normalize=FALSE)")
    R("expr.mmseq <- exp(mmseq.out$log_mu)")
    R("expr.mmseq[ is.na(expr.mmseq) ] <- 0.0")
    # Calculate matching vector of true expression values
    R.assign("expr.truth", cuffRes.clTranscriptome.getExpression([ s for s in R("row.names(mmseq.out$log_mu)") ]))
    R("expr.truth <- unlist(expr.truth)")
    R.assign("expr.groundTruth", transcriptPool.getExpressionAsList(transcriptPool.lTruthSet))
    R("expr.groundTruth <- unlist(expr.groundTruth)")
    R.assign("mTranscriptTP", cuffRes.clTranscriptome.getTranscriptTPness([ s for s in R("row.names(mmseq.out$log_mu)") ]))
    R("mTranscriptTP <- unlist(mTranscriptTP)")

    R("""
        # hist_clReconstructed <- hist(log(expr.mmseq[mTranscriptTP]), breaks=50)
        hist_truth <- hist(log(expr.truth[mTranscriptTP]), breaks=30)
        hist_gtruth <- hist(log(expr.groundTruth), breaks=30)
        # print(sum(hist_truth$counts)/sum(hist_gtruth$counts))                
        pdf("_logExpHis.pdf")
        plot(hist_gtruth, main="", xlab="Expression level on log scale", col="royalblue", xlim=c(2,10), ylim=c(0,130))
        plot(hist_truth, col="red", xlim=c(2,10), ylim=c(0,130), add=T)
        legend("topright", c("All transcripts", "Correctly reconstructed by Cufflinks"), pch=c(22,22), col=c("royalblue", "red"))
        dev.off()
    """)

    # Hist plots of correctly reconstructed transcripts (red) and incorrectly reconstructed transcripts (blue)
    R("""
        hist_clCorrectlyReconstructed <- hist(log(expr.mmseq[mTranscriptTP]), breaks=30)
        hist_clIncorrectlyReconstructed <- hist(log(expr.mmseq[!mTranscriptTP]), breaks=30)
        pdf("_logExpHis2.pdf", width=6, height=10)
        par(mfrow=c(2,1))
        library(scales)
        plot(hist_clCorrectlyReconstructed, main="Histograms of reconstructed transcripts", xlab="", col=alpha("royalblue", 0.5), xlim=c(2,10), ylim=c(0,60))
        legend("topright", "Correctly reconstructed by Cufflinks", pch=c(22,22), col=alpha("royalblue", 0.5))
        plot(hist_clIncorrectlyReconstructed, main="", xlab="Expression level on log scale", col=alpha("red", 0.5), xlim=c(2,10), ylim=c(0,60))
        legend("topright", "Incorrectly reconstructed by Cufflinks", pch=c(22,22), col=alpha("red", 0.5))
        dev.off()
    """)

    # Plot annotated FP expression with RABT-reconstructed FP expression
    r = dRes[(1.0, 0.4)]

    R.assign("rb.fn.expr.mmseq", r.clMmseqExpr.fpExprMmseq)
    R("rb.mmseq.out <- readmmseq(mmseq_files=rb.fn.expr.mmseq, normalize=FALSE)")
    R("rb.expr.mmseq <- exp(rb.mmseq.out$log_mu)")
    R("rb.expr.mmseq[ is.na(rb.expr.mmseq) ] <- 0.0")
    # Calculate matching vector of true expression values
    R.assign("rb.expr.truth", r.clTranscriptome.getExpression([ s for s in R("row.names(rb.mmseq.out$log_mu)") ]))
    R("rb.expr.truth <- unlist(rb.expr.truth)")
    R.assign("rb.mTranscriptTP", r.clTranscriptome.getTranscriptTPness([ s for s in R("row.names(rb.mmseq.out$log_mu)") ]))
    R("rb.mTranscriptTP <- unlist(rb.mTranscriptTP)")

    R.assign("an.fn.expr.mmseq", r.mmseqExpr.fpExprMmseq)
    R("an.mmseq.out <- readmmseq(mmseq_files=an.fn.expr.mmseq, normalize=FALSE)")
    R("an.expr.mmseq <- exp(an.mmseq.out$log_mu)")
    R("an.expr.mmseq[ is.na(an.expr.mmseq) ] <- 0.0")
    # Calculate matching vector of true expression values
    R.assign("an.expr.truth", r.clTranscriptome.getExpression([ s for s in R("row.names(an.mmseq.out$log_mu)") ]))
    R("an.expr.truth <- unlist(an.expr.truth)")
    R.assign("an.mTranscriptTP", r.clTranscriptome.getTranscriptTPness([ s for s in R("row.names(an.mmseq.out$log_mu)") ]))
    R("an.mTranscriptTP <- unlist(an.mTranscriptTP)")

    R("""
           anHis <- hist(log(an.expr.mmseq[!an.mTranscriptTP]), breaks=100)               
           rbHis <- hist(log(rb.expr.mmseq[!rb.mTranscriptTP]), breaks=100)                
        
           pdf("_AnRbFPexpComp.pdf")
           plot(anHis, main="", xlab="FP expression level on log scale (s=0.8, p=0.4)", col=rgb(0,0,1,1/4), xlim=c(0,9), ylim=c(0,30))
           plot(rbHis, col=rgb(1,0,0,1/4), add=T, xlim=c(0,9), ylim=c(0,30))
           legend("topright", c("Annotations", "RABT"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lwd=4)

           dev.off()
    """)

    # Q-Q plot and density plot of est. annoFPs aginst est. newFPs of RABT under different s and p
    R("""
        pdf("_rabtFPQQ.pdf")
        #pdf("_logRabtFPQQ.pdf")
        #pdf("_rabtFPDen.pdf")
        #pdf("_logRabtFPDen.pdf")
        par(mfrow=c(5,5), mar=c(1,1,2,1),oma=c(4,1,2,0),mgp=c(3, 0.5, 0))
        #par(mfrow=c(5,5), mar=c(5.1,4.1,4.1,2.1), oma=c(2,2,2,2), mgp=c(3, 0.5, 0))
    """)
    for (s, p) in dRes.keys():
            r = dRes[(s, p)]
            R.assign("sens", s)
            R.assign("prec", p)
            R.assign("rb.fn.expr.mmseq", r.clMmseqExpr.fpExprMmseq)
            R("rb.mmseq.out <- readmmseq(mmseq_files=rb.fn.expr.mmseq, normalize=FALSE)")
            R("rb.expr.mmseq <- exp(rb.mmseq.out$log_mu)")
            R("rb.expr.mmseq[ is.na(rb.expr.mmseq) ] <- 0.0")

            annoFPs = set(r.simulatedTranscriptome.lFP)
            totals = set(R("row.names(rb.mmseq.out$log_mu)")) # total transcripts
            matches = set(r.clTranscriptome.knownSetTranscriptsMatch) # total TPs
            totalTPs = totals.intersection(matches)
            totalFPs = totals.difference(totalTPs)
            newFPs = totalFPs.difference(annoFPs)

            annoFPness = [ s in annoFPs for s in R("row.names(rb.mmseq.out$log_mu)") ]
            newFPness = [ s in newFPs for s in R("row.names(rb.mmseq.out$log_mu)") ]
            R.assign("annoFPs", annoFPness)
            R.assign("newFPs", newFPness)
            R("annoFPs <- unlist(annoFPs)")
            R("newFPs <- unlist(newFPs)")
            R("""
                qqplot(rb.expr.mmseq[annoFPs], rb.expr.mmseq[newFPs], main=paste("s=", sens, ", p=", prec), xlab="annoFPs", ylab="newFPs")
                #qqplot(log(rb.expr.mmseq[annoFPs]), log(rb.expr.mmseq[newFPs]), main=paste("s=", sens, ", p=", prec), xlab="annoFPs", ylab="newFPs")
                #plot(density(rb.expr.mmseq[annoFPs]),col="blue",lwd=2, main=paste("s=", sens, ", p=", prec), xlab="RABT FPs", ylab="")
                #lines(density(rb.expr.mmseq[newFPs]),lwd=2,col="red")
                #plot(density(log(rb.expr.mmseq[annoFPs])),col="blue",lwd=2, main=paste("s=", sens, ", p=", prec), xlab="RABT FPs", ylab="")
                #lines(density(log(rb.expr.mmseq[newFPs])),lwd=2,col="red")

            """)
    R("dev.off()")

    # one examplar s=0.6, p=0.4
    R("""
        pdf("_logRabtFPDenExamplar.pdf")
    """)
    s = 0.6
    p = 0.4
    r = dRes[(s, p)]
    R.assign("sens", s)
    R.assign("prec", p)
    R.assign("rb.fn.expr.mmseq", r.clMmseqExpr.fpExprMmseq)
    R("rb.mmseq.out <- readmmseq(mmseq_files=rb.fn.expr.mmseq, normalize=FALSE)")
    R("rb.expr.mmseq <- exp(rb.mmseq.out$log_mu)")
    R("rb.expr.mmseq[ is.na(rb.expr.mmseq) ] <- 0.0")

    annoFPs = set(r.simulatedTranscriptome.lFP)
    totals = set(R("row.names(rb.mmseq.out$log_mu)")) # total transcripts
    matches = set(r.clTranscriptome.knownSetTranscriptsMatch) # total TPs
    totalTPs = totals.intersection(matches)
    totalFPs = totals.difference(totalTPs)
    newFPs = totalFPs.difference(annoFPs)

    annoFPness = [ s in annoFPs for s in R("row.names(rb.mmseq.out$log_mu)") ]
    newFPness = [ s in newFPs for s in R("row.names(rb.mmseq.out$log_mu)") ]
    R.assign("annoFPs", annoFPness)
    R.assign("newFPs", newFPness)
    R("annoFPs <- unlist(annoFPs)")
    R("newFPs <- unlist(newFPs)")
    R("""
        plot(density(log(rb.expr.mmseq[annoFPs])),col="blue",lwd=2, main=paste("s=", sens, ", p=", prec), xlab="RABT FPs", ylab="")
        lines(density(log(rb.expr.mmseq[newFPs])),lwd=2,col="red")

    """)
    R("dev.off()")    

    # Density plot of est. Cufflinks, Oases, Anno, RABT FPs
    # Cufflinks
    R.assign("cuff.fn.expr.mmseq", cuffRes.clMmseqExpr.fpExprMmseq)
    R("cuff.mmseq.out <- readmmseq(mmseq_files=cuff.fn.expr.mmseq, normalize=FALSE)")
    R("cuff.expr.mmseq <- exp(cuff.mmseq.out$log_mu)")
    R("cuff.expr.mmseq[ is.na(cuff.expr.mmseq) ] <- 0.0")
    R.assign("cuff.mTranscriptTP", cuffRes.clTranscriptome.getTranscriptTPness([ s for s in R("row.names(cuff.mmseq.out$log_mu)") ]))
    R("cuff.mTranscriptTP <- unlist(cuff.mTranscriptTP)")
    R("cuffFPs <- cuff.expr.mmseq[!cuff.mTranscriptTP]")

    # Oases
    R.assign("oases.fn.expr.mmseq", oasesRes.oaMmseqExpr.fpExprMmseq)
    R("oases.mmseq.out <- readmmseq(mmseq_files=oases.fn.expr.mmseq, normalize=FALSE)")
    R("oases.expr.mmseq <- exp(oases.mmseq.out$log_mu)")
    R("oases.expr.mmseq[ is.na(oases.expr.mmseq) ] <- 0.0")
    R.assign("oases.mTranscriptTP", oasesRes.oasesTranscriptome.getTranscriptTPness([ s for s in R("row.names(oases.mmseq.out$log_mu)") ]))
    R("oases.mTranscriptTP <- unlist(oases.mTranscriptTP)")
    R("oasesFPs <- oases.expr.mmseq[!oases.mTranscriptTP]")

    R("""
        pdf("_FPDen.pdf")
        #pdf("_logFPDen.pdf")
        par(mfrow=c(5,5), mar=c(1,1,2,1),oma=c(4,1,2,0),mgp=c(3, 0.5, 0))
        #par(mfrow=c(5,5), mar=c(5.1,4.1,4.1,2.1), oma=c(2,2,2,2), mgp=c(3, 0.5, 0))
    """)
    for (s, p) in dRes.keys():
            r = dRes[(s, p)]
            R.assign("sens", s)
            R.assign("prec", p)

            # Anno
            R.assign("an.fn.expr.mmseq", r.mmseqExpr.fpExprMmseq)
            R("an.mmseq.out <- readmmseq(mmseq_files=an.fn.expr.mmseq, normalize=FALSE)")
            R("an.expr.mmseq <- exp(an.mmseq.out$log_mu)")
            R("an.expr.mmseq[ is.na(an.expr.mmseq) ] <- 0.0")
            # Calculate matching vector of true expression values
            R.assign("an.mTranscriptTP", r.clTranscriptome.getTranscriptTPness([ s for s in R("row.names(an.mmseq.out$log_mu)") ]))
            R("an.mTranscriptTP <- unlist(an.mTranscriptTP)")
            R("annoFPs <- an.expr.mmseq[!an.mTranscriptTP]")

            # RABT
            R.assign("rb.fn.expr.mmseq", r.clMmseqExpr.fpExprMmseq)
            R("rb.mmseq.out <- readmmseq(mmseq_files=rb.fn.expr.mmseq, normalize=FALSE)")
            R("rb.expr.mmseq <- exp(rb.mmseq.out$log_mu)")
            R("rb.expr.mmseq[ is.na(rb.expr.mmseq) ] <- 0.0")
            R.assign("rb.mTranscriptTP", r.clTranscriptome.getTranscriptTPness([ s for s in R("row.names(rb.mmseq.out$log_mu)") ]))
            R("rb.mTranscriptTP <- unlist(rb.mTranscriptTP)")
            R("rabtFPs <- rb.expr.mmseq[!rb.mTranscriptTP]")

            R("""
                plot(density(cuffFPs),col="blue",lwd=2, main=paste("s=", sens, ", p=", prec), xlab="est FPs", ylab="")
                lines(density(oasesFPs),lwd=2,col="yellow")
                lines(density(annoFPs),lwd=2,col="green")
                lines(density(rabtFPs),lwd=2,col="red")

                # plot(density(log(cuffFPs)),col="blue",lwd=2, main=paste("s=", sens, ", p=", prec), xlab="est FPs", ylab="")
                # lines(density(log(oasesFPs)),lwd=2,col="yellow")
                # lines(density(log(annoFPs)),lwd=2,col="green")
                # lines(density(log(rabtFPs)),lwd=2,col="red")
            """)
    R("dev.off()")

    # Plot annotations and RABT for TP+FP and TP only transcripts    
    '''
    R.pdf("_fig3a_SensPrec.pdf", width=6, height=10)
    R("par(mfrow=c(2,1))")
    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n", xlab="sensitivity", ylab="precision", main="Expression correlation of TP transcripts")
    R.text(plotAnSensitivity, plotAnPrecision, labels=plotAnExpressionTP, cex=1.0, col=1)
    R.legend("bottomright", R.c("annotations"), pch=R.c(0), col=R.c(1))
    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n", xlab="sensitivity", ylab="precision", main="Expression correlation of TP transcripts")
    R.arrows(rpy2.robjects.vectors.FloatVector(plotAnSensitivity), rpy2.robjects.vectors.FloatVector(plotAnPrecision), 
             rpy2.robjects.vectors.FloatVector(plotRbSensitivity), rpy2.robjects.vectors.FloatVector(plotRbPrecision),
             length=0.05, col="black")
    # R.text(cuffRes.clTranscriptome.getSensitivity(), cuffRes.clTranscriptome.getPrecision(), labels=plotClExpressionTP, cex=1.0, col=2)
    # R.text(plotRbSensitivity, plotRbPrecision, labels=plotRbExpressionTP, cex=1.0, col=3)
    # R.text(oasesRes.oasesTranscriptome.getSensitivity(), oasesRes.oasesTranscriptome.getPrecision(), labels=plotOaExpressionTP, cex=1.0, col=4)
    # R.legend("bottomright", R.c("Cufflinks", "RABT+Cufflinks", "Oases"), pch=R.c(0,0,0), col=R.c(2,3,4))
    R("dev.off()")

    R.pdf("_fig3b_SensPrec.pdf", width=6, height=10)
    R("par(mfrow=c(2,1))")
    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n",xlab="sensitivity", ylab="precision", main="Normalized mean FP expression")
    R.text(plotAnSensitivity, plotAnPrecision, labels=plotAnExpressionFP, cex=1.0, col=1)
    R.legend("bottomright", R.c("annotations"), pch=R.c(0), col=R.c(1))
    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n",xlab="sensitivity", ylab="precision", main="Normalized mean FP expression")
    R.arrows(rpy2.robjects.vectors.FloatVector(plotAnSensitivity), rpy2.robjects.vectors.FloatVector(plotAnPrecision), 
             rpy2.robjects.vectors.FloatVector(plotRbSensitivity), rpy2.robjects.vectors.FloatVector(plotRbPrecision),
             length=0.05, col="black")
    # R.text(cuffRes.clTranscriptome.getSensitivity(), cuffRes.clTranscriptome.getPrecision(), labels=plotClExpressionFP, cex=1.0, col=2)
    # R.text(plotRbSensitivity, plotRbPrecision, labels=plotRbExpressionFP, cex=1.0, col=3)
    # R.text(oasesRes.oasesTranscriptome.getSensitivity(), oasesRes.oasesTranscriptome.getPrecision(), labels=plotOaExpressionFP, cex=1.0, col=4)
    # R.legend("bottomright", R.c("Cufflinks", "RABT+Cufflinks", "Oases"), pch=R.c(0,0,0), col=R.c(2,3,4))
    R("dev.off()")
    #ip.embed()
    '''
    R.pdf("_fig_SensPrec.pdf", width=12, height=10)
    R("par(mfrow=c(2,2))")
    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n", xlab="sensitivity", ylab="precision", main="Expression correlation of TP transcripts")
    R.text(plotAnSensitivity, plotAnPrecision, labels=plotAnExpressionTP, cex=1.0, col=1)
    R.legend("bottomright", R.c("annotations"), pch=R.c(0), col=R.c(1))

    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n",xlab="sensitivity", ylab="precision", main="Normalized mean FP expression")
    R.text(plotAnSensitivity, plotAnPrecision, labels=plotAnExpressionFP, cex=1.0, col=1)
    R.legend("bottomright", R.c("annotations"), pch=R.c(0), col=R.c(1))

    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n", xlab="sensitivity", ylab="precision", main="Expression correlation of TP transcripts")
    R.arrows(rpy2.robjects.vectors.FloatVector(plotAnSensitivity), rpy2.robjects.vectors.FloatVector(plotAnPrecision), 
             rpy2.robjects.vectors.FloatVector(plotRbSensitivity), rpy2.robjects.vectors.FloatVector(plotRbPrecision),
             length=0.05, col="black")

    R.plot(x=R.c(0,1), y=R.c(0,1), xlim=R.c(0,1.05), ylim=R.c(0.05,0.5), type="n",xlab="sensitivity", ylab="precision", main="Normalized mean FP expression")
    R.arrows(rpy2.robjects.vectors.FloatVector(plotAnSensitivity), rpy2.robjects.vectors.FloatVector(plotAnPrecision), 
             rpy2.robjects.vectors.FloatVector(plotRbSensitivity), rpy2.robjects.vectors.FloatVector(plotRbPrecision),
             length=0.05, col="black")
    R("dev.off()")