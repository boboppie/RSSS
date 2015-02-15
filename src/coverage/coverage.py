#!/usr/bin/env python

import sys
sys.path.append("core")
sys.path.append("coverage")
import ductape as tp
from steps import *
from coverage_pipeline import *
import rpy2
import rpy2.robjects

if __name__ == '__main__':
    #unittest.main(exit=True, verbosity=2)
    organism = "human";pActiveGene=0.49;pExpressedTrans=0.39;title="Human"
    # organism = "mouse";pActiveGene=0.64;pExpressedTrans=0.61;title="Mouse"
    # organism = "worm";pActiveGene=0.4;pExpressedTrans=0.86;title="Worm"

    # Compute/load (sens, prec)-independant computations
    transcriptPool = TranscriptPool2(flSubset=True, seedVal=1, pActiveGene=pActiveGene, pExpressedTrans=pExpressedTrans, organism=organism)
    lTranscripts = list(Bio.SeqIO.parse(transcriptPool.fpTruthSet, "fasta"))
    sumTruthSetLen = sum(map(lambda x: len(x), lTranscripts))

    # Compute/load (sens, prec)-dependant computations
    #lCoverage = map(int, np.logspace(1, 3, num=7)) #30min?

    # Return numbers spaced evenly on a log scale
    # Ref http://docs.scipy.org/doc/numpy/reference/generated/numpy.logspace.html
    lCoverage = map(int, np.logspace(1, 4, num=10)) #~4h
    #coverage = #simulated fragments * 2 * 75 / ( sum_{t in 1..N_e} l_t ), where l_t is the length of transcript t
    lNewCoverage = map(lambda x: x*len(lTranscripts)*2*75/sumTruthSetLen, lCoverage)
    dRes = {}
    for coverage in lCoverage:
        dRes[coverage] = pipelineCoverage(transcriptPool=transcriptPool, coverage=coverage, organism=organism)

    # Calculate results
    R = rpy2.robjects.r
    plotClSens = [dRes[coverage].clTranscriptome2.getSensitivity() for coverage in lCoverage]
    plotClPrec = [dRes[coverage].clTranscriptome2.getPrecision() for coverage in lCoverage]
    plotOaSens = [dRes[coverage].oasesTranscriptome.getSensitivity() for coverage in lCoverage]
    plotOaPrec = [dRes[coverage].oasesTranscriptome.getPrecision() for coverage in lCoverage]
    plotTrSens = [dRes[coverage].trinityTranscriptome.getSensitivity() for coverage in lCoverage]
    plotTrPrec = [dRes[coverage].trinityTranscriptome.getPrecision() for coverage in lCoverage]
    rangeSensY = R.range(0, plotClSens, plotOaSens, plotTrSens)
    rangePrecY = R.range(0, plotClPrec, plotOaPrec, plotTrPrec)

    # ip.embed()
    cExfact = 1.5
    # Plot results
    R.pdf("_figCoverageVsReconstruction.pdf", width=6, height=9)    
    R("par(mfrow=c(2,1))")
    R("par(mar=c(3.9,4.5,4,2) + 0.1, cex=1.1)")
    R.plot(lNewCoverage, plotClSens, type="b", main=title, xlab="Coverage", ylab="Sensitivity", log="x", col=2, ylim=R.c(0,1), pch=2, cex=cExfact)
    R.points(lNewCoverage, plotOaSens, type="b", col=4, pch=4, cex=cExfact)
    R.points(lNewCoverage, plotTrSens, type="b", col=5, pch=4, cex=cExfact)
    # R.abline(h=R.c(0.2, 0.4, 0.6, 0.8), col=1, lty=3, lwd=cExfact)
    R.par(bg="white")
    R.legend("topleft", R.c("Cufflinks", "Oases", "Trinity"), pch=R.c(2,4), col=R.c(2, 4, 5), lty=R.c(0, 0, 0), lwd=R.c(1.2, 1.2, 1.2))
    R("par(mar=c(5.3,4.5,2.5,2) + 0.1)")
    R.plot(lNewCoverage, plotClPrec, type="b", main="", xlab="Coverage", ylab="Precision", log="x", col=2, ylim=R.c(0,1), pch=2, cex=cExfact)
    R.points(lNewCoverage, plotOaPrec, type="b", col=4, pch=4, cex=cExfact)
    R.points(lNewCoverage, plotTrPrec, type="b", col=5, pch=4, cex=cExfact)
    # R.abline(h=R.c(0.2, 0.4, 0.6, 0.8), col=1, lty=3, lwd=cExfact)
    R.par(bg="white")
    R.legend("topleft", R.c("Cufflinks", "Oases", "Trinity"), pch=R.c(2,4), col=R.c(2, 4, 5), lty=R.c(0, 0, 0), lwd=R.c(1.2, 1.2, 1.2))
    R("dev.off()")

    # for mouse plot
    '''
    cExfact = 1.5
    # Plot results
    R.pdf("_figCoverageVsReconstruction.pdf", width=14, height=6)
    R("par(mfrow=c(1,2), cex=1.1, oma=c(0,0,2,0))")
    R("par(mar=c(5,4.5,2,2) + 0.1)")

    #tp.log(lnewCoverage)
    R.plot(lnewCoverage, plotClSens, type="b", xlab="Coverage", ylab="Sensitivity", log="x", col=2, ylim=R.c(0,1), pch=2, cex=cExfact)
    R.points(lnewCoverage, plotOaSens, type="b", col=4, pch=4, cex=cExfact)
    # R.abline(h=R.c(0.2, 0.4, 0.6, 0.8), col=1, lty=3, lwd=cExfact)
    R.par(bg="white")
    R.legend("topleft", R.c("Cufflinks", "Oases"), pch=R.c(2,4), col=R.c(2, 4), lty=R.c(0, 0), lwd=R.c(1.2, 1.2))
    R("par(mar=c(5,4,2,2) + 0.1)")
    R.plot(lnewCoverage, plotClPrec, type="b", xlab="Coverage", ylab="Precision", log="x", col=2, ylim=R.c(0,1), pch=2, cex=cExfact)
    R.points(lnewCoverage, plotOaPrec, type="b", col=4, pch=4, cex=cExfact)
    # R.abline(h=R.c(0.2, 0.4, 0.6, 0.8), col=1, lty=3, lwd=cExfact)
    R.par(bg="white")
    R.legend("topleft", R.c("Cufflinks", "Oases"), pch=R.c(2,4), col=R.c(2, 4), lty=R.c(0, 0), lwd=R.c(1.2, 1.2))
    R("title(main='Mouse',outer=T)")
    R("dev.off()")
    '''
