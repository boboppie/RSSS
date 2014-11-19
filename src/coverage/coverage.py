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
    organism = "human";pActiveGene=0.49;pExpressedTrans=0.39
    # organism = "mouse";pActiveGene=0.64;pExpressedTrans=0.61
    # organism = "worm";pActiveGene=0.4;pExpressedTrans=0.86

    # Compute/load (sens, prec)-independant computations
    transcriptPool = TranscriptPool2(flSubset=True, seedVal=1, pActiveGene=pActiveGene, pExpressedTrans=pExpressedTrans, organism=organism)

    # Compute/load (sens, prec)-dependant computations
    #lCoverage = map(int, np.logspace(1, 3, num=7)) #30min?

    # Return numbers spaced evenly on a log scale
    # Ref http://docs.scipy.org/doc/numpy/reference/generated/numpy.logspace.html
    lCoverage = map(int, np.logspace(1, 4, num=10)) #~4h
    dRes = {}
    for coverage in lCoverage:
        dRes[coverage] = pipelineCoverage(transcriptPool=transcriptPool, coverage=coverage, organism=organism)

    # Calculate results
    R = rpy2.robjects.r
    plotClSens = [dRes[coverage].clTranscriptome2.getSensitivity() for coverage in lCoverage]
    plotClPrec = [dRes[coverage].clTranscriptome2.getPrecision() for coverage in lCoverage]
    plotOaSens = [dRes[coverage].oasesTranscriptome.getSensitivity() for coverage in lCoverage]
    plotOaPrec = [dRes[coverage].oasesTranscriptome.getPrecision() for coverage in lCoverage]
    rangeSensY = R.range(0, plotClSens, plotOaSens)
    rangePrecY = R.range(0, plotClPrec, plotOaPrec)

    # ip.embed()
    cExfact = 1.5
    # Plot results
    R.pdf("_figCoverageVsReconstruction.pdf", width=6, height=8)    
    R("par(mfrow=c(2,1))")
    R.plot(lCoverage, plotClSens, type="b", main="Human", xlab="coverage", ylab="sensitivity", log="x", col=2, ylim=R.c(0,1), pch=2, cex=cExfact)
    R.points(lCoverage, plotOaSens, type="b", col=4, pch=4, cex=cExfact)
    # R.abline(h=R.c(0.2, 0.4, 0.6, 0.8), col=1, lty=3, lwd=cExfact)
    R.par(bg="white")
    R.legend("topleft", R.c("Cufflinks", "Oases"), pch=R.c(2,4), col=R.c(2, 4), lty=R.c(0, 0), lwd=R.c(0,0))
    R.plot(lCoverage, plotClPrec, type="b", main="", xlab="coverage", ylab="precision", log="x", col=2, ylim=R.c(0,1), pch=2, cex=cExfact)
    R.points(lCoverage, plotOaPrec, type="b", col=4, pch=4, cex=cExfact)
    # R.abline(h=R.c(0.2, 0.4, 0.6, 0.8), col=1, lty=3, lwd=cExfact)
    R.par(bg="white")
    R.legend("topleft", R.c("Cufflinks", "Oases"), pch=R.c(2,4), col=R.c(2, 4), lty=R.c(0, 0), lwd=R.c(0, 0))
    R("dev.off()")
