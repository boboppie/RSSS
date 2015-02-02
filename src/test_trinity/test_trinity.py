#!/usr/bin/env python

import sys
sys.path.append("core")
sys.path.append("coverage")
from steps import *

if __name__ == '__main__':
    #unittest.main(exit=True, verbosity=2)
    transcriptPool = TranscriptPool2(flSubset=True, seedVal=1, organism='worm', pActiveGene=0.1, pExpressedTrans=0.86)
    simulatedReads = SimulatedReads2(nTotalFrags=1E5, transcriptPool=transcriptPool)#, recompute=True)
    simulatedTranscriptome = SimulatedTranscriptome2(transcriptPool=transcriptPool, sensitivity=0.3, precision=0.2)#, recompute=True)
    trinity = Trinity(simulatedReads = simulatedReads)
    trinityTranscriptome = OasesTranscriptome(oases=trinity, transcriptPool=transcriptPool)#, recompute=True)
    print (trinityTranscriptome.getSensitivity(), trinityTranscriptome.getPrecision())

    #bowtieIndex = BowtieIndex(simulatedTranscriptome=simulatedTranscriptome)#,recompute=True)
    #bowtieAlign = BowtieAlign(simulatedReads=simulatedReads, bowtieIndex=bowtieIndex)#,recompute=True)
    #hitsfile = Hitsfile(bowtieAlign=bowtieAlign,simulatedTranscriptome=simulatedTranscriptome)#,recompute=True)
    #mmseqExpr = Mmseq(hitsfile=hitsfile)#,recompute=True)
    #plotMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=simulatedTranscriptome,mmseqExpr=mmseqExpr,transcriptPool=transcriptPool)#,recompute=True)
    #print "Correlation between predictions and ground truth: %f" % (plotMmseqPrediction.correlation,)


