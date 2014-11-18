#!/usr/bin/env python

import sys; sys.path.append("core")
import ductape as tp
from steps import *

#return np.linspace(pMin, pMax, num)
class pipelineCufflinks():
    def __init__(self, transcriptPool, simulatedReads, tophatBowtieIndex, tophatAlign, coverage, sensitivity, precision,  RABTassembly):
        self.transcriptPool = transcriptPool
        self.simulatedReads = simulatedReads
        self.simulatedTranscriptome = SimulatedTranscriptome2(transcriptPool=self.transcriptPool, sensitivity=sensitivity, precision=precision)#, recompute=True)
        # Estimate expression from annotated transcripts
        self.bowtieIndex = BowtieIndex(simulatedTranscriptome=self.simulatedTranscriptome)#,recompute=True)
        self.bowtieAlign = BowtieAlign(simulatedReads=self.simulatedReads, bowtieIndex=self.bowtieIndex)#,recompute=True)
        self.hitsfile = Hitsfile(bowtieAlign=self.bowtieAlign,simulatedTranscriptome=self.simulatedTranscriptome)#,recompute=True)
        self.mmseqExpr = Mmseq(hitsfile=self.hitsfile)#,recompute=True)
        self.plotMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=self.simulatedTranscriptome,mmseqExpr=self.mmseqExpr,transcriptPool=self.transcriptPool)#,recompute=True)
        
        # Estimate transcriptome with Cufflinks/RABT
        self.tophatBowtieIndex = tophatBowtieIndex
        self.tophatAlign = tophatAlign
        self.cufflinks = Cufflinks(tophatAlign=self.tophatAlign, simulatedTranscriptome=self.simulatedTranscriptome, RABTassembly=RABTassembly)#, recompute=True)
        self.clTranscriptome = CufflinksTranscriptome2(cufflinks=self.cufflinks, simulatedTranscriptome=self.simulatedTranscriptome)#, recompute=True)
        # Estimate expression with Cufflinks/RABT transcripts
        self.clBowtieIndex = BowtieIndex(simulatedTranscriptome=self.clTranscriptome)#, recompute=True)
        self.clBowtieAlign = BowtieAlign(simulatedReads=self.simulatedReads,bowtieIndex=self.clBowtieIndex)#, recompute=True)
        self.clHitsfile = Hitsfile(bowtieAlign=self.clBowtieAlign,simulatedTranscriptome=self.clTranscriptome)#, recompute=True)
        self.clMmseqExpr = Mmseq(hitsfile=self.clHitsfile)#, recompute=True)
        self.plotClMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=self.clTranscriptome,mmseqExpr=self.clMmseqExpr,transcriptPool=self.transcriptPool)#, recompute=True)

class pipelineOases():
    def __init__(self, transcriptPool, simulatedReads):
        self.transcriptPool = transcriptPool
        self.simulatedReads = simulatedReads
        self.oases = Oases(simulatedReads = self.simulatedReads)#,recompute=True)
        self.oasesTranscriptome = OasesTranscriptome(oases=self.oases, transcriptPool=self.transcriptPool)#, recompute=True)
        self.oaBowtieIndex = BowtieIndex(simulatedTranscriptome=self.oasesTranscriptome)#,recompute=True)
        self.oaBowtieAlign = BowtieAlign(simulatedReads=self.simulatedReads, bowtieIndex=self.oaBowtieIndex)#,recompute=True)
        self.oaHitsfile = Hitsfile(bowtieAlign=self.oaBowtieAlign,simulatedTranscriptome=self.oasesTranscriptome,enstFormat=False)#,recompute=True)
        self.oaMmseqExpr = Mmseq(hitsfile=self.oaHitsfile)#,recompute=True)
        self.oaPlotMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=self.oasesTranscriptome,mmseqExpr=self.oaMmseqExpr,transcriptPool=self.transcriptPool)#,recompute=True)
