#!/usr/bin/env python

import sys; sys.path.append("core")
import ductape as tp
from steps import *

class pipelineCoverage():
  def __init__(self, transcriptPool, coverage, organism):
      self.transcriptPool = transcriptPool
      
      ## Compute/load (sens, prec)-independant computations
      # self.transcriptPool = TranscriptPool2(flSubset=True, seedVal=1)
      self.simulatedReads = SimulatedReads2(nTotalFrags=len(self.transcriptPool.lTruthSet)*coverage, transcriptPool=self.transcriptPool)#, recompute=True)

      ## Annotations
      self.simulatedTranscriptome = SimulatedTranscriptome2(transcriptPool=transcriptPool, sensitivity=0.29, precision=0.16)#, recompute=True)
      #self.bowtieIndex = BowtieIndex(simulatedTranscriptome=simulatedTranscriptome)#,recompute=True)
      #self.bowtieAlign = BowtieAlign(simulatedReads=simulatedReads, bowtieIndex=bowtieIndex)#,recompute=True)
      #self.hitsfile = Hitsfile(bowtieAlign=bowtieAlign,simulatedTranscriptome=simulatedTranscriptome)#,recompute=True)
      #self.mmseqExpr = Mmseq(hitsfile=hitsfile)#,recompute=True)
      #self.plotMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=simulatedTranscriptome,mmseqExpr=mmseqExpr)#,recompute=True)

      ## Cufflinks: reconstruct transcriptome
      self.tophatBowtieIndex = TophatBowtieIndex(transcriptPool=self.transcriptPool)#,recompute=True)
      self.tophatAlign = TophatAlign(simulatedReads=self.simulatedReads, tophatBowtieIndex=self.tophatBowtieIndex, organism=organism)#,recompute=True)
      self.cufflinks = Cufflinks(tophatAlign=self.tophatAlign, simulatedTranscriptome=None, RABTassembly=False, organism=organism)
      #self.clTranscriptome1 = CufflinksTranscriptome1(cufflinks=self.cufflinks, simulatedTranscriptome=self.simulatedTranscriptome)
      self.clTranscriptome2 = CufflinksTranscriptome2(cufflinks=self.cufflinks, simulatedTranscriptome=self.simulatedTranscriptome)
      #self.clBowtieIndex = BowtieIndex(simulatedTranscriptome=self.clTranscriptome2)#,recompute=True)
      #self.clBowtieAlign = BowtieAlign(simulatedReads=self.simulatedReads, bowtieIndex=self.clBowtieIndex)#,recompute=True)
      #self.clHitsfile = Hitsfile(bowtieAlign=self.clBowtieAlign,simulatedTranscriptome=self.clTranscriptome2)#,recompute=True)
      #self.clMmseqExpr = Mmseq(hitsfile=self.clHitsfile)#,recompute=True)
      #self.clPlotMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=self.clTranscriptome2,mmseqExpr=self.clMmseqExpr)#,recompute=True)

      ## Oases
      self.oases = Oases(simulatedReads = self.simulatedReads)#,recompute=True)
      self.oasesTranscriptome = OasesTranscriptome(oases=self.oases, transcriptPool=self.transcriptPool)#, recompute=True)
      #self.oaBowtieIndex = BowtieIndex(simulatedTranscriptome=self.oasesTranscriptome)#,recompute=True)
      #self.oaBowtieAlign = BowtieAlign(simulatedReads=self.simulatedReads, bowtieIndex=self.oaBowtieIndex)#,recompute=True)
      #self.oaHitsfile = Hitsfile(bowtieAlign=self.oaBowtieAlign,simulatedTranscriptome=self.oasesTranscriptome,enstFormat=False)#,recompute=True)
      #self.oaMmseqExpr = Mmseq(hitsfile=self.oaHitsfile)#,recompute=True)
      #self.oaPlotMmseqPrediction = PlotMmseqPrediction2(simulatedTranscriptome=self.oasesTranscriptome,mmseqExpr=self.oaMmseqExpr)#,recompute=True)
