#!/usr/bin/env python

import cPickle
import os
import hashlib
import gzip
import random
from pprint import pprint
from collections import OrderedDict
import tempfile
import subprocess
import re
import shutil
import unittest
import multiprocessing
import operator
import numpy as np
import numpy.random
import scipy.stats as sps
import matplotlib; matplotlib.use("PDF")
import matplotlib.pyplot as plt
import Bio
import Bio.SeqIO
import IPython as ip
import rpy2
import rpy2.robjects
import ductape as tp

def grepGtfForTranscriptIDs(fpInp, fpOut, lTranscriptID):
    sTranscriptID = set(lTranscriptID)
    fhInp = open(fpInp, "r")
    fhOut = open(fpOut, "w")
    for l in fhInp: 
        transcriptid =  l.split("\t")[8].split(";")[1].split(" ")[2][1:-1]
        if (transcriptid in sTranscriptID):
            fhOut.write(l)
    fhInp.close()
    fhOut.close()

def grepFasta(fpInp, fpOut, fuFilter):
    fhInp = open(fpInp, "r")
    fhOut = open(fpOut, "w")
    for record in Bio.SeqIO.parse(fhInp, "fasta"):
        if fuFilter(record):
            Bio.SeqIO.write(record, fhOut, "fasta")
    fhInp.close()
    fhOut.close()

def grepFastaForTranscriptID(fpInp, fpOut, lTranscriptID):
    sTranscriptID = set(lTranscriptID)
    fhInp = open(fpInp, "r")
    fhOut = open(fpOut, "w")
    fMatch = False
    for l in fhInp:
        if (l[0] == '>'):
            if (l.split()[0][1:] in sTranscriptID):
                fMatch = True
            else:
                fMatch = False
        if (fMatch): fhOut.write(l)
    fhInp.close()
    fhOut.close()

def plotCoverageVsClReconstruction():
    class pipeline: pass
    sensitivity = 0.5
    precision = 1.0
    nTruthSet = 100
    transcriptPool = TranscriptPool() #tp.log("transcript pool size: %d, first transcript id: %s" % (len(transcriptPool.lTranscripts), transcriptPool.lTranscripts[0].id))
    simulatedTranscriptome = SimulatedTranscriptome(nTruthSet=nTruthSet,completeness=sensitivity,purity=precision,transcriptPool=transcriptPool) #tp.log(truthSet.iTranscriptPoolInTruthSet)
    #lCoverage = list(np.logspace(2, 3, 3)) # development range
    #lCoverage = list(np.logspace(1,5,9)) # sparse sampling over whole "reasonable" range
    lCoverage = list(np.logspace(1,3,9)) # demonstrates convergence
    #lCoverage = list(np.logspace(1,5,17)) # density same as above, for whole range
    dRes = {coverage: pipeline() for coverage in lCoverage}
    for (coverage, r) in dRes.iteritems():
        r.simulatedReads = SimulatedReads(nTotalFrags=nTruthSet*coverage,simulatedTranscriptome=simulatedTranscriptome,insertSizeMean=250,insertSizeSd=30)
        # Reconstruct transcriptome with Cufflinks
        r.tophatBowtieIndex = TophatBowtieIndex(transcriptPool=transcriptPool)
        r.tophatAlign = TophatAlign(simulatedReads=r.simulatedReads,tophatBowtieIndex=r.tophatBowtieIndex)
        r.cufflinks = Cufflinks(tophatAlign=r.tophatAlign)#,recompute=True)
        r.clTranscriptome = CufflinksTranscriptome(cufflinks=r.cufflinks,simulatedTranscriptome=simulatedTranscriptome)#,recompute=True)
        # Reconstruct transcriptome with Cufflinks+RABT
        r.cufflinksRABT = Cufflinks(tophatAlign=r.tophatAlign,RABTassembly=True)#,recompute=True)
        r.rbTranscriptome = CufflinksTranscriptome(cufflinks=r.cufflinksRABT, simulatedTranscriptome=simulatedTranscriptome)#,recompute=True)

    nClTP = [float(dRes[coverage].clTranscriptome.nTP) for coverage in lCoverage]
    nClFP = [float(dRes[coverage].clTranscriptome.nFP) for coverage in lCoverage]
    nClFN = [float(dRes[coverage].clTranscriptome.nFN) for coverage in lCoverage]
    lClSensitivity = list(np.array(nClTP) / (np.array(nClTP)+np.array(nClFN)))
    lClPrecision = list(np.array(nClTP) / (np.array(nClTP)+np.array(nClFP)))

    nRbTP = [float(dRes[coverage].rbTranscriptome.nTP) for coverage in lCoverage]
    nRbFP = [float(dRes[coverage].rbTranscriptome.nFP) for coverage in lCoverage]
    nRbFN = [float(dRes[coverage].rbTranscriptome.nFN) for coverage in lCoverage]
    lRbSensitivity = list(np.array(nRbTP) / (np.array(nRbTP)+np.array(nRbFN)))
    lRbPrecision = list(np.array(nRbTP) / (np.array(nRbTP)+np.array(nRbFP)))

    R = rpy2.robjects.r
    R.pdf("plotCoverageVsClReconstruction.pdf", width=7, height=7)
    R("par(mfrow=c(2,2))")
    # Plot TP & FP \n sensitivity & precision
    R.plot(lCoverage, nClTP, main="Cufflinks reconstruction sensitivity", xlab="mean no of fragments from a transcript", ylab="TP transcripts", log="x", ylim=R.c(0.0, max(np.max(nClTP), np.max(nRbTP))))
    R.points(lCoverage, nRbTP, pch=2)
    R.plot(lCoverage, nClFP, main="Cufflinks reconstruction precision", xlab="mean no of fragments from a transcript", ylab="FP transcripts", log="x", ylim=R.c(0.0, max(np.max(nClFP), np.max(nRbFP))))
    R.points(lCoverage, nRbFP, pch=2)
    R.plot(lCoverage, lClSensitivity, main="Cufflinks reconstruction sensitivity", xlab="mean no of fragments from a transcript", ylab="sensitivity", log="x", ylim=R.c(0.0, 1.0))
    R.points(lCoverage, lRbSensitivity, pch=2)
    R.plot(lCoverage, lClPrecision, main="Cufflinks reconstruction precision", xlab="mean no of fragments from a transcript", ylab="precision", log="x", ylim=R.c(0.0, 1.0))
    R.points(lCoverage, lRbPrecision, pch=2)

    #R.plot(lCoverage, lSensitivity, main="Cufflinks reconstruction sensitivity", xlab="mean no of fragments from a transcript", ylab="sensitivity", log="x", ylim=R.c(0.0, 1.0))
    #R.axis(1, at=R.c(10, 100, 1000, 10000, 100000), labels=R.c("10^1", "10^2","10^3","10^4", "10^5"))
    #R.axis(1, at=R.c(10, 100, 1000, 10000, 100000), labels=R.c("10^1", "10^2","10^3","10^4", "10^5"))
    R("dev.off()")

    """
    #ks = readGTF(simulatedTranscriptome.fpKnownSetGTF)
    ks = readGTF(simulatedTranscriptome.fpTruthSetGTF)
    kstid = set([ ks[i]['attribute']['transcript_id'] for i in range(len(ks)) ])
    cl = readGTF(dRes[100].clTranscriptome.fpKnownSetGTF)
    cltid = set([ cl[i]['attribute']['transcript_id'] for i in range(len(cl)) ])

    tp.log("KnownSet transcripts not in CL:")
    for kst in kstid: 
        if not(kst in cltid): print "Not in cltid: %s" % (kst,)

    tp.log("CL transcripts not in KnownSet:")
    for clt in cltid: 
        if not(clt in kstid): print "Not in kstid: %s" % (clt,)
    """

class TranscriptPool2(tp.NaivelyCachedComputation):
    def compute(self, flSubset=True, pActiveGene=0.49, pExpressedTrans=0.39, seedVal=1, organism="human"):

        self.rawDataDict = {'human': {'paSeqname': '21|22',
                              'pathGenomeRawFa': 'ensembl_human/release-66/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.66.dna.toplevel.fa.gz', 
                              'pathTranscriptomeRawFa': 'ensembl_human/release-66/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.66.cdna.all.fa.gz',
                              'pathTranscriptomeRawGtf': 'ensembl_human/release-66/gtf/homo_sapiens/Homo_sapiens.GRCh37.66.gtf.gz',
                              'fastagrep': 'supercontig|GRCh37:Y:1:10000|GRCh37:[^1-9XMY]',
                              'build': 'GRCh37'}, 
                            'mouse': {'paSeqname': '2|3',
                                      'pathGenomeRawFa': 'ensembl_mouse/release-66/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.66.dna.toplevel.fa.gz', 
                                      'pathTranscriptomeRawFa': 'ensembl_mouse/release-66/fasta/mus_musculus/cdna/Mus_musculus.NCBIM37.66.cdna.all.fa.gz',
                                      'pathTranscriptomeRawGtf': 'ensembl_mouse/release-66/gtf/mus_musculus/Mus_musculus.NCBIM37.66.gtf.gz',
                                      'fastagrep': 'scaffold|PATCH',
                                      'build': 'NCBIM37'},
                            'worm': {'paSeqname': 'III|V',
                                      'pathGenomeRawFa': 'ensembl_worm/release-66/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WS220.66.dna.toplevel.fa.gz', 
                                      'pathTranscriptomeRawFa': 'ensembl_worm/release-66/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WS220.66.cdna.all.fa.gz',
                                      'pathTranscriptomeRawGtf': 'ensembl_worm/release-66/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WS220.66.gtf.gz',
                                      'fastagrep': 'WS220:[^IVXM]',
                                      'build': 'WS220'},          
                            }

        self.organism = organism
        if self.organism not in self.rawDataDict:
            tp.log("Can not find organism %s in rawDataDict, use human as default" % self.organism)
            self.organism = "human"

        # Set random seed(s) for replication
        np.random.seed(seedVal)
        random.seed(seedVal)
        R = rpy2.robjects.r
        R("set.seed(%d)" % (seedVal,))

        # Make subset
        if (flSubset):
            self.paSeqname = self.rawDataDict[self.organism]['paSeqname']
        else:
            self.paSeqname = "*"

        # Preprocess Ensembl genes and transcripts
        self.unzipAndFilter()

        # List of (ENST, ENSG, expr) of all transcripts (gene/transcript pool)
        cmd = 'grep "^>" %s | perl -p -a -e \'$_ = substr($F[0], 1) . " " . substr($F[3], 5) . "\n";\'' % (self.fpTranscriptomeSubFa,)
        self.lTranGeneEx = [ l.split() for l in subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).split("\n") if l]
        for tge in self.lTranGeneEx: tge.append(0.0)

        # Generate a set of active genes (names) from the gene pool
        self.pActiveGene = pActiveGene # f_g
        self.sGene = set([transcript[1] for transcript in self.lTranGeneEx])
        self.sActiveGene = set(random.sample(self.sGene, int(round(self.pActiveGene * len(self.sGene)))))

        # Generate a list of expressed transcripts (TruthSet) from the transcript pool
        self.pExpressedTrans = pExpressedTrans # f_t
        self.sTrans = set([transcript[0] for transcript in self.lTranGeneEx])
        # at lease 1 transcript of each active gene has to be expressed
        # put one transcript from each active gene to the expressed transcript pool
        self.sActiveGeneTrans = set()

        # create a dict key - active gene, value - list of transcripts of the key
        geneTransDict = dict()

        for g in self.sActiveGene:
            for tge in self.lTranGeneEx:
                if g == tge[1]:
                    self.sActiveGeneTrans.add(tge[0])
                    if g in geneTransDict:
                        geneTransDict[g].append(tge[0])
                    else:
                        geneTransDict[g] = [tge[0]]
       
        # create a set containing one transcript of each active gene
        sTranBase = set()
        
        for g, lT in geneTransDict.items():
            sTranBase.add(random.sample(lT, 1)[0])
                      
        # a new set which will be used to sample               
        sTranExtraPool = self.sActiveGeneTrans.difference(sTranBase)   
        sampleSize = int(round(self.pExpressedTrans * len(self.sActiveGeneTrans) - len(sTranBase)))
        if sampleSize < 0:
            raise Exception("f_t is too small!")
        sTranExtra = set(random.sample(sTranExtraPool, sampleSize))

        # the TruthSet
        self.sExpressedTrans = sTranBase.union(sTranExtra)
        self.lTruthSet = list(self.sExpressedTrans)
        self.lUnexpressedSet = list(self.sTrans.difference(self.sExpressedTrans))

        for tge in self.lTranGeneEx:
            if (tge[0] in self.lTruthSet): tge[2] = np.random.gamma(shape = 1.2, scale = 1 / 0.001)

        # Write a .fa of turth set transcripts
        self.makeTruthSet()

    def unzipAndFilter(self):
        tp.log("Unzipping raw genome .fasta")
        self.fnGenomeRawFa = "genome.raw.fa"
        self.fpGenomeRawFa = self.fp("genome.raw.fa")
        tp.run("gunzip -c '%(inp)s' > '%(out)s'" % {
                "inp": datafp(self.rawDataDict[self.organism]['pathGenomeRawFa']),
                "out": self.fpGenomeRawFa
                })

        tp.log("Filtering alternative haplotypes/supercontigs from genome .fasta")
        self.fnGenomeHapFa = "genome.hap.fa"
        self.fpGenomeHapFa = self.fp("genome.hap.fa")

        tp.run("fastagrep.sh -v '%(pattern)s' %(inp)s > %(out)s" % {
                "pattern" : self.rawDataDict[self.organism]['fastagrep'],
                "inp": self.fnGenomeRawFa,
                "out": self.fnGenomeHapFa
                }, work_dir = self.fp(""))

        tp.log("Filtering sequence subset from the genome .fasta")
        self.fnGenomeSubFa = "genome.sub.fa"
        self.fpGenomeSubFa = self.fp("genome.sub.fa")
        tp.run("fastagrep.sh '%(pat)s' '%(inp)s' > '%(out)s'" % {
                "pat": "^>(%s)" % (self.paSeqname,),
                "inp": self.fnGenomeHapFa,
                "out": self.fnGenomeSubFa
                }, work_dir = self.fp(""))

        tp.log("Unzipping raw transcriptome .fasta")
        self.fnTranscriptomeRawFa = "transcriptome.raw.fa"
        self.fpTranscriptomeRawFa = self.fp("transcriptome.raw.fa")
        tp.run("gunzip -c '%(inp)s' > '%(out)s'" % {
                "inp": datafp(self.rawDataDict[self.organism]['pathTranscriptomeRawFa']),
                "out": self.fpTranscriptomeRawFa
                })

        tp.log("Filtering alternative haplotypes/supercontigs from transcriptome .fasta")
        self.fnTranscriptomeHapFa = "transcriptome.hap.fa"
        self.fpTranscriptomeHapFa = self.fp("transcriptome.hap.fa")
        tp.run("fastagrep.sh -v '%(pattern)s' %(inp)s > %(out)s" % {
                "pattern" : self.rawDataDict[self.organism]['fastagrep'],
                "inp": self.fnTranscriptomeRawFa,
                "out": self.fnTranscriptomeHapFa
                }, work_dir = self.fp(""))

        tp.log("Filtering sequence subset (and really short transcripts) from the transcriptome .fasta")
        self.fnTranscriptomeSubFa = "transcriptome.sub.fa"
        self.fpTranscriptomeSubFa = self.fp("transcriptome.sub.fa")
        prog = re.compile("chromosome:%s:(%s)" % (self.rawDataDict[self.organism]['build'], self.paSeqname,))
        grepFasta(fpInp = self.fpTranscriptomeHapFa,
                  fpOut = self.fpTranscriptomeSubFa,
                  fuFilter = lambda record: prog.search(record.description) and len(record.seq) > 75
                  )
        """
        tp.run("fastagrep.sh '%(pat)s' %(inp)s > %(out)s" % {
                "pat": "chromosome:GRCh37:(%s)" % (self.paSeqname,),
                "inp": self.fnTranscriptomeHapFa,
                "out": self.fnTranscriptomeSubFa
                }, work_dir = self.fp(""))
        """

        tp.log("Unzipping raw transcriptome .gtf")
        self.fpTranscriptomeRawGtf = self.fp("transcriptome.raw.gtf")
        tp.run("gunzip -c %(inp)s > %(out)s" % {
                "inp": datafp(self.rawDataDict[self.organism]['pathTranscriptomeRawGtf']),
                "out": self.fpTranscriptomeRawGtf
                })

        tp.log("Filtering alternative haplotypes/supercontigs from transcriptome .gtf")
        self.fpTranscriptomeHapGtf = self.fp("transcriptome.hap.gtf")
        grepGtfForTranscriptIDs(
            fpInp = self.fpTranscriptomeRawGtf, 
            fpOut = self.fpTranscriptomeHapGtf,
            lTranscriptID = subprocess.check_output("grep '^>' %s | awk '{ print substr($1, 2) }'" % (self.fp("transcriptome.hap.fa"),), stderr=subprocess.STDOUT, shell=True).split()
            )

        tp.log("Filtering sequence subset from the transcriptome .gtf")
        self.fpTranscriptomeSubGtf = self.fp("transcriptome.sub.gtf")
        grepGtfForTranscriptIDs(
            fpInp = self.fpTranscriptomeHapGtf, 
            fpOut = self.fpTranscriptomeSubGtf,
            lTranscriptID = subprocess.check_output("grep '^>' %s | awk '{ print substr($1, 2) }'" % (self.fp("transcriptome.sub.fa"),), stderr=subprocess.STDOUT, shell=True).split()
            )

        self.fpGenomeFa = self.fpGenomeSubFa

    def makeTruthSet(self):
        # Generate .fa and .gtf for expressed and active transcripts
        self.fpTruthSet = self.fp("truthSet.fa") 
        grepFastaForTranscriptID(fpInp=self.fpTranscriptomeSubFa, fpOut=self.fpTruthSet, lTranscriptID=self.lTruthSet)
        self.fpTruthSetGTF = self.fp("truthSet.gtf") 
        grepGtfForTranscriptIDs(fpInp=self.fpTranscriptomeSubGtf, fpOut=self.fpTruthSetGTF, lTranscriptID=self.lTruthSet)

        # Filter expression for truth set
        self.dTruthSetExpression = self.getExpressionAsDict(self.lTruthSet)

        # Calculate expression levels for truth set and known set
        tp.logList("TruthSet expression levels", list(self.dTruthSetExpression.values()))

        # Plot TruthSet generated expression levels
        R = rpy2.robjects.r
        self.fpPlotTruthSetExpression = self.fp("plotTruthSetExpression.pdf")
        R.pdf(file=self.fpPlotTruthSetExpression, width=6, height=6)
        R.hist(x=R.unlist(list(self.dTruthSetExpression.values())), main="plotTruthSetExpression.pdf", xlab="expression", breaks=40)
        R("dev.off()")
        tp.log("Wrote TruthSet expression level plot: %s" % (self.fpPlotTruthSetExpression,))

    def getExpressionAsDict(self, lTranscript):
        #non-redundant
        dTransEx = { l[0] : l[2] for l in self.lTranGeneEx }
        od = OrderedDict()
        for tid in lTranscript:
           od[tid] = dTransEx.get(tid, 0.0)
        return od

    def getExpressionAsList(self, lTranscript):
        dTransEx = { l[0] : l[2] for l in self.lTranGeneEx }
        return [ dTransEx.get(tid, 0.0) for tid in lTranscript ]
  
    def getExpressionById(self, transcriptId):
        dTransEx = { l[0] : l[2] for l in self.lTranGeneEx }
        return dTransEx.get(transcriptId, 0.0)


class SimulatedReads2(tp.NaivelyCachedComputation):
    def compute(self, nTotalFrags, transcriptPool, insertSizeMean=250, insertSizeSd=30, readLength=75):
        """
        Simulate reads from a TruthSet
            nTotalFrags = no of total read pairs
        Result:
            self.fpReads is the file with the fasta reads
        """
        self.nTotalFrags = nTotalFrags
        self.transcriptPool = transcriptPool
        self.insertSizeMean = insertSizeMean
        self.insertSizeSd = insertSizeSd

        # load transcripts from truth set fasta
        lTranscripts = list(Bio.SeqIO.parse(transcriptPool.fpTruthSet, "fasta"))
        # generate simulated expression values
        self.dExpression = transcriptPool.dTruthSetExpression

        # calculate total transcriptome size in bp
        transcriptomeSize = sum(map(lambda transcript: len(transcript.seq), lTranscripts))
        tp.log("total size of the transcriptome: %d" % (transcriptomeSize,))

        self.fnReads1 = "reads1.fq"
        self.fpReads1 = self.fp(self.fnReads1)
        self.fnReads2 = "reads2.fq"
        self.fpReads2 = self.fp(self.fnReads2)

        SimulatedReads.simulateTranscriptomeReads(lTranscripts, self.dExpression, nTotalFrags, 
                                   read_length=readLength, insert_size_mean=insertSizeMean, insert_size_sd=insertSizeSd,
                                   fn_output1 = self.fpReads1, fn_output2 = self.fpReads2, v=True,
                                   )

class SimulatedReads(tp.NaivelyCachedComputation):
    def compute(self, nTotalFrags, simulatedTranscriptome, insertSizeMean=250, insertSizeSd=30, readLength=75):
        """
        Simulate reads from a TruthSet
            nTotalFrags = no of total read pairs
        Result:
            self.fpReads is the file with the fasta reads
        """
        self.insertSizeMean = insertSizeMean
        self.insertSizeSd = insertSizeSd

        self.simulatedTranscriptome = simulatedTranscriptome
        # load transcripts from truth set fasta
        lTranscripts = list(Bio.SeqIO.parse(simulatedTranscriptome.fpTruthSet, "fasta"))
        # generate simulated expression values
        self.lExpression = simulatedTranscriptome.lTruthSetExpression

        # calculate total transcriptome size in bp
        transcriptomeSize = sum(map(lambda transcript: len(transcript.seq), lTranscripts))
        tp.log("total size of the transcriptome: %d" % (transcriptomeSize,))

        self.fnReads1 = "reads1.fq"
        self.fpReads1 = self.fp(self.fnReads1)
        self.fnReads2 = "reads2.fq"
        self.fpReads2 = self.fp(self.fnReads2)

        SimulatedReads.simulateTranscriptomeReads(lTranscripts, self.lExpression, nTotalFrags, 
                                   read_length=readLength, insert_size_mean=insertSizeMean, insert_size_sd=insertSizeSd,
                                   fn_output1 = self.fpReads1, fn_output2 = self.fpReads2, v=True,
                                   )

    @staticmethod
    def simulateTranscriptomeReads(l_transcript, d_expression, n_frags, read_length, insert_size_mean, insert_size_sd, fn_output1, fn_output2, addPoissonNoise=True, v=False):
        """
        Simulate reads from paired-end sequencing:
            l_transcript is a list of Seq objects representing the transcript sequences
            d_expression is a dictionary of transcrip ids and expression values proportional to the transcript concentration in the sample
            n_frags is the number of read pairs (fragments)
        """

        # allocate reads to individual transcripts

        """
        ### Test effective length calculation
        max_transcript_len = len(max(l_transcript, key=len))
        min_transcript_len = len(min(l_transcript, key=len))
        print "max length of a transcript: %d" % (max_transcript_len)
        print "min length of a transcript: %d" % (min_transcript_len)
        import random
        smpl_size = 10
        rand_smpl = [l_transcript[i] for i in sorted(random.sample(xrange(len(l_transcript)), smpl_size))]
        for transcript in rand_smpl:
            effective_len = sum(map(lambda x: (sps.norm.cdf(x + 0.5, loc=insert_size_mean, scale=insert_size_sd) - sps.norm.cdf(x - 0.5, loc=insert_size_mean, scale=insert_size_sd)) * (len(transcript) - x + 1), range(1, len(transcript))))
            print "transcript length %d - effective_len %d" % (len(transcript), effective_len,)
        ### Finish test effective length calculation
        """
        
        tot = 0.0
        effective_len_vec = []
        for transcript in l_transcript:
            expression = d_expression.get(transcript.id)
            # fragment (a.k.a insert) is from normal distribution
            effective_len = sum(map(lambda x: (sps.norm.cdf(x + 0.5, loc=insert_size_mean, scale=insert_size_sd) - sps.norm.cdf(x - 0.5, loc=insert_size_mean, scale=insert_size_sd)) * (len(transcript) - x + 1), range(1, len(transcript))))
            effective_len_vec.append(effective_len * expression)
            tot += effective_len * expression

        # determine the amount of reads to generate for a transcript
        p_vec = map(lambda x: x/float(tot), effective_len_vec)
        n_transcript_reads_vec = np.random.multinomial(n_frags, p_vec)

        # loop over all transcripts
        fh1 = open(fn_output1, "w")
        fh2 = open(fn_output2, "w")
        for (transcript, progress, n_transcript_reads) in zip(l_transcript, range(len(l_transcript)), n_transcript_reads_vec):
            # determine the amount of reads to generate for a transcript
            # n_transcript_reads = int(round((len(transcript) - read_length) * expression / float(tot) * n_frags))
        
            # add Poisson noise to transcript expression
            if (addPoissonNoise): n_transcript_reads = int(round(np.random.poisson(n_transcript_reads))) 
            if v: print "simulateTranscriptomeReads: %d / %d transcript_length=%d frags=%d" % (progress + 1, len(l_transcript), len(transcript.seq), n_transcript_reads)
            for i in range(n_transcript_reads):
                simulated_read = SimulatedReads.simulateTranscriptRead(transcript, read_length, insert_size_mean, insert_size_sd, "SIMREAD.%s.read%d" % (transcript.id, i + 1))
                if simulated_read is None: continue
                Bio.SeqIO.write(simulated_read[0], fh1, "fastq")
                Bio.SeqIO.write(simulated_read[1], fh2, "fastq")
        fh1.close()
        fh2.close()

    @staticmethod
    def simulateTranscriptRead(transcript, read_length=36, insert_size_mean=250, insert_size_sd=30, readid="NON_UNIQUE_READ_ID"):
        """
        Seq metadata is mapped to a .fasta "read description line" as:
            >read.id read.description
        """

        # The "insert" encompasses read1 and read2 as well as the unknown gap between them
        # Ref http://thegenomefactory.blogspot.co.uk/2013/08/paired-end-read-confusion-library.html
        # insert means fragment
        insert_size = int(round(np.random.normal(loc=insert_size_mean, scale=insert_size_sd)))
        if (insert_size < read_length + 1) or (insert_size >= len(transcript)): return(None)

        ##i_read1_first position: replace read_length to len(transcript) - insert_size?
        i_read1_first = np.random.randint(0, len(transcript) - insert_size) # read start position
        
        # if (insert_size < read_length + 1): return(None) # paired reads have to have at least one non-overlapping base
        if (i_read1_first + insert_size > len(transcript)): return(None) # discard insert sizes that aren't applicable to the transcript
        i_read1_endp1 = i_read1_first + read_length
        i_read2_first = i_read1_first + insert_size - read_length
        i_read2_endp1 = i_read2_first + read_length
        read1 = transcript[i_read1_first:i_read1_endp1]
        read2 = transcript[i_read2_first:i_read2_endp1].reverse_complement() # paired-end reads
        # either bowtie or bam2hits mangles readid's with a space in them
        # both read_id's have to be EXACTLY the same beyond a /1 /2 in the end (for e.g. samtools to work?)
        readid_head = "%s:Lfirst=%d:Lendp1=%d:Rfirst=%d:Rendp1=%d:#0" % (readid, i_read1_first, i_read1_endp1,i_read2_first, i_read2_endp1) 
        read1.id = "%s/1" % (readid_head, ) 
        read2.id = "%s/2" % (readid_head, )
        read1.description = "" 
        read2.description = ""
        read1.letter_annotations["phred_quality"] = [40] * len(read1.seq)
        read2.letter_annotations["phred_quality"] = [40] * len(read2.seq)
        return(read1, read2)
        """
        insert_size = int(round(np.random.normal(loc=insert_size_mean, scale=insert_size_sd)))
        if (insert_size > len(transcript)): return(None) # discard insert sizes that aren't applicable to the transcript
        n_start_pos = len(transcript) - insert_size + 1
        if (i_read1_first is None): i_read1_first = np.random.randint(0, n_start_pos) # uniform sampling of read location on the transcript
        i_read1_endp1 = i_read1_first + read_length
        i_read2_first = i_read1_first + insert_size - read_length
        i_read2_endp1 = i_read2_first + read_length
        read1 = transcript[i_read1_first:i_read1_endp1]
        read2 = transcript[i_read2_first:i_read2_endp1].reverse_complement() # paired-end reads
        # either bowtie or bam2hits mangles readid's with a space in them
        read1.id = "%s:first=%d:endp1=%d#0/1" % (readid, i_read1_first, i_read1_endp1) 
        read2.id = "%s:first=%d:endp1=%d#0/2" % (readid, i_read2_first, i_read2_endp1)
        read1.description = "" 
        read2.description = ""
        read1.letter_annotations["phred_quality"] = [40] * len(read1.seq)
        read2.letter_annotations["phred_quality"] = [40] * len(read2.seq)
        return(read1, read2)
        """
        
class SimulatedTranscriptome2(tp.NaivelyCachedComputation):
    def compute(self, transcriptPool, sensitivity, precision):
        """
        Result:
            transcriptPool TP transcripts in truthSet and knownSet
            transcriptPool FN transcripts in truthSet and not in knownSet
            transcriptPool FP transcripts in knownSet and not in truthSet
            
            fpTruthSet is a .fasta of truth set transcripts (for simulating reads)
            fpKnownSet is a .fasta of known set transcripts (for annotation-based expression estimation)
        """
        self.sensitivity = sensitivity
        self.precision = precision
        self.transcriptPool = transcriptPool

        # Transcripts from active genes => TP,FN according to sensitivity
        lExpressedSet = transcriptPool.lTruthSet
        random.shuffle(lExpressedSet)
        self.lTP = lExpressedSet[:int(sensitivity * len(lExpressedSet))]
        self.lFN = lExpressedSet[int(sensitivity * len(lExpressedSet)):]
        assert(len(lExpressedSet) == len(self.lTP) + len(self.lFN))

        # Calculate number of false positives /nFP/ needed to achieve precision /p/
        #     FP = FP_a + FP_s
        #     p = TP / (TP + FP) 
        #     FP = TP / p - TP
        self.nFP = int(len(self.lTP) / self.precision - len(self.lTP))

        # Shuffle FP transcripts and randomly select /fFP/ transcripts
        self.lFP = random.sample(transcriptPool.lUnexpressedSet, self.nFP)

        # Calculate inferred sensitivity and precision
        self.inferredSensitivity = len(self.lTP) / float(len(self.lTP) + len(self.lFN))
        self.inferredPrecision = len(self.lTP) / float(len(self.lTP) + len(self.lFP))

        # Compare {real, effective} {sensitivity, precision}
        tp.log("desired sensitivity: %.8f" % (self.sensitivity))
        tp.log("inferrd sensitivity: %.8f" % (self.inferredSensitivity))
        tp.log("desired precision: %.8f" % (self.precision))
        tp.log("inferrd precision: %.8f" % (self.inferredPrecision))

        # Generate a list of transcripts for TruthSet and KnownSet
        #self.lTruthSet = transcriptPool.lTruthSet
        self.lKnownSet = self.lTP + self.lFP

        # Generate {Truth,Known}Set .fa files
        #self.fpTruthSet = transcriptPool.fpTruthSet
        self.fpKnownSet = self.fp("knownSet.fa")
        #grepFastaForTranscriptID(fpInp=transcriptPool.fpTranscriptomeSubFa, fpOut=self.fpTruthSet, lTranscriptID=self.lTruthSet)
        grepFastaForTranscriptID(fpInp=transcriptPool.fpTranscriptomeSubFa, fpOut=self.fpKnownSet, lTranscriptID=self.lKnownSet)
        
        # Generate {Truth,Known}Set .gtf files
        #self.fpTruthSetGTF = transcriptPool.fpTruthSetGTF
        self.fpKnownSetGTF = self.fp("knownSet.gtf")
        #grepGtfForTranscriptIDs(fpInp=transcriptPool.fpTranscriptomeSubGtf, fpOut=self.fpTruthSetGTF, lTranscriptID=self.lTruthSet)
        grepGtfForTranscriptIDs(fpInp=transcriptPool.fpTranscriptomeSubGtf, fpOut=self.fpKnownSetGTF, lTranscriptID=self.lKnownSet)

        # Calculate expression levels for truth set and known set
        self.lKnownSetExpression = transcriptPool.getExpressionAsList(self.lKnownSet)
        tp.logList("KnownSet expression levels", self.lKnownSetExpression)

        """
        # Plot KnownSet generated expression levels
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(self.lKnownSetExpression, 100, normed=True)
        self.fpPlotKnownSetExpression = self.fp("plotKnownSetExpression.pdf")
        fig.savefig(self.fpPlotKnownSetExpression)
        tp.log("Wrote TruthSet expression level plot: %s" % (self.fpPlotKnownSetExpression,))
        """
    def getExpression(self, lTranscript):
        return self.transcriptPool.getExpressionAsList(lTranscript)

    def getTranscriptTPness(self, lTranscript):
        return [s in self.transcriptPool.lTruthSet for s in lTranscript]

    def getSensitivity(self):
        return(self.inferredSensitivity)

    def getPrecision(self):
        return(self.inferredPrecision)

class PlotMmseqPrediction2(tp.NaivelyCachedComputation):
    def compute(self,simulatedTranscriptome,mmseqExpr,transcriptPool):
        R = rpy2.robjects.r
        R.source(os.path.join(os.environ["RSSS_MMSEQ_DIR"], "mmseq.R"))
        # Read mmseq expression estimates
        R.assign("fn.expr.mmseq", mmseqExpr.fpExprMmseq)
        R("mmseq.out <- readmmseq(mmseq_files=fn.expr.mmseq, normalize=FALSE)")
        R("expr.mmseq <- exp(mmseq.out$log_mu)")
        R("expr.mmseq[ is.na(expr.mmseq) ] <- 0.0")

        # Plot TruthSet generated expression levels
        self.fpPlotMmseqExpression = self.fp("plotMmseqExpression.pdf")
        R.pdf(file=self.fpPlotMmseqExpression, width=6, height=6)
        R("hist(expr.mmseq, xlab=\"expression\", breaks=40)")
        R("dev.off()")
        tp.log("Wrote Mmseq expression level plot: %s" % (self.fpPlotMmseqExpression))

        # Calculate matching vector of true expression values
        R.assign("expr.truth", simulatedTranscriptome.getExpression([ s for s in R("row.names(mmseq.out$log_mu)") ]))
        R("expr.truth <- unlist(expr.truth)")
        # Calculate expression correlation
        R.assign("mTranscriptTP", simulatedTranscriptome.getTranscriptTPness([ s for s in R("row.names(mmseq.out$log_mu)") ]))
        R("mTranscriptTP <- unlist(mTranscriptTP)")
        self.correlation = R("cor(expr.truth, expr.mmseq)")[0]
        self.correlationTP = R("cor(expr.truth[mTranscriptTP], expr.mmseq[mTranscriptTP])")[0]
        self.logCorrelationTP = R("cor(log(expr.truth[mTranscriptTP & (expr.mmseq > 0)]), log(expr.mmseq[mTranscriptTP & (expr.mmseq > 0)]))")[0]
        #self.meanFPExprFrac = R("mean(expr.mmseq[!mTranscriptTP]) / mean(expr.mmseq[mTranscriptTP])")[0]
        self.meanFPExprFrac = R("sum(expr.mmseq[!mTranscriptTP]) / sum(expr.mmseq[mTranscriptTP])")[0]
        # Plot true vs estimated
        R.assign("fn.cor.plot", self.fp("plothCorrTruthMmseq.pdf"))
        R.assign("sens", simulatedTranscriptome.getSensitivity())
        R.assign("prec", simulatedTranscriptome.getPrecision())
        R("""
            pdf(file=fn.cor.plot, width=6, height=6)                                                           
            par(mfrow=c(1,1), pty='s')                                                                                  
            plot(expr.truth[!mTranscriptTP], expr.mmseq[!mTranscriptTP],
                main=sprintf("sens=%.2f, prec=%.2f, corAll=%.2f, corTP=%.2f, fr=%.8f",
                             sens, prec, cor(expr.truth, expr.mmseq), cor(expr.truth[mTranscriptTP], expr.mmseq[mTranscriptTP]), 
                             mean(expr.mmseq[!mTranscriptTP]) / mean(expr.mmseq[mTranscriptTP])),
                xlab="true expression",
                ylab="estimated expression",
                xlim=range(0, max(expr.truth)), ylim=range(0, max(expr.mmseq)),
                col="blue", pch="."
            )
            points(expr.truth[mTranscriptTP], expr.mmseq[mTranscriptTP], col='green', pch=".")
            #abline(a=0, b=1, col="gray")
            dev.off()
        """)

class TranscriptPool(tp.NaivelyCachedComputation):
    def compute(self):
        # Create python objects with transcript pool transcript sequences
        self.fpRawTranscripts = datafp("ensembl_human/release-66/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.66.cdna.all.fa.gz")
        fhRawTranscript = gzip.open(self.fpRawTranscripts)
        lRawTranscript = list(Bio.SeqIO.parse(fhRawTranscript, "fasta"))
        fhRawTranscript.close()
        self.iTranscripts = filter(lambda i: lRawTranscript[i].description.split()[2].split(":")[2] == '22', range(len((lRawTranscript))))
        self.lTranscripts = [lRawTranscript[i] for i in self.iTranscripts]
        # Create a .gtf with transcript pool annotations
        self.fpTranscriptsGtf = self.fp("transcriptPool.gtf")
        grepFile("^22", datafp("ensembl_human/release-66/gtf/homo_sapiens/Homo_sapiens.GRCh37.66.gtf.gz"), self.fpTranscriptsGtf)
        # Create a .fa with transcript pool genome sequences
        self.fpGenomeFa = datafp("ensembl_human/release-66/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.66.dna.chromosome.22.fa")
        """
        # look at what's in the description lines of the fasta
        d_status_count = {}
        d_chromosome_count = {}
        for transcript in l_transcript_pool:
            #print transcript.description
            transcript_status = transcript.description.split()[1].split(":")[1]
            d_status_count[transcript_status] = d_status_count.get(transcript_status, 0) + 1
            transcript_chromosome = transcript.description.split()[2].split(":")[2]
            d_chromosome_count[transcript_chromosome] = d_chromosome_count.get(transcript_chromosome, 0) + 1
    
        pprint.pprint(d_status_count)
        pprint.pprint(d_chromosome_count)
        """

class SimulatedTranscriptome(tp.NaivelyCachedComputation):
    def compute(self, transcriptPool, nTruthSet, completeness, purity, truthSetDensity=0.75):#, fpTranscriptPoolGTF="../../data/ensembl_human/Homo_sapiens.GRCh37.66.chr22.igvfix.gtf"):
        """
        Parameters:
            transcriptPool is the set of sequences to generate the truth set and known set
            nTruthSet is the number of transcripts in the truth set
            completeness is the fraction of truth set present in the current knowledge
            purity is the fraction of current knowledge present in the truth set

        Result:
            iTP transcripts in truthSet and knownSet
            iFN transcripts in truthSet and not in knownSet
            iFP transcripts in knownSet and not in truthSet
            
            fpTruthSet is the fasta file containing the truth set transcripts (for simulating reads)
            fpKnownSet is the fasta file containing the known set transcripts (for annotation-based expression estimation)
        """
        self.transcriptPool = transcriptPool

        # Sort a list of transcriptPool indices by the gene id
        def geneid(transcript): 
            return transcript.description.split()[3].split(":")[1]
        # Shuffle transcript such that transcripts with the same geneid are adjacent
        iTranscriptPool = range(len(transcriptPool.lTranscripts))
        iTranscriptPool.sort(cmp=lambda i, j: cmp(geneid(transcriptPool.lTranscripts[i]), geneid(transcriptPool.lTranscripts[j])))

        # Truncate and shuffle iTranscriptPool
        self.truthSetDensity = truthSetDensity
        self.truthSetCutoff = int(round(nTruthSet / self.truthSetDensity)) + 1
        iTranscriptPool = iTranscriptPool[:self.truthSetCutoff]#tp.log(iTranscriptPool)
        random.shuffle(iTranscriptPool)#tp.log(iTranscriptPool)

        # Calculate sizes for the various transcript groups
        self.nTruthSet = nTruthSet
        self.completeness = completeness
        self.purity = purity
        self.nTP = int(round(completeness * nTruthSet))
        self.nFN = nTruthSet - self.nTP
        self.nFP = int(round(self.nTP / purity - self.nTP))
        nTotal = self.nTP + self.nFP + self.nFN

        # Select TruthSet and KnownSet indices from the truncated transcript pool
        self.iTruthSet = iTranscriptPool[:self.nTP + self.nFN]
        self.iKnownSet = iTranscriptPool[:self.nTP] + iTranscriptPool[self.nTP+self.nFN:self.nTP+self.nFN+self.nFP]
 
        # Paranoid sanity check: no of TP / FN / FP should match KnownSet / TruthSet numbers
        tp.log(self.nTP + self.nFP + self.nFN)
        tp.log(len(self.iTruthSet) + len(self.iKnownSet) - self.nTP)
        assert self.nTP + self.nFP + self.nFN == len(self.iTruthSet) + len(self.iKnownSet) - self.nTP, "Check that purity large enough for given truthSetDensity?"
        
        # Write TruthSetFasta
        self.fpTruthSet = self.fp("TruthSet.fasta")
        fhTruthSet = open(self.fpTruthSet, "w")
        cRecords = 0
        for i in self.iTruthSet:
            cRecords += Bio.SeqIO.write(transcriptPool.lTranscripts[i], fhTruthSet, "fasta")
        assert(cRecords == len(self.iTruthSet))
        fhTruthSet.close()

        # Write KnownSetFasta
        self.fpKnownSet = self.fp("KnownSet.fasta")
        fhKnownSet = open(self.fpKnownSet, "w")
        cRecords = 0
        for i in self.iKnownSet:
            cRecords += Bio.SeqIO.write(transcriptPool.lTranscripts[i], fhKnownSet, "fasta")
        assert(cRecords == len(self.iKnownSet))
        fhKnownSet.close()

        # Write {Truth,Known}Set GTF files
        def transcriptid(transcript): 
            return transcript.description.split()[0]
        self.truthSetTranscripts = map(transcriptid, [transcriptPool.lTranscripts[i] for i in self.iTruthSet]) #tp.log(self.truthSetTranscripts)
        self.fpTruthSetGTF = self.fp("TruthSet.gtf")
        grepGtfForTranscriptIDs(fpInp=transcriptPool.fpTranscriptsGtf, fpOut=self.fpTruthSetGTF, lTranscriptID=self.truthSetTranscripts)
        self.knownSetTranscripts = map(transcriptid, [transcriptPool.lTranscripts[i] for i in self.iKnownSet]) #tp.log(self.knownSetTranscripts)
        self.fpKnownSetGTF = self.fp("KnownSet.gtf")
        grepGtfForTranscriptIDs(fpInp=transcriptPool.fpTranscriptsGtf, fpOut=self.fpKnownSetGTF, lTranscriptID=self.knownSetTranscripts)

        # Calculate expression levels (gamma distribution for truth set, zero for false positives in known set)
        self.lTruthSetExpression = np.random.gamma(shape = 1.2, scale = 1 / 0.001, size = len(self.truthSetTranscripts))
        self.lKnownSetExpression = np.concatenate((self.lTruthSetExpression[:self.nTP], np.zeros(shape=self.nFP)))
        tp.logList("TruthSet expression levels", self.lTruthSetExpression)
        ip.embed()
        #tp.logList("KnownSet expression levels", self.lKnownSetExpression)
        # Truncate true expression level distribution
        #threshold = 4000 # TWEAKME high value cutoff for the gamma distribution
        #while(any(self.lExpression > threshold)):
        #    self.lExpression[self.lExpression > threshold] = np.random.gamma(shape = 1.2, scale = 1 / 0.001, size = sum(self.lExpression > threshold))

        # Plot TruthSet generated expression levels
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.hist(self.lTruthSetExpression, 100, normed=True)
        #self.fpPlotTruthSetExpression = self.fp("plotTruthSetExpression.pdf")
        #fig.savefig(self.fpPlotTruthSetExpression)
        #tp.log("Wrote TruthSet expression level plot: %s" % (self.fpPlotTruthSetExpression,))

        # Plot KnownSet generated expression levels
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.hist(self.lKnownSetExpression, 100, normed=True)
        #self.fpPlotKnownSetExpression = self.fp("plotKnownSetExpression.pdf")
        #fig.savefig(self.fpPlotKnownSetExpression)
        #tp.log("Wrote TruthSet expression level plot: %s" % (self.fpPlotKnownSetExpression,))

    @staticmethod # no of FN transcripts to add to the truth set to get the specified completeness
    def n_fn_to_truthset(n_tp, p_completeness): return(int(round(n_tp / p_completeness - n_tp)))

    @staticmethod # no of FP transcripts to add to the known set to get the specified contamination
    def n_fp_to_knownset(n_tp, p_contamination): return(int(round(n_tp / (1 - p_contamination) - n_tp)))

class SimulatedTranscriptomeTests(unittest.TestCase):
    def test_n_fn_to_truthset_et_al(self):
        self.assertEqual(SimulatedTranscriptome.n_fn_to_truthset(100, 0.8), 25) # 25 <= 100 / (100 + 25) = 0.8
        self.assertEqual(SimulatedTranscriptome.n_fn_to_truthset(100, 0.5), 100) # 100 <= 100 / (100 + 100) = 0.5
        self.assertEqual(SimulatedTranscriptome.n_fp_to_knownset(100, 0.2), 25) # 25, as 25 / (100 + 25) = 0.2
        self.assertEqual(SimulatedTranscriptome.n_fp_to_knownset(100, 0.5), 100) # 100, as 100 / (100 + 100) = 0.5

class SimulatedReadsTests(unittest.TestCase):

    def test_simulate_transcript_reads_1(self):
        # simulate fragment with adjacent reads read pair
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGCTAGC"), id="test_simulate_transcript_reads")
        while(True):
            reads = SimulatedReads.simulateTranscriptRead(seq, read_length=4, insert_size_mean=8, insert_size_sd=10e-10)
            if not (reads is None):
                (read1, read2) = reads
                self.assertEqual(str(read1.seq), "ATGC")
                self.assertEqual(str(read2.seq), "GCTA")
                break

    def test_simulate_transcript_reads_2(self):
        # simulate fragment with one common base pair
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCTAGC"), id="test_simulate_transcript_reads")
        while(True):
            reads = SimulatedReads.simulateTranscriptRead(seq, read_length=4, insert_size_mean=7, insert_size_sd=10e-10)
            if not (reads is None):
                (read1, read2) = reads
                self.assertEqual(str(read1.seq), "ATCT")
                self.assertEqual(str(read2.seq), "GCTA")
                break

    def test_simulate_transcript_reads_3(self):
        # simulate fragment with a large dispersion in read pairs, look for specific cases
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCCCCAT"), id="test_simulate_transcript_reads")
        while(True):
            reads = SimulatedReads.simulateTranscriptRead(seq, read_length=2, insert_size_mean=8, insert_size_sd=4)
            if not (reads is None):
                (read1, read2) = reads
                if (str(read1.seq) == "AT") and (str(read2.seq) == "AT"):
                    self.assertTrue(True)
                    break

    def test_simulate_transcript_reads_4(self):
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("CCCATGCC"), id="test_simulate_transcript_reads")
        while(True):
            reads = SimulatedReads.simulateTranscriptRead(seq, read_length=2, insert_size_mean=4, insert_size_sd=4)
            if not (reads is None):
                (read1, read2) = reads #print str(read1.seq), str(read2.seq)
                if (str(read1.seq) == "AT") and (str(read2.seq) == "GC"):
                    self.assertTrue(True)
                    break

    """
    # modify for paired reads
    def test_simulate_transcriptome_reads_1(self):
        # simulate transcriptome with two transcripts of different length
        fh_testfile1 = tempfile.NamedTemporaryFile()
        fn_testfile1 = fh_testfile1.name
        fh_testfile2 = tempfile.NamedTemporaryFile()
        fn_testfile2 = fh_testfile2.name
        SimulatedReads.simulateTranscriptomeReads(l_transcript = [
                Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("AAAA"), id="test1_transcript1"),
                Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("GGGGGGGG"), id="test1_transcript2")
            ],
            l_expression = [6, 6],
            n_frags = 3,
            read_length = 2,
            insert_size_mean = 4,
            insert_size_sd = 10e-10,
            fn_output1 = fn_testfile1,
            fn_output2 = fn_testfile2,
            addPoissonNoise=False
        )
        l_transcript_r1 = list(Bio.SeqIO.parse(fn_testfile1, "fastq"))
        n_read_from_transcript1 = sum(map(lambda seq: 1 if (str(seq.seq).find("A") > -1) else 0, l_transcript_r1))
        n_read_from_transcript2 = sum(map(lambda seq: 1 if (str(seq.seq).find("G") > -1) else 0, l_transcript_r1))

        print 2 * n_read_from_transcript1, n_read_from_transcript2
        
        self.assertTrue(2 * n_read_from_transcript1 == n_read_from_transcript2)
        l_transcript_r2 = list(Bio.SeqIO.parse(fn_testfile2, "fastq"))
        n_read_from_transcript1 = sum(map(lambda seq: 1 if (str(seq.seq).find("T") > -1) else 0, l_transcript_r2))
        n_read_from_transcript2 = sum(map(lambda seq: 1 if (str(seq.seq).find("C") > -1) else 0, l_transcript_r2))
        self.assertTrue(2 * n_read_from_transcript1 == n_read_from_transcript2)
        fh_testfile1.close()
        fh_testfile2.close()
    """

    """    
    def test_simulate_transcriptome_reads_2(self):
        # simulate transcriptome with two transcripts of different expression
        fh_testfile = tempfile.NamedTemporaryFile()
        fn_testfile = fh_testfile.name
        SimulatedReads.simulateTranscriptomeReads(l_transcript = [
                Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("CCCC"), id="test2_transcript1"),
                Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("TTTT"), id="test2_transcript2")
            ],
            l_expression = [7, 14],
            n_frags = 3,
            read_length = 2,
            insert_size_mean = 2,
            insert_size_sd = 1,
            fn_output = fn_testfile,
            addPoissonNoise=False
        )
        l_transcript_res = list(Bio.SeqIO.parse(fn_testfile, "fastq"))
        n_read_from_transcript1 = sum(map(lambda seq: 1 if (str(seq.seq).find("C") > -1) else 0, l_transcript_res))
        n_read_from_transcript2 = sum(map(lambda seq: 1 if (str(seq.seq).find("T") > -1) else 0, l_transcript_res))
        self.assertTrue(2 * n_read_from_transcript1 == n_read_from_transcript2)
        fh_testfile.close()
    """

class BowtieIndex(tp.NaivelyCachedComputation):
    def compute(self, simulatedTranscriptome):
        self.fpBowtieIndex = self.fp("bowtie_index")
        cmd = "bowtie-build --offrate 3 -f '%(transcriptome)s' '%(bowtieIndex)s'"
        tp.run(cmd % {
                      'transcriptome': simulatedTranscriptome.fpKnownSet,
                      'bowtieIndex': self.fpBowtieIndex
                     })

class BowtieAlign(tp.NaivelyCachedComputation):
    def compute(self, simulatedReads, bowtieIndex):
        self.fpAlignment = self.fp("BowtieAlign.namesorted")
        cmd = "bowtie --chunkmbs 256 -p %(nproc)d -a --best --strata --fullref --sam -m 100 -I 0 -X 500 '%(bowtie_index)s' -1 '%(reads1)s' -2 '%(reads2)s' | samtools view -F 0xC -bS - | samtools sort -n - '%(namesorted)s'"
        #cmd = "bowtie -a --best --strata --fullref --sam '%(bowtie_index)s' '%(reads1)s' | samtools view -bS - | samtools sort -n - '%(namesorted)s'"
        tp.run(cmd % {"nproc": nproc(),
                      "bowtie_index": bowtieIndex.fpBowtieIndex,
                      "reads1": simulatedReads.fpReads1,
                      "reads2": simulatedReads.fpReads2,
                      "namesorted": self.fpAlignment,
                      })
        self.fpAlignment = self.fpAlignment + ".bam"                          

class TophatBowtieIndex(tp.NaivelyCachedComputation):
    def compute(self, transcriptPool):
        self.transcriptPool = transcriptPool
        self.fpGenomeFa = transcriptPool.fpGenomeFa
        self.fpBowtieIndex = self.fp("bowtieIndex")
        cmd = "bowtie-build --offrate 3 -f '%(referenceSeq)s' '%(bowtieIndex)s'"
        tp.run(cmd % {'referenceSeq': self.fpGenomeFa,
                      'bowtieIndex': self.fpBowtieIndex
                      })

class TophatAlign(tp.NaivelyCachedComputation):
    def compute(self, simulatedReads, tophatBowtieIndex, organism="human"):
        self.simulatedReads = simulatedReads
        self.tophatBowtieIndex = tophatBowtieIndex
        self.fpAlignment = self.fp("tophat_out")
        cmd = "tophat --bowtie1 -p %(nproc)d -o '%(tophatOut)s' '%(bowtieIndex)s' '%(reads1)s' '%(reads2)s'"
        #cmd = "tophat --no-novel-juncs -G '%(annotGTF)s' -o '%(tophatOut)s' '%(bowtieIndex)s' '%(reads1)s' '%(reads2)s'"

        # parameter setting suggested by RGASP, 50 for human as default
        if organism is "warm":
            cmd = "tophat -i 30 --min-coverage-intron 30 --min-segment-intron 30 -p %(nproc)d -o '%(tophatOut)s' '%(bowtieIndex)s' '%(reads1)s' '%(reads2)s'"

        if organism is "fly":
            cmd = "tophat -i 40 --min-coverage-intron 40 --min-segment-intron 40 -p %(nproc)d -o '%(tophatOut)s' '%(bowtieIndex)s' '%(reads1)s' '%(reads2)s'"
        
        tp.run(cmd % {"nproc": nproc(),
                      "tophatOut": self.fpAlignment,
                      "bowtieIndex": tophatBowtieIndex.fpBowtieIndex,
                      "reads1": simulatedReads.fpReads1,
                      "reads2": simulatedReads.fpReads2#,
                      #"annotGTF": simulatedReads.simulatedTranscriptome.fpTruthSetGTF
                      })

class Cufflinks(tp.NaivelyCachedComputation):
    def compute(self, tophatAlign, simulatedTranscriptome, RABTassembly = False, organism="human"):
        self.fpClOut = self.fp("cl_out")
        #cmd = "cufflinks -o '%(clOut)s' '%(thAcceptedHits)s'" # estimate insert size from data (estimates seem to be reliable)
        #cmd = "cufflinks -G '%(annotGTF)s' -o '%(clOut)s' '%(thAcceptedHits)s' -m 250 -s 30" #integrate directly from simulatedReads
        #cmd = "cufflinks -o '%(clOut)s' '%(thAcceptedHits)s' -m 250 -s 30" #integrate directly from simulatedReads

        cmd = "cufflinks -p %(nproc)d -o '%(clOut)s' %(RABTassembly)s '%(thAcceptedHits)s' -m %(insertSizeMean)s -s %(insertSizeSd)s"

        # parameter setting suggested by RGASP, 50 for human as default
        if organism is "warm":
            cmd = "cufflinks --min-intron-length 30 -p %(nproc)d -o '%(clOut)s' %(RABTassembly)s '%(thAcceptedHits)s' -m %(insertSizeMean)s -s %(insertSizeSd)s"

        if organism is "fly":
           cmd = "cufflinks --min-intron-length 40 -p %(nproc)d -o '%(clOut)s' %(RABTassembly)s '%(thAcceptedHits)s' -m %(insertSizeMean)s -s %(insertSizeSd)s"

        tp.run(cmd % {"nproc": nproc(),
                      "clOut": self.fpClOut,
                      "thAcceptedHits": tophatAlign.fpAlignment + "/accepted_hits.bam",
                      "RABTassembly": "-g '%s'" % (simulatedTranscriptome.fpKnownSetGTF,) if RABTassembly else "",
                      "insertSizeMean": tophatAlign.simulatedReads.insertSizeMean,
                      "insertSizeSd": tophatAlign.simulatedReads.insertSizeSd
                      })

class CufflinksTranscriptome1(tp.NaivelyCachedComputation):
    def compute(self, cufflinks, simulatedTranscriptome):
        self.transcriptPool = simulatedTranscriptome.transcriptPool

        # Load TruthSet and Cufflinks annotations
        self.fpKnownSetGTF = cufflinks.fpClOut + "/transcripts.gtf"
        dTrCufflinks = aggregateExonBoundaries(readGTF(self.fpKnownSetGTF))
        dTrTruthSet = aggregateExonBoundaries(readGTF(simulatedTranscriptome.transcriptPool.fpTruthSetGTF))

        # Generate .fasta for transcripts
        self.fpGenomeFa = simulatedTranscriptome.transcriptPool.fpGenomeFa
        self.fpKnownSet = self.fp("KnownSet.fasta")
        cmd = "extract_transcripts '%(referenceSeq)s' '%(transcriptsGTF)s' > '%(transcriptsFasta)s'"
        tp.run(cmd % {'referenceSeq': self.fpGenomeFa,
                      'transcriptsGTF': self.fpKnownSetGTF,
                      'transcriptsFasta': self.fpKnownSet
                      })

        # Match cufflinks transcripts to truth set transcripts, allowing for mismatches in the beginning/end
        lFullMatch = [] # list of 2-tuples for full matches (exon boundaries and start regions)
        lExonMatch = [] # list of 2-tuples for exon boundary matches that aren't full matches
        for clTid in dTrCufflinks:
            for trTid in dTrTruthSet:
                """
                if clTid == trTid and (dTrCufflinks[clTid]["coords"] != dTrTruthSet[trTid]["coords"]):
                    tp.log("Mismatch!")
                    tp.log(",".join(map(str, dTrCufflinks[clTid]["coords"])))
                    tp.log(",".join(map(str, dTrTruthSet[trTid]["coords"])))
                """
                if (dTrCufflinks[clTid]["coords"] == dTrTruthSet[trTid]["coords"]):
                    lFullMatch.append((clTid, trTid))
                elif (dTrCufflinks[clTid]["coords"][1:-1] == dTrTruthSet[trTid]["coords"][1:-1]):
                    lExonMatch.append((clTid, trTid))

        # Classify Cufflinks transcripts as TP or FP, generate expected expression values for transcripts
        self.nTP = 0
        self.nFP = 0
        self.knownSetTranscripts = dTrCufflinks.keys()
        self.knownSetTranscriptsMatch = { i: i for i in dTrCufflinks.keys() } # matches CUFF* to ENST* if possible
        self.lKnownSetExpression = np.zeros(len(dTrCufflinks.keys()))
        for (i, clTid) in enumerate(self.knownSetTranscripts):
            fullMatchFw = [k for (j, k) in lFullMatch if j == clTid]
            fullMatchBw = [j for (j, k) in lFullMatch if k == fullMatchFw[0]] if len(fullMatchFw) == 1 else []
            exonMatchFw = [k for (j, k) in lExonMatch if j == clTid]
            exonMatchBw = [j for (j, k) in lExonMatch if k == exonMatchFw[0]] if len(exonMatchFw) == 1 else []

            anyMatch = False
            if len(fullMatchFw) == 1 and len(fullMatchBw) == 1:
                self.nTP += 1
                self.knownSetTranscriptsMatch[clTid] = fullMatchFw[0]
                self.lKnownSetExpression[i] = simulatedTranscriptome.transcriptPool.getExpression([fullMatchFw[0]])[0]
                #self.lKnownSetExpression[i] = simulatedTranscriptome.lTruthSetExpression[simulatedTranscriptome.truthSetTranscripts.index(fullMatchFw[0])]
                anyMatch = True
                #tp.log("Full match: %s %s" % (fullMatchFw[0],fullMatchBw[0]))

            if len(exonMatchFw) == 1 and len(exonMatchBw) == 1:
                anyMatch = True
                if len(fullMatchBw) == 0:
                    self.nTP += 1
                    self.knownSetTranscriptsMatch[clTid] = exonMatchFw[0]
                    self.lKnownSetExpression[i] = simulatedTranscriptome.transcriptPool.getExpression([exonMatchFw[0]])[0]
                    #self.lKnownSetExpression[i] = simulatedTranscriptome.lTruthSetExpression[simulatedTranscriptome.truthSetTranscripts.index(exonMatchFw[0])]
                    #tp.log("Exon match: %s %s" % (exonMatchFw[0],exonMatchBw[0]))
                else:
                    self.nFP += 1
                    self.lKnownSetExpression[i] = 0.0
                    #tp.log("Exon match => FP, since full match exists!")

            if not anyMatch:
                self.nFP += 1
                self.lKnownSetExpression[i] = 0.0
                #tp.log("No unique matches found => FP")

        # False negatives: truth set transcripts that can't be matched with anything in the cufflinks output
        self.nFN = simulatedTranscriptome.transcriptPool.nActiveExTranscript - self.nTP
        tp.log("nTP: %f\tnFP: %f\tnFN: %f" % (self.nTP, self.nFP, self.nFN))

        # compatibility with plotMmseq
        self.lKnownSet = self.knownSetTranscripts

    def getExpression(self, lTranscript):
        return self.transcriptPool.getExpression( [ self.knownSetTranscriptsMatch[trid] for trid in lTranscript ] )

    def writeIGVSessionFile():
        self.fnIGVSessionXML = "CufflinksTruthSetGenome.xml"
        fhIGVSessionXML = open(self.fp(self.fnIGVSessionXML), "w")
        fhIGVSessionXML.write(
'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="%(fpSequenceFa)s" locus="chr22:45251945-45460798" version="4">
    <Resources>
        <Resource path="%(fpCufflinksGTF)s"/>
        <Resource path="%(fpTruthSetGTF)s"/>
    </Resources>
    <Panel height="778" name="FeaturePanel" width="1663">
        <Track altColor="0,0,178" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>
        <Track altColor="0,0,178" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" height="45" id="%(fpTruthSetGTF)s" name="Truth set transcripts" renderer="GENE_TRACK" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
        <Track altColor="0,0,178" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" height="45" id="%(fpCufflinksGTF)s" name="Cufflinks prediction" renderer="GENE_TRACK" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
    </Panel>
<PanelLayout dividerFractions="0.0076045627376425855"/>
</Session>''' % {
               'fpSequenceFa': os.path.abspath("../../data/ensembl_human/chr22.fa"), 
               'fpTruthSetGTF': os.path.abspath(simulatedTranscriptome.fpTruthSetGTF),
               'fpCufflinksGTF': os.path.abspath(fpClTranscripts)
        })
        fhIGVSessionXML.close()


class CufflinksTranscriptome2(tp.NaivelyCachedComputation):
    def compute(self, cufflinks, simulatedTranscriptome, maxLengthDiffFrac = 0.2, minCommonSeqFrac=0.8):
        self.transcriptPool = simulatedTranscriptome.transcriptPool

        # Load TruthSet and Cufflinks annotations
        self.fpKnownSetGTF = cufflinks.fpClOut + "/transcripts.gtf"
        dKnownSet = aggregateExonBoundaries(readGTF(self.fpKnownSetGTF))
        dTruthSet = aggregateExonBoundaries(readGTF(simulatedTranscriptome.transcriptPool.fpTruthSetGTF))
        self.dKnownSet = dKnownSet
        self.dTruthSet = dTruthSet

        # Generate .fasta for transcripts
        self.fpGenomeFa = simulatedTranscriptome.transcriptPool.fpGenomeFa
        self.fpKnownSet = self.fp("KnownSet.fasta")
        cmd = "extract_transcripts '%(referenceSeq)s' '%(transcriptsGTF)s' > '%(transcriptsFasta)s'"
        tp.run(cmd % {'referenceSeq': self.fpGenomeFa,
                      'transcriptsGTF': self.fpKnownSetGTF,
                      'transcriptsFasta': self.fpKnownSet
                      })

        # Calculate the length of a transcript from an array of exon coordinates
        def trLength(coords): return(sum([ abs(coords[i] - coords[i+1]) + 1 for i in range(0, len(coords), 2) ]))

        # Generate match scores by iterating over all (known, true) transcript pairs
        lMatchCandidates = []
        nFullMatches = 0
        for (knCount, knId) in enumerate(dKnownSet.keys()):
            tp.log("matching transcripts: %d / %d" % (knCount+1, len(dKnownSet.keys())))
            for (trCount, trId) in enumerate(dTruthSet.keys()):
                knCoord = dKnownSet[knId]['coords']
                trCoord = dTruthSet[trId]['coords']
                # Check strands
                if (dKnownSet[knId]['strand'] != dTruthSet[trId]['strand']): continue
                # Filter out "merged" mis-assemblies (=transcript lengths are different)
                lengthDiffFrac = abs(trLength(knCoord) - trLength(trCoord)) / float(trLength(trCoord))
                if lengthDiffFrac > maxLengthDiffFrac: continue
                # Check that longest common substring is long enough
                lCommon = min(trLength(knCoord), trLength(trCoord))
                lTrue = trLength(trCoord)
                if (float(lCommon) / float(lTrue)) < minCommonSeqFrac: continue
                # Single-exon transcripts: check that they overlap
                if len(knCoord)==2 and len(trCoord)==2 and ((max(knCoord) < min(trCoord)) or (max(trCoord) < min(knCoord))): continue
                # Check inner exon boundary match
                if (knCoord[1:-1] == trCoord[1:-1]):
                    # Calculate "mismatch" score
                    lMismatch = abs(min(knCoord) - min(trCoord))
                    rMismatch = abs(max(knCoord) - max(trCoord))
                    score = lMismatch**2 + rMismatch**2
                    lMatchCandidates.append((knId, trId, score))

        # Filter candidate pairs according to increasing score (to get a one-to-one mapping)
        lMatchCandidates.sort(key=operator.itemgetter(2))
        self.knownSetTranscriptsMatch = {}
        # While there are still match candidates left
        while(len(lMatchCandidates) > 0):
            iMatch = lMatchCandidates.pop(0)
            # Add lMatchCandidate to lMatches
            self.knownSetTranscriptsMatch[iMatch[0]] = iMatch[1]
            # Remove match candidates containing transcripts from the current accepted match
            lMatchCandidates = filter(lambda x: x[0] != iMatch[0] and x[1] != iMatch[1], lMatchCandidates)

    def getExpression(self, lTranscript):
        return self.transcriptPool.getExpressionAsList( [ self.knownSetTranscriptsMatch.get(trid, "CUFFL_FP_TRANSCRIPT") for trid in lTranscript ] )

    def getTranscriptTPness(self, lTranscript):
        return [s in self.knownSetTranscriptsMatch.keys() for s in lTranscript]

    def getSensitivity(self):
        self.nTP = len(self.knownSetTranscriptsMatch.keys())
        self.nFP = len(self.dKnownSet) - self.nTP
        self.nFN = len(self.dTruthSet) - self.nTP
        return(float(self.nTP) / (self.nTP + self.nFN))

    def getPrecision(self):
        self.nTP = len(self.knownSetTranscriptsMatch.keys())
        self.nFP = len(self.dKnownSet) - self.nTP
        self.nFN = len(self.dTruthSet) - self.nTP
        return(float(self.nTP) / (self.nTP + self.nFP))

    def writeIGVSessionFile():
        self.fnIGVSessionXML = "CufflinksTruthSetGenome.xml"
        fhIGVSessionXML = open(self.fp(self.fnIGVSessionXML), "w")
        fhIGVSessionXML.write(
'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="%(fpSequenceFa)s" locus="chr22:45251945-45460798" version="4">
    <Resources>
        <Resource path="%(fpCufflinksGTF)s"/>
        <Resource path="%(fpTruthSetGTF)s"/>
    </Resources>
    <Panel height="778" name="FeaturePanel" width="1663">
        <Track altColor="0,0,178" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>
        <Track altColor="0,0,178" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" height="45" id="%(fpTruthSetGTF)s" name="Truth set transcripts" renderer="GENE_TRACK" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
        <Track altColor="0,0,178" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" height="45" id="%(fpCufflinksGTF)s" name="Cufflinks prediction" renderer="GENE_TRACK" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
    </Panel>
<PanelLayout dividerFractions="0.0076045627376425855"/>
</Session>''' % {
               'fpSequenceFa': os.path.abspath("../../data/ensembl_human/chr22.fa"), 
               'fpTruthSetGTF': os.path.abspath(simulatedTranscriptome.fpTruthSetGTF),
               'fpCufflinksGTF': os.path.abspath(fpClTranscripts)
        })
        fhIGVSessionXML.close()

class Oases(tp.NaivelyCachedComputation):
    def compute(self, simulatedReads):
        # from: shuffleSequences_fasta.py in velvet
        self.fpShuffledFa = self.fp("shuffledReads.fq")

        fh1 = open(simulatedReads.fpReads1)
        fh2 = open(simulatedReads.fpReads2)
        def interleave(iter1, iter2):                                                                                          
            while True:
                yield iter1.next()
                yield iter2.next()
        outfile = open(self.fpShuffledFa, 'w')
        records = interleave(Bio.SeqIO.parse(fh1,"fastq"), Bio.SeqIO.parse(fh2, "fastq"))
        count = Bio.SeqIO.write(records, outfile, "fasta")
        outfile.close()

        cmd = "%(oases_pipeline)s -m 21 -M 23 -o %(fpOut)s -p '-ins_length %(insert_size)s' -d %(fas)s"
        tp.run(cmd % {"oases_pipeline": "OMP_THREAD_LIMIT=%d; oases_pipeline.py" % (nproc(),),
                      "fas": self.fpShuffledFa,
                      "insert_size": simulatedReads.insertSizeMean,
                      "fas": self.fpShuffledFa,
                      "fpOut": self.fp("")
                      })

class OasesTranscriptome(tp.NaivelyCachedComputation):
    def compute(self, oases, transcriptPool):
        self.transcriptPool = transcriptPool
        self.oases = oases
        self.fpKnownSet = self.oases.fp("Merged/transcripts.fa")
        fhKnownSet = open(self.fpKnownSet, "r")
        self.lKnownSet = list(Bio.SeqIO.parse(fhKnownSet,"fasta"))
        fhKnownSet.close()
        self.fpTruthSet = self.transcriptPool.fpTruthSet
        fhTruthSet = open(self.fpTruthSet, "r")
        self.lTruthSet = list(Bio.SeqIO.parse(fhTruthSet,"fasta"))
        fhTruthSet.close()
        # Generate match scores by iterating over all (known, true) transcript pairs
        maxLengthDiffFrac = 0.2
        lMatchCandidates = []
        dTruthSetExons = aggregateExonBoundaries(readGTF(transcriptPool.fpTruthSetGTF))
        nFullMatches = 0
        for (knCount, knRec) in enumerate(self.lKnownSet):
            tp.log("matching transcripts: %d / %d" % (knCount+1, len(self.lKnownSet)))
            for (trCount, trRec) in enumerate(self.lTruthSet):
                #tp.log("matching transcripts: %d / %d, %d / %d" % (knCount, len(self.lKnownSet), trCount, len(self.lTruthSet)))
                if (str(knRec.seq) == str(trRec.seq)): nFullMatches = nFullMatches + 1
                nExons = len(dTruthSetExons[trRec.id]["coords"]) / 2
                if nExons < 3:
                    trueLeExonLen = None
                    trueRiExonLen = None
                else:
                    trueLeExonLen = abs(dTruthSetExons[trRec.id]["coords"][1] - dTruthSetExons[trRec.id]["coords"][0])
                    trueRiExonLen = abs(dTruthSetExons[trRec.id]["coords"][-1] - dTruthSetExons[trRec.id]["coords"][-2])
                score = OasesTranscriptome.transcriptMatch(str(knRec.seq), str(trRec.seq), trueLeExonLen=trueLeExonLen, trueRiExonLen=trueRiExonLen, maxLengthDiffFrac=maxLengthDiffFrac)
                if score < float('inf'):
                    lMatchCandidates.append((knRec.id, trRec.id, score))
        # Filter candidate pairs according to increasing score (to get a one-to-one mapping)
        lMatchCandidates.sort(key=operator.itemgetter(2))
        self.knownSetTranscriptsMatch = {}
        # While there are still match candidates left
        while(len(lMatchCandidates) > 0):
            iMatch = lMatchCandidates.pop(0)
            # Add lMatchCandidate to lMatches
            self.knownSetTranscriptsMatch[iMatch[0]] = iMatch[1]
            # Remove match candidates containing transcripts from the current accepted match
            lMatchCandidates = filter(lambda x: x[0] != iMatch[0] and x[1] != iMatch[1], lMatchCandidates)

    def getExpression(self, lTranscript):
        return self.transcriptPool.getExpressionAsList( [ self.knownSetTranscriptsMatch.get(trid, "OASES_FP_TRANSCRIPT") for trid in lTranscript ] )

    def getTranscriptTPness(self, lTranscript):
        return [s in self.knownSetTranscriptsMatch.keys() for s in lTranscript]

    def getSensitivity(self):
        self.nTP = len(self.knownSetTranscriptsMatch.keys())
        self.nFP = len(self.lKnownSet) - self.nTP
        self.nFN = len(self.lTruthSet) - self.nTP
        return(float(self.nTP) / (self.nTP + self.nFN))

    def getPrecision(self):
        self.nTP = len(self.knownSetTranscriptsMatch.keys())
        self.nFP = len(self.lKnownSet) - self.nTP
        self.nFN = len(self.lTruthSet) - self.nTP
        return(float(self.nTP) / (self.nTP + self.nFP))

    @staticmethod
    def longestCommonSubstrings(s, t):
        """
        Returns starting indices of the longest common substrings, both in s and t
        http://en.wikipedia.org/wiki/Longest_common_substring_problem
        """
        z = 0
        ret = []
        lp = [0] * len(t)
        ln = [0] * len(t)
        for i in range(len(s)):
            for j in range(len(t)):
                if s[i] == t[j]:
                    if (i == 0) or (j == 0): 
                        ln[j] = 1
                    else:
                        ln[j] = lp[j - 1] + 1
                    if ln[j] > z:
                        z = ln[j]
                        ret = []
                    if ln[j] == z:
                        ret = ret + [(i + 1 - z, j + 1 - z)]
            else:
                ln[j] = 0
            lp = ln
            ln = [0] * len(t)
        return (ret, z)

    @staticmethod
    def transcriptMatch(candSeq, trueSeq, trueLeExonLen=None, trueRiExonLen=None, maxLengthDiffFrac=0.2, minCommonSeqFrac=0.8):
        """
            1) candidate transcript length within 20% of true transcript length
            2) longest common substring of the candidate sequence and the true sequence covers 80% of the true transcript
            3) if true sequence has inner exons:
            4) longest common substring of the candidate sequence and the true sequence covers all inner exons
            5) if several LCS fit the criteria above, prefer the "match" closest to the centre
        """
        # Speedup: discard obvious mismatches quickly
        if not(trueLeExonLen is None) and not(trueRiExonLen is None):
            if not(trueSeq[trueLeExonLen:-trueRiExonLen] in candSeq): return float('inf')
        else:
            maxCutoff = int(len(trueSeq) * (1 - minCommonSeqFrac) / 2)
            if not(trueSeq[maxCutoff:-maxCutoff] in candSeq): return float('inf')

        # Filter out "merged" mis-assemblies (=transcript lengths are different)
        lengthDiffFrac = abs(len(candSeq) - len(trueSeq)) / float(len(trueSeq))
        if lengthDiffFrac > maxLengthDiffFrac: return float('inf')

        # Find longest common substring, check out that longest common substring is long enough
        (lStart, lCommon) = OasesTranscriptome.longestCommonSubstrings(candSeq, trueSeq)
        if (float(lCommon) / len(trueSeq)) < minCommonSeqFrac: return float("inf")
        # Iterate over all longest common substrings
        minScore = float('inf')
        for start in lStart:
            # Check that all inner exon bases are covered by the LCS
            if (not(trueLeExonLen is None) and not(start[1] <= trueLeExonLen)) or \
               (not(trueRiExonLen is None) and not(len(trueSeq) - trueRiExonLen) <= (start[1] + lCommon)):
                continue
            # Calculate "mismatch" score
            lMismatch = max(start[0], start[1])
            rMismatch = max(len(candSeq) - start[0] - lCommon, len(trueSeq) - start[1] - lCommon)
            thisScore = lMismatch**2 + rMismatch**2
            if (thisScore < minScore):
                minScore = thisScore
        return(minScore)

class OasisTranscriptomeTests(unittest.TestCase):
    def test_transcriptMatch1(self):
        """
        Test longestCommonSubstrings calculation with maxLengthDiffFrac and minCommonSeqFrac "disabled".
        """
        self.assertEqual(OasesTranscriptome.transcriptMatch("TGC",   "ATGCGCATGC", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 37)
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGC",  "ATGCGCATGC", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 36)
        self.assertEqual(OasesTranscriptome.transcriptMatch("GCGCA", "ATGCGCATGC", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 13)
        self.assertEqual(OasesTranscriptome.transcriptMatch("GC",    "ATGCGCATGC", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 32)
        self.assertEqual(OasesTranscriptome.transcriptMatch("GGTGC", "ATGCGCATGC", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 40)

    def test_transcriptMatch2(self):
        """
        Test trueLeExonLen and trueRiExonLen calculation
        True:  AG|TGCGGA|TA
               AG TGCGGA TA    => 0
                C TGCGGA G     => 8
                  TGCGGA       => 8
                  TGCGGA AAA   => 13
             CCCC TGCGGA       => 20
                G TGCGGA T     => 2
                G AGCGGA T     => inf
                G TGCGGA T     => 2
                G TGCGGC T     => inf
        1"""
        self.assertEqual(OasesTranscriptome.transcriptMatch(  "AGTGCGGATA",  "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 0)
        self.assertEqual(OasesTranscriptome.transcriptMatch(   "CTGCGGAG",   "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 8)
        self.assertEqual(OasesTranscriptome.transcriptMatch(    "TGCGGA",    "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 8)
        self.assertEqual(OasesTranscriptome.transcriptMatch(    "TGCGGAAAA", "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 13)
        self.assertEqual(OasesTranscriptome.transcriptMatch("CCCCTGCGGA",    "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 20)
        self.assertEqual(OasesTranscriptome.transcriptMatch(   "GTGCGGAT",   "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 2)
        self.assertEqual(OasesTranscriptome.transcriptMatch(   "GAGCGGAT",   "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), float('inf'))
        self.assertEqual(OasesTranscriptome.transcriptMatch(   "GTGCGGAT",   "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), 2)
        self.assertEqual(OasesTranscriptome.transcriptMatch(   "GTGCGGCT",   "AGTGCGGATA", trueLeExonLen=2, trueRiExonLen=2, maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0), float('inf'))

    def test_transcriptMatch3(self):
        """
        Test maxLengthDiffFrac
        True:  ATGCATGCAT, maxLengthDiffFrac=0.2
               ATGCATG       => float('inf')
               ATGCATGC      => 4
               ATGCATGCA     => 1
               ATGCATGCAT    => 0
               ATGCATGCATG   => 1
               ATGCATGCATGC  => 4
               ATGCATGCATGCA => float('inf')
        """
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATG",       "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), float('inf'))
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATGC",      "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), 4)
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATGCA",     "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), 1)
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATGCAT",    "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), 0)
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATGCATG",   "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), 1)
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATGCATGC",  "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), 4)
        self.assertEqual(OasesTranscriptome.transcriptMatch("ATGCATGCATGCA", "ATGCATGCAT", maxLengthDiffFrac=0.2, minCommonSeqFrac=0), float('inf'))

    def test_transcriptMatch4(self):
        """
        Test minCommonSeqFrac:
        True:  AAAAAAAAAA, maxLengthDiffFrac=0.2
               AAAAAAATTT      => float('inf')
               AAAAAAAATT      => 4
               AAAAAAAAAT      => 4
               AAAAAAAAAA      => 4
               AAAAAAAAAAA     => 4
               AAAAAAAAAAAA    => 4
        """
        self.assertEqual(OasesTranscriptome.transcriptMatch("AAAAAAATTT",   "AAAAAAAAAA", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0.8), float('inf'))
        self.assertEqual(OasesTranscriptome.transcriptMatch("AAAAAAAATT",   "AAAAAAAAAA", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0.8), 4)
        self.assertEqual(OasesTranscriptome.transcriptMatch("AAAAAAAAAT",   "AAAAAAAAAA", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0.8), 1)
        self.assertEqual(OasesTranscriptome.transcriptMatch("AAAAAAAAAA",   "AAAAAAAAAA", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0.8), 0)
        self.assertEqual(OasesTranscriptome.transcriptMatch("AAAAAAAAAAA",  "AAAAAAAAAA", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0.8), 1)
        self.assertEqual(OasesTranscriptome.transcriptMatch("AAAAAAAAAAAA", "AAAAAAAAAA", maxLengthDiffFrac=float('inf'), minCommonSeqFrac=0.8), 2)

def datafp(s):
    # This should be changed to e.g. reading the path from an environmental variable...
    return os.path.join(os.environ["RSSS_DATA_DIR"], s)

def jobid():
    if 'PBS_ARRAYID' in os.environ:
        return(int(os.environ['PBS_ARRAYID']))
    else:
        return None

def nproc():
    """ 
    For login nodes, use 4 cores for testing.
    For cluster runs, use n-1 cores.
    """
    if os.uname()[1][0:5] == "login": 
        return(4)
    else:
        return multiprocessing.cpu_count() - 1

def grepFile(regex, fpInp, fpOut):
    cmd = "%(grep)s \"%(regex)s\" '%(fpInp)s' > '%(fpOut)s'"
    tp.run(cmd % {"grep": ("zgrep" if fpInp[-3:] == ".gz" else "grep"),
                  "regex": regex,
                  "fpInp": fpInp,
                  "fpOut": fpOut
                  })

def readGTF(fp):
    "Parse an track-line-free Ensembl/Cufflinks GTF file into a python list of dictionaries."
    lGtfRecord = []
    dExonBoundaries = {}
    for l in open(fp, "r"):
        tokens = l.split("\t")
        gtfRecord = {
            "seqname": tokens[0],
            "source": tokens[1],
            "feature": tokens[2],
            "start": int(tokens[3]),
            "end": int(tokens[4]),
            "score": tokens[5],#float(tokens[5]),
            "strand": tokens[6],
            "frame": tokens[7],
            "attribute": OrderedDict([((token.lstrip(" ").partition(" "))[0], token.lstrip(" ").partition(" ")[2].strip('\"')) 
                                      for token in tokens[8].rstrip(';\n').split(";")])
        }
        lGtfRecord.append(gtfRecord)
    return(lGtfRecord)

def readGTFTest():
    coverage = 1000
    nTruthSet = 100
    completeness = 1.0
    purity = 1.0
    transcriptPool = TranscriptPool()
    simulatedTranscriptome = SimulatedTranscriptome(nTruthSet=nTruthSet, completeness=completeness, purity=purity, transcriptPool=transcriptPool) #tp.log(truthSet.iTranscriptPoolInTruthSet)
    simulatedReads = SimulatedReads(nTotalFrags=nTruthSet*coverage,simulatedTranscriptome=simulatedTranscriptome,insertSizeMean=250, insertSizeSd=30)
    tophatBowtieIndex = TophatBowtieIndex(transcriptPool=transcriptPool)
    tophatAlign = TophatAlign(simulatedReads=simulatedReads, tophatBowtieIndex=tophatBowtieIndex)
    cufflinks = Cufflinks(tophatAlign=tophatAlign)
    pprint(readGTF(cufflinks.fpClOut + "/transcripts.gtf"))

def aggregateExonBoundaries(lGtfRecord):
    "Create a dictionary of {transcript_id: {coords, strand}} from a GTF table where coords are exon boundaries."
    dTranscript = {}
    for gtfRecord in lGtfRecord:
        if (gtfRecord["feature"] != "exon"): continue
        transcriptId = gtfRecord["attribute"]["transcript_id"]
        if not(transcriptId in dTranscript.keys()):
            dTranscript[transcriptId] = {"coords": [gtfRecord["start"], gtfRecord["end"]]}
            dTranscript[transcriptId]["strand"] = gtfRecord["strand"]
        else:
            dTranscript[transcriptId]["coords"].append(gtfRecord["start"])
            dTranscript[transcriptId]["coords"].append(gtfRecord["end"])
    for transcriptVal in dTranscript.values():
        transcriptVal["coords"].sort()
    return(dTranscript)

def aggregateExonBoundariesTest():
    coverage = 1000
    nTruthSet = 100
    completeness = 1.0
    purity = 1.0
    transcriptPool = TranscriptPool()
    simulatedTranscriptome = SimulatedTranscriptome(nTruthSet=nTruthSet, completeness=completeness, purity=purity, transcriptPool=transcriptPool) #tp.log(truthSet.iTranscriptPoolInTruthSet)
    simulatedReads = SimulatedReads(nTotalFrags=nTruthSet*coverage,simulatedTranscriptome=simulatedTranscriptome,insertSizeMean=250, insertSizeSd=30)
    tophatBowtieIndex = TophatBowtieIndex(transcriptPool=transcriptPool)
    tophatAlign = TophatAlign(simulatedReads=simulatedReads, tophatBowtieIndex=tophatBowtieIndex)
    cufflinks = Cufflinks(tophatAlign=tophatAlign)
    #pprint(readGTF(cufflinks.fpClOut + "/transcripts.gtf"))
    pprint(aggregateExonBoundaries(readGTF(cufflinks.fpClOut + "/transcripts.gtf")))

class Hitsfile(tp.NaivelyCachedComputation):
    def compute(self, simulatedTranscriptome, bowtieAlign, fpBam2hits="bam2hits", enstFormat=True):
        self.fpHitsfile = self.fp("hitsfile")
        cmd = 'OMP_NUM_THREADS=%(nproc)d; %(bam2hits)s %(match_regex)s \'%(fa_transcripts)s\' \'%(bam_reads)s\' > \'%(hitsfile)s\''
        # Ensembl-specific transcripts
        #cmd = '%(bam2hits)s -m "(E\S+).*gene:(E\S+).*" 1 2 \'%(fa_transcripts)s\' \'%(bam_reads)s\' > \'%(hitsfile)s\''
        # manual insert size filtering
        #cmd = '%(bam2hits)s -m "(E\S+).*gene:(E\S+).*" 1 2 -i 244 41 \'%(fa_transcripts)s\' \'%(bam_reads)s\' > \'%(hitsfile)s\''
        ## TODO para match_regex?
        tp.run(cmd % {"nproc": nproc(),
                      "bam2hits": fpBam2hits,
                      "fa_transcripts": simulatedTranscriptome.fpKnownSet,
                      "bam_reads": bowtieAlign.fpAlignment,
                      "hitsfile": self.fpHitsfile,
                      "match_regex": '-m "(\S+).*gene:(\S+).*" 1 2' if enstFormat else '-m "(.*)" 1 1',
                      })             

class Mmseq(tp.NaivelyCachedComputation):
    def compute(self, hitsfile, fpMmseq="mmseq"):
        cmd = "OMP_NUM_THREADS=%(nproc)d; %(mmseq)s -alpha 1.2 -beta 0.001 '%(hitsfile)s' '%(output_base)s'"
        #cmd = "%(mmseq)s -alpha 1.2 -beta 0.001 -gibbs_iter 1048576 '%(hitsfile)s' '%(output_base)s'"
        self.fpOutput = self.fp("expr")
        tp.run(cmd % {"nproc": nproc(),
                      "mmseq": fpMmseq,
                      "hitsfile": hitsfile.fpHitsfile, 
                      "output_base": self.fpOutput
                      })
        self.fpExprMmseq = self.fpOutput + ".mmseq"

if __name__ == '__main__':
    tp.log("Running unit tests")
    unittest.main(exit=False)
    random.seed(1)
    numpy.random.seed(1)
    #plotCoverageVsClReconstruction()
