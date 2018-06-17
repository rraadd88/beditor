#!usr/bin/python

import logging
import subprocess
import logging
import shutil
import re
import sys
import os
import logging, operator

from collections import defaultdict
from os.path import join, isfile, basename, dirname, isdir, abspath, exists
from os import makedirs

from Bio import SeqIO

def parseFasta(fileObj): #FIXME replace with biopython
    " parse a fasta file, where each seq is on a single line, return dict id -> seq "
    seqs = {}
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if seqId!=None:
                seqs[seqId]  = "".join(parts)
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs[seqId]  = "".join(parts)
    return seqs

def pamIsCpf1(pam):
    " if you change this, also change bin/filterFaToBed! "
    return (pam in ["TTN", "TTTN", "TYCV", "TATV"])
        
def pamIsSaCas9(pam):
    " only used for notes and efficiency scores, unlike its Cpf1 cousin function "
    return (pam in ["NNGRRT", "NNNRRT"])

def pamIsXCas9(pam):
    " "
    return (pam in ["NGK", "NGN"])
    
def pamIsSpCas9(pam):
    " only used for notes and efficiency scores, unlike its Cpf1 cousin function "
    return (pam in ["NGG", "NGA", "NGCG"])

def setupPamInfo(pam,addGenePlasmids,scoreNames):
    " modify a few globals based on the current pam "
    PAMLEN = len(pam)
    if pamIsCpf1(pam):
        logging.debug("switching on Cpf1 mode, guide length is 23bp")
        GUIDELEN = 23
        cpf1Mode = True
        scoreNames = cpf1ScoreNames
    elif pam=="NNGRRT" or pam=="NNNRRT":
        logging.debug("switching on S. aureus mode, guide length is 21bp")
        addGenePlasmids = addGenePlasmidsAureus
        GUIDELEN = 21
        cpf1Mode = False
        scoreNames = saCas9ScoreNames
    else:
        GUIDELEN = 20
        cpf1Mode = False
    return GUIDELEN,cpf1Mode,addGenePlasmids,PAMLEN,scoreNames

def runCmd(cmd):
    from beditor.lib.global_vars import dirs2ps 
    cmd = cmd.replace("$BIN", dirs2ps['binDir'])
    cmd = cmd.replace("$PYTHON", dirs2ps['pyp'])
    cmd = cmd.replace("$SCRIPT", dirs2ps['scriptDir'])
#     print(cmd)
    err=subprocess.call(cmd,shell=True)
    if err!=0:
        print('bash command error: {}\n{}\n'.format(err,cmd))
        sys.exit(1)
           

#---
revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M',
    "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n", "V" : "B", "v":"b", 
    "B" : "V", "b": "v", "W" : "W", "w" : "w"}

def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

#--
def flankSeqIter(seq, startDict, pamLen, doFilterNs,GUIDELEN,pamPlusLen,cpf1Mode=False):
    """ given a seq and dictionary of pamPos -> strand and the length of the pamSite
    yield tuples of (name, pamStart, guideStart, strand, flankSeq, pamSeq)

    flankSeq is the guide sequence (=flanking the PAM).

    if doFilterNs is set, will not return any sequences that contain an N character
    pamPlusSeq are the 5bp after the PAM. If not enough space, pamPlusSeq is None
    """

    startList = sorted(startDict.keys())
    for pamStart in startList:
        strand = startDict[pamStart]

        pamPlusSeq = None
        if cpf1Mode: # Cpf1: get the sequence to the right of the PAM
            if strand=="+":
                guideStart = pamStart+pamLen
                flankSeq = seq[guideStart:guideStart+GUIDELEN]
                pamSeq = seq[pamStart:pamStart+pamLen]
                if pamStart-pamPlusLen >= 0:
                    pamPlusSeq = seq[pamStart-pamPlusLen:pamStart]
            else: # strand is minus
                guideStart = pamStart-GUIDELEN
                flankSeq = revComp(seq[guideStart:pamStart])
                pamSeq = revComp(seq[pamStart:pamStart+pamLen])
                if pamStart+pamLen+pamPlusLen < len(seq):
                    pamPlusSeq = revComp(seq[pamStart+pamLen:pamStart+pamLen+pamPlusLen])
        else: # common case: get the sequence on the left side of the PAM
            if strand=="+":
                guideStart = pamStart-GUIDELEN
                flankSeq = seq[guideStart:pamStart]
                pamSeq = seq[pamStart:pamStart+pamLen]
                if pamStart+pamLen+pamPlusLen < len(seq):
                    pamPlusSeq = seq[pamStart+pamLen:pamStart+pamLen+pamPlusLen]
            else: # strand is minus
                guideStart = pamStart+pamLen
                flankSeq = revComp(seq[guideStart:guideStart+GUIDELEN])
                pamSeq = revComp(seq[pamStart:pamStart+pamLen])
                if pamStart-pamPlusLen >= 0:
                    pamPlusSeq = revComp(seq[pamStart-pamPlusLen:pamStart])

        if "N" in flankSeq and doFilterNs:
            continue

        yield "s%d%s" % (pamStart, strand), pamStart, guideStart, strand, flankSeq, pamSeq, pamPlusSeq
        
def writePamFlank(seq, startDict, pam, faFname,GUIDELEN,pamPlusLen):
    " write pam flanking sequences to fasta file, optionally with versions where each nucl is removed "
    #print "writing pams to %s<br>" % faFname
    faFh = open(faFname, "w")
    for pamId, pamStart, guideStart, strand, flankSeq, pamSeq, pamPlusSeq in flankSeqIter(seq, startDict, len(pam), True,GUIDELEN,pamPlusLen):
        faFh.write(">%s\n%s\n" % (pamId, flankSeq))
    faFh.close()
    
def annotateBedWithPos(inBed, outBed):
    """
    given an input bed4 and an output bed filename, add an additional column 5 to the bed file
    that is a descriptive text of the chromosome pos (e.g. chr1:1.23 Mbp).
    """
    ofh = open(outBed, "w")
    for line in open(inBed):
        chrom, start = line.split("\t")[:2]
        start = int(start)

        if start>1000000:
            startStr = "%.2f Mbp" % (float(start)/1000000)
        else:
            startStr = "%.2f Kbp" % (float(start)/1000)
        desc = "%s %s" % (chrom, startStr)

        ofh.write(line.rstrip("\n"))
        ofh.write("\t")
        ofh.write(desc)
        ofh.write("\n")
    ofh.close()
# def findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname,maxMMs,genomep,
#         GUIDELEN,
#         MFAC,
#         MAXOCC,
#         offtargetPams,
#         ALTPAMMINSCORE,
#         DEBUG):    
#     " align faFname to genome and create matchedBedFname "
#     print('queue')
#     print(queue)
#     print('batchId')
#     print(batchId)
#     print('batchBase')
#     print(batchBase)
#     print('faFname')
#     print(faFname)
#     print('genome')
#     print(genome)
#     print('pam')
#     print(pam)
#     print('bedFname')
#     print(bedFname)


#     matchesBedFname = batchBase+".matches.bed"
#     saFname = batchBase+".sa"
#     pamLen = len(pam)
#     genomeDir = dirname(genomep) # make var local, see below

#     open(matchesBedFname, "w") # truncate to 0 size

#     # increase MAXOCC if there is only a single query, but only in CGI mode

#     maxDiff = maxMMs
#     queue.startStep(batchId, "bwa", "Alignment of potential guides, mismatches <= %d" % maxDiff)
#     convertMsg = "Converting alignments"
#     seqLen = GUIDELEN

#     bwaM = MFAC*MAXOCC # -m is queue size in bwa
#     print(maxDiff)
#     print(seqLen)
#     print('ssssssssssssssss')
#     cmd = "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomep)s %(faFname)s > %(saFname)s" % locals()
#     runCmd(cmd)

#     queue.startStep(batchId, "saiToBed", convertMsg)
#     maxOcc = MAXOCC # create local var from global
#     # EXTRACTION OF POSITIONS + CONVERSION + SORT/CLIP
#     # the sorting should improve the twoBitToFa runtime
    
#     cmd = "$BIN/bwa samse -n %(maxOcc)d %(genomep)s %(saFname)s %(faFname)s | $SCRIPT/xa2multi.pl | $SCRIPT/samToBed %(pam)s %(seqLen)d | sort -k1,1 -k2,2n | $BIN/bedClip stdin %(genomep)s.sizes stdout >> %(matchesBedFname)s " % locals()
#     runCmd(cmd)

#     # arguments: guideSeq, mainPat, altPats, altScore, passX1Score
#     filtMatchesBedFname = batchBase+".filtMatches.bed"
#     queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
#     altPats = ",".join(offtargetPams.get(pam, ["na"]))
#     bedFnameTmp = bedFname+".tmp"
#     altPamMinScore = str(ALTPAMMINSCORE)
#     cmd = "bedtools getfasta -name -s -fi {} -bed {} -fo {}.fasta".format(genomep,matchesBedFname,matchesBedFname)
#     runCmd(cmd)
#     cmd = "python3 $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d %(matchesBedFname)s.fasta > %(filtMatchesBedFname)s" % locals()
#     runCmd(cmd)

#     segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

#     # if we have gene model segments, annotate them, otherwise just use the chrom position
#     if isfile(segFname):
#         queue.startStep(batchId, "genes", "Annotating matches with genes")
#         cmd = "cat %(filtMatchesBedFname)s | $BIN/overlapSelect %(segFname)s stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 > %(bedFnameTmp)s " % locals()
#         runCmd(cmd)
#     else:
#         queue.startStep(batchId, "chromPos", "Annotating matches with chromosome position")
#         annotateBedWithPos(filtMatchesBedFname, bedFnameTmp)

#     # make sure the final bed file is never in a half-written state, 
#     # as it is our signal that the job is complete
#     shutil.move(bedFnameTmp, bedFname)
#     queue.startStep(batchId, "done", "Job completed")

#     # remove the temporary files
#     tempFnames = [saFname, matchesBedFname, filtMatchesBedFname]
#     for tfn in tempFnames:
#         if isfile(tfn) and not DEBUG:
#             os.remove(tfn)
#     return bedFname

# def processSubmission(faFname, genome, pam, bedFname, batchBase, batchId, queue,maxMMs,genomep,params_findofftargetbwa,doEffScoring=False, cpf1Mode=False):
#     """ search fasta file against genome, filter for pam matches and write to bedFName 
#     optionally write status updates to work queue. Remove faFname.
#     """
#     if doEffScoring and not cpf1Mode:
#         queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
#         createBatchEffScoreTable(batchId)
        
# #     if useBowtie:
# #         findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pam, bedFname)
# #     else:
#     findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname,maxMMs,genomep=genomep,**params_findofftargetbwa)
#     print(faFname)
#     return bedFname

def getOfftargets(seq, org, pam, batchId, startDict, queue,TEMPDIR,GUIDELEN,pamPlusLen,maxMMs,genomep,
                MFAC,
                MAXOCC,
                offtargetPams,
                ALTPAMMINSCORE,
                DEBUG):
    """ write guides to fasta and run bwa or use cached results.
    Return name of the BED file with the matches.
    Write progress status updates to queue object.
    """
#     print(params_findofftargetbwa)
    batchBase = join(TEMPDIR, batchId)
    otBedFname = batchBase+".bed"
    bedFname=otBedFname
    print(otBedFname)
    faFname = batchBase+".fa"

#     if not isfile(otBedFname):
#     # write potential PAM sites to file 
#         faFname = batchBase+".fa"
#         writePamFlank(seq, startDict, pam, faFname,GUIDELEN,pamPlusLen)
#         processSubmission(faFname, org, pam, otBedFname, batchBase, batchId, queue,maxMMs,genomep,params_findofftargetbwa)
#         if doEffScoring and not cpf1Mode:
#             queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
#             createBatchEffScoreTable(batchId)

#         findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname,maxMMs,genomep=genomep,**params_findofftargetbwa)
#     print(faFname)
# def findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pam, bedFname,maxMMs,genomep,
#         GUIDELEN,
#         MFAC,
#         MAXOCC,
#         offtargetPams,
#         ALTPAMMINSCORE,
#         DEBUG):    
#     " align faFname to genome and create matchedBedFname "

    matchesBedFname = batchBase+".matches.bed"
    saFname = batchBase+".sa"
    pamLen = len(pam)
    genomeDir = dirname(genomep) # make var local, see below

    open(matchesBedFname, "w") # truncate to 0 size

    # increase MAXOCC if there is only a single query, but only in CGI mode

    maxDiff = maxMMs
    queue.startStep(batchId, "bwa", "Alignment of potential guides, mismatches <= %d" % maxDiff)
    convertMsg = "Converting alignments"
    seqLen = GUIDELEN

    bwaM = MFAC*MAXOCC # -m is queue size in bwa
#     print(maxDiff)
#     print(seqLen)
    print('ssssssssssssssss')
    cmd = "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomep)s %(faFname)s > %(saFname)s" % locals()
    runCmd(cmd)

    queue.startStep(batchId, "saiToBed", convertMsg)
    maxOcc = MAXOCC # create local var from global
    # EXTRACTION OF POSITIONS + CONVERSION + SORT/CLIP
    # the sorting should improve the twoBitToFa runtime
    
    cmd = "$BIN/bwa samse -n %(maxOcc)d %(genomep)s %(saFname)s %(faFname)s | $SCRIPT/xa2multi.pl | $SCRIPT/samToBed %(pam)s %(seqLen)d | sort -k1,1 -k2,2n | $BIN/bedClip stdin %(genomep)s.sizes stdout >> %(matchesBedFname)s " % locals()
    runCmd(cmd)

    # arguments: guideSeq, mainPat, altPats, altScore, passX1Score
#     filtMatchesBedFname = batchBase+".filtMatches.bed"
#     queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
#     altPats = ",".join(offtargetPams.get(pam, ["na"]))
#     bedFnameTmp = bedFname+".tmp"
#     altPamMinScore = str(ALTPAMMINSCORE)
    cmd = "bedtools getfasta -name -s -fi {} -bed {} -fo {}.fasta".format(genomep,matchesBedFname,matchesBedFname)
    runCmd(cmd)
#     cmd = "python3 $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d %(matchesBedFname)s.fasta > %(filtMatchesBedFname)s" % locals()
#     runCmd(cmd)

#     genome=org #FIXME
#     segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

#     # if we have gene model segments, annotate them, otherwise just use the chrom position
#     if isfile(segFname):
#         queue.startStep(batchId, "genes", "Annotating matches with genes")
#         cmd = "cat %(filtMatchesBedFname)s | $BIN/overlapSelect %(segFname)s stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 > %(bedFnameTmp)s " % locals()
#         runCmd(cmd)
#     else:
#         queue.startStep(batchId, "chromPos", "Annotating matches with chromosome position")
#         annotateBedWithPos(filtMatchesBedFname, bedFnameTmp)

#     # make sure the final bed file is never in a half-written state, 
#     # as it is our signal that the job is complete
#     shutil.move(bedFnameTmp, bedFname)
#     queue.startStep(batchId, "done", "Job completed")

#     # remove the temporary files
#     tempFnames = [saFname, matchesBedFname, filtMatchesBedFname]
#     for tfn in tempFnames:
#         if isfile(tfn) and not DEBUG:
#             os.remove(tfn)
# #     return bedFname

    return matchesBedFname
#--

def intToExtPamId(pamId,pamIdRe):
    " convert the internal pam Id like s20+ to the external one, like 21Forw "
    pamPos, strand, rest = pamIdRe.match(pamId).groups()
    if strand=="+":
        strDesc = 'forw'
    else:
        strDesc = 'rev'
    guideDesc = str(int(pamPos)+1)+strDesc
    return guideDesc

def iterGuideRows(guideData, guideHeaders,scoreNames,
                     cpf1Mode,scoreDescs,pamIdRe,addHeaders=False, seqId=None, minSpec=None, minFusi=None):
    "yield rows from guide data. Need to know if for Cpf1 or not "
    headers, tableScoreNames = makeGuideHeaders(guideHeaders,scoreNames,
                     cpf1Mode,scoreDescs)

#     if satMutOpt:
#         headers.append("Oligonucleotide")
#         headers.append("AdapterHandle+PrimerFw")
#         headers.append("AdapterHandle+PrimerRev")
#         oligoPrefix, oligoSuffix, primerFwPrefix, primerRevPrefix, batchId, genome, position, ampLen, tm = satMutOpt

#         batchBase = join(TEMPDIR, batchId)
#         otBedFname = batchBase+".bed"
#         otMatches = parseOfftargets(otBedFname)

#         guideData.sort(key=operator.itemgetter(3)) # sort by position, makes more sense here

# #     if seqId != None:
# #         headers.insert(0, "#seqId")
#     else:
#         headers[0] = "#"+headers[0]

    if addHeaders:
        yield headers

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount = guideRow
        if minSpec and guideScore < minSpec:
            continue
        if minFusi and effScores["fusi"] < minFusi:
            continue

        otCount = 0
        if otData!=None:
            otCount = len(otData)

        guideDesc = intToExtPamId(pamId,pamIdRe)

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)
        row = [guideDesc, fullSeq, guideScore, otCount, ontargetDesc]

#         for scoreName in tableScoreNames:
#             row.append(effScores.get(scoreName, "NotEnoughFlankSeq"))

        row = [str(x) for x in row]
#         if seqId != None:
#             row.insert(0, seqId)
        yield row
    
def makeGuideHeaders(guideHeaders,scoreNames,
                     cpf1Mode,scoreDescs):
    " return list of the headers of the guide output file "
    headers = list(tuple(guideHeaders)) # make a copy of the list

    logging.debug("active scoreNames: %s" % scoreNames)
    tableScoreNames = list(tuple(scoreNames))
    if not cpf1Mode:
        tableScoreNames.append("oof")

    for scoreName in tableScoreNames:
        headers.append(scoreDescs[scoreName][0]+"-Score")
    return headers, tableScoreNames

def parseOfftargets(bedFname):
    """ parse a bed file with annotataed off target matches from overlapSelect,
    has two name fields, one with the pam position/strand and one with the
    overlapped segment 
    
    return as dict pamId -> editDist -> (chrom, start, end, seq, strand, segType, segName, x1Score)
    segType is "ex" "int" or "ig" (=intergenic)
    if intergenic, geneNameStr is two genes, split by |
    """
    # example input:
    # chrIV 9864393 9864410 s41-|-|5|ACTTGACTG|0    chrIV   9864303 9864408 ex:K07F5.16
    # chrIV   9864393 9864410 s41-|-|5|ACTGTAGCTAGCT|9999    chrIV   9864408 9864470 in:K07F5.16
    logging.info("reading offtargets from %s" % bedFname)

    # first sort into dict (pamId,chrom,start,end,editDist,strand) 
    # -> (segType, segName) 
    pamData = {}
    for line in open(bedFname):
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, name, segment = fields
        # hg38: ignore alternate chromosomes otherwise the 
        # regions on the main chroms look as if they could not be 
        # targeted at all with Cas9
        if chrom.endswith("_alt"):
            continue
        nameFields = name.split("|")
        pamId, strand, editDist, seq = nameFields[:4]

        if len(nameFields)>4:
            x1Count = int(nameFields[4])
        else:
            x1Count = 0
        editDist = int(editDist)
        # some gene models include colons
        if ":" in segment:
            segType, segName = segment.split(":", maxsplit=1)
        else:
            segType, segName = "", segment
        start, end = int(start), int(end)
        otKey = (pamId, chrom, start, end, editDist, seq, strand, x1Count)

        # if a offtarget overlaps an intron/exon or ig/exon boundary it will
        # appear twice; in this case, we only keep the exon offtarget
        if otKey in pamData and segType!="ex":
            continue
        pamData[otKey] = (segType, segName)

    # index by pamId and edit distance
    indexedOts = defaultdict(dict)
    for otKey, otVal in pamData.items():
        pamId, chrom, start, end, editDist, seq, strand, x1Score = otKey
        segType, segName = otVal
        otTuple = (chrom, start, end, seq, strand, segType, segName, x1Score)
        indexedOts[pamId].setdefault(editDist, []).append( otTuple )

    return indexedOts
#--

def concatGuideAndPam(guideSeq, pamSeq,cpf1Mode=False, pamPlusSeq=""):
    " return guide+pam or pam+guide, depending on cpf1Mode "
    if cpf1Mode:
        return pamPlusSeq+pamSeq+guideSeq
    else:
        return guideSeq+pamSeq+pamPlusSeq
def makePosList(org, countDict, guideSeq, pam, inputPos, maxMMs,cpf1Mode,
               GUIDELEN,MFAC,MAXOCC,offtargetPams,ALTPAMMINSCORE,DEBUG):
    """     
    for a given guide sequence, return a list of tuples that
    describes the offtargets sorted by score and a string to describe the offtargets in the
    format x/y/z/w of mismatch counts
    inputPos has format "chrom:start-end:strand". All 0MM matches in this range
    are ignored from scoring ("ontargets")
    Also return the same description for just the last 12 bp and the score
    of the guide sequence (calculated using all offtargets).
    """
    inChrom, inStart, inEnd, inStrand = parsePos(inputPos)
#     count = 0
    otCounts = []
    posList = []
    mitOtScores = []
    cfdScores = []
    last12MmCounts = []
    ontargetDesc = ""
    subOptMatchCount = 0

    # for each edit distance, get the off targets and iterate over them
    foundOneOntarget = False
    for editDist in range(0, maxMMs+1):
        #print countDict,"<p>"
        matches = countDict.get(editDist, [])

        #print otCounts,"<p>"
        last12MmOtCount = 0

        # create html and score for every offtarget
        otCount = 0
        for chrom, start, end, otSeq, strand, segType, geneNameStr, x1Count in matches:
            # skip on-targets
            if segType!="":
                segTypeDesc = segTypeConv[segType]
                geneDesc = segTypeDesc+":"+geneNameStr
                geneDesc = geneDesc.replace("|", "-")
            else:
                geneDesc = geneNameStr

            # is this the on-target?
            # if we got a genome position, use it. Otherwise use a random off-target with 0MMs
            # as the on-target ("auto-ontarget")
            if editDist==0 and \
                ((chrom==inChrom and start >= inStart and end <= inEnd and x1Count < MAXOCC) \
                or (inChrom=='' and foundOneOntarget==False and x1Count < MAXOCC)):
                foundOneOntarget = True
                ontargetDesc = geneDesc
                continue

            otCount += 1
            guideNoPam = guideSeq[:len(guideSeq)-len(pam)]
            otSeqNoPam = otSeq[:len(otSeq)-len(pam)]

            if len(otSeqNoPam)==19:
                otSeqNoPam = "A"+otSeqNoPam # should not change the score a lot, weight0 is very low
                guideNoPam = "A"+guideNoPam

#             if pamIsCpf1(pam):
                # Cpf1 has no scores yet
            mitScore=0.0
            cfdScore=0.0
#             else:
#                 # MIT score must not include the PAM
# #                 mitScore = calcHitScore(guideNoPam, otSeqNoPam)
#                 # this is a heuristic based on the guideSeq data where alternative
#                 # PAMs represent only ~10% of all cleaveage events.
#                 # We divide the MIT score by 5 to make sure that these off-targets 
#                 # are not ranked among the top but still appear in the list somewhat
#                 if pam=="NGG" and otSeq[-2:]!="GG":
#                     mitScore = mitScore * 0.2

                # CFD score must include the PAM
#                 cfdScore = calcCfdScore(guideSeq, otSeq)

            mitOtScores.append(mitScore)
            if cfdScore != None:
                cfdScores.append(cfdScore)

            posStr = "%s:%d-%s:%s" % (chrom, int(start)+1,end, strand)
            if (chrom==inChrom):
                dist = abs(start-inStart)
            else:
                dist = None

            alnHtml, hasLast12Mm = makeAlnStr(org, guideSeq, otSeq, pam, mitScore, cfdScore, posStr, dist,cpf1Mode)
            if not hasLast12Mm:
                last12MmOtCount+=1
            posList.append( (otSeq, mitScore, cfdScore, editDist, posStr, geneDesc, alnHtml) )
            # taking the maximum is probably not necessary, 
            # there should be only one offtarget for X1-exceeding matches
            subOptMatchCount = max(int(x1Count), subOptMatchCount)

        last12MmCounts.append(str(last12MmOtCount))
        # create a list of number of offtargets for this edit dist
        otCounts.append( str(otCount) )

    # calculate the guide scores
    guideScore = 0
    guideCfdScore = 0

#     if pamIsCpf1(pam):
#         guideScore = -1
#         guideCfdScore = -1
#     else:
#         if subOptMatchCount > MAXOCC:
#             guideScore = 0
#             guideCfdScore = 0
#         else:
#             guideScore = calcMitGuideScore(sum(mitOtScores))
#             guideCfdScore = calcMitGuideScore(sum(cfdScores))

    # obtain the off-target info: coordinates, descriptions and off-target counts
    if subOptMatchCount > MAXOCC:
        posList = []
        ontargetDesc = ""
        last12DescStr = ""
        otDescStr = ""
    else:
        otDescStr = "&thinsp;-&thinsp;".join(otCounts)
        last12DescStr = "&thinsp;-&thinsp;".join(last12MmCounts)

    if pamIsCpf1(pam):
        # sort by edit dist if using Cfp1
        posList.sort(key=operator.itemgetter(3))
    else:
        # sort by CFD score if we have it
        posList.sort(reverse=True, key=operator.itemgetter(2))

    return posList, otDescStr, guideScore, guideCfdScore, last12DescStr, \
        ontargetDesc, subOptMatchCount

def parsePos(text):
    """ parse a string of format chr:start-end:strand and return a 4-tuple
    Strand defaults to + and end defaults to start+23
    """
    print(text)
    if text!=None and len(text)!=0 and text!="?":
        fields = text.split(":")
        if len(fields)==2:
            chrom, posRange = fields
            strand = "+"
        else:
            chrom, posRange, strand = fields
        posRange = posRange.replace(",","")
        if "-" in posRange:
            start, end = posRange.split("-")
            start, end = int(start), int(end)
        else:
            # if the end position is not specified (as by default done by UCSC outlinks), use start+23
            start = int(posRange)
            end = start+23
    else:
        chrom, start, end, strand = "", 0, 0, "+"
    return chrom, start, end, strand

def mergeGuideInfo(seq, startDict, pamPat, otMatches, inputPos, 
                   GUIDELEN,pamPlusLen,maxMMs,cpf1Mode,params_findofftargetbwa,effScores=None, sortBy=None, org=None,sort=False):
    """
    merges guide information from the sequence, the efficiency scores and the off-targets.
    creates rows with too many fields. Probably needs refactoring.


    for each pam in startDict, retrieve the guide sequence next to it and score it
    sortBy can be "effScore", "mhScore", "oofScore" or "pos"
    """
#     allEnzymes = readEnzymes()

    guideData = []
    guideScores = {}
    hasNotFound = False
    pamIdToSeq = {}

    pamSeqs = list(flankSeqIter(seq.upper(), startDict, len(pamPat), True,GUIDELEN,pamPlusLen))

    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            pamMatches = otMatches[pamId]
            guideSeqFull = concatGuideAndPam(guideSeq, pamSeq)
            mutEnzymes = []#matchRestrEnz(allEnzymes, guideSeq, pamSeq, pamPlusSeq)
            posList, otDesc, guideScore, guideCfdScore, last12Desc, ontargetDesc, \
               subOptMatchCount = \
                   makePosList(org, pamMatches, guideSeqFull, pamPat, inputPos,maxMMs,cpf1Mode,**params_findofftargetbwa)

        # no off-targets found?
        else:
            posList, otDesc, guideScore = None, "Not found", None
            guideCfdScore = None
            last12Desc = ""
            hasNotFound = True
            mutEnzymes = []
            ontargetDesc = ""
            subOptMatchCount = False
            seq34Mer = None

#         guideRow = [guideScore, guideCfdScore, effScores.get(pamId, {}), pamStart, guideStart, strand, pamId, guideSeq, pamSeq, posList, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount]
        guideRow = [guideScore, guideCfdScore, effScores, pamStart, guideStart, strand, pamId, guideSeq, pamSeq, posList, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount]
        guideData.append( guideRow )
        guideScores[pamId] = guideScore
        pamIdToSeq[pamId] = guideSeq

    if sortBy == "pos":
        sortFunc = (lambda row: row[3])
        reverse = False
    elif sortBy is not None and sortBy!="spec":
        sortFunc = (lambda row: row[2].get(sortBy, 0))
        reverse = True
    else:
        sortFunc = operator.itemgetter(0)
        reverse = True
    if sort:
        guideData.sort(reverse=reverse, key=sortFunc)
    
    return guideData, guideScores, hasNotFound, pamIdToSeq
#--
def makeAlnStr(org, seq1, seq2, pam, mitScore, cfdScore, posStr, chromDist,cpf1Mode):
    " given two strings of equal length, return a html-formatted string that highlights the differences "
    lines = [ [], [], [] ]
    last12MmCount = 0

    if cpf1Mode:
        lines[0].append("<i>"+seq1[:len(pam)]+"</i> ")
        lines[1].append("<i>"+seq2[:len(pam)]+"</i> ")
        lines[2].append("".join([" "]*(len(pam)+1)))

    if cpf1Mode:
        guideStart = len(pam)
        guideEnd = len(seq1)
    else:
        guideStart = 0
        guideEnd = len(seq1)-len(pam)

    for i in range(guideStart, guideEnd):
        if seq1[i]==seq2[i]:
            lines[0].append(seq1[i])
            lines[1].append(seq2[i])
            lines[2].append(" ")
        else:
            lines[0].append("<b>%s</b>"%seq1[i])
            lines[1].append("<b>%s</b>"%seq2[i])
            lines[2].append("*")
            if i>7:
                last12MmCount += 1

    if not cpf1Mode:
        lines[0].append(" <i>"+seq1[-len(pam):]+"</i>")
        lines[1].append(" <i>"+seq2[-len(pam):]+"</i>")
    lines = ["".join(l) for l in lines]

    if len(posStr)>1 and posStr[0].isdigit():
        posStr = "chr"+posStr

    htmlText1 = "<small><pre>guide:      %s<br>off-target: %s<br>            %s</pre>" \
        % (lines[0], lines[1], lines[2])
    if cpf1Mode:
        htmlText2 = "CPf1: No off-target scores available</small>"
    else:
        if cfdScore==None:
            cfdStr = "Cannot calculate CFD score on non-ACTG characters"
        else:
            cfdStr = "%f" % cfdScore
        htmlText2 = "CFD Off-target score: %s<br>MIT Off-target score: %.2f<br>Position: %s</small>" % (cfdStr, mitScore, posStr)
        if chromDist!=None:
            htmlText2 += "<br><small>Distance from target: %.3f Mbp</small>" % (float(chromDist)/1000000.0)
#             if org.startswith("mm") or org.startswith("hg") or org.startswith("rn"):
#                 if chromDist > 4000000:
#                     htmlText2 += "<br><small>&gt;4Mbp = unlikely to be in linkage with target</small>"
#                 else:
#                     htmlText2 += "<br><small>&lt;4Mbp= likely to be in linkage with target!</small>"

    hasLast12Mm = last12MmCount>0
    return htmlText1+htmlText2, hasLast12Mm
def highlightMismatches(guide, offTarget, pamLen,cpf1Mode=False):
    " return a string that marks mismatches between guide and offtarget with * "
    if cpf1Mode:
        offTarget = offTarget[pamLen:]
    else:
        offTarget = offTarget[:-pamLen]
    assert(len(guide)==len(offTarget))

    s = []
    for x, y in zip(guide, offTarget):
        if x==y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)

def iterOfftargetRows(guideData, offtargetHeaders,pamIdRe,addHeaders=False, skipRepetitive=True, seqId=None):
    " yield bulk offtarget rows for the tab-sep download file "
    otRows = []
    print(offtargetHeaders)
    headers = list(offtargetHeaders) # clone list
#     if seqId:
#         headers.insert(0, "seqId")

    if addHeaders:
        otRows.append(headers)

    skipCount = 0

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, \
            ontargetDesc, subOptMatchCount = guideRow

        if otData!=None:
            otCount = len(otData)
            if otCount > 5000 and skipRepetitive:
                skipCount += otCount
                continue

            for otSeq, mitScore, cfdScore, editDist, pos, gene, alnHtml in otData:
                gene = gene.replace(",", "_").replace(";","-")
                chrom, start, end, strand = parsePos(pos)
                guideDesc = intToExtPamId(pamId,pamIdRe)
                mismStr = highlightMismatches(guideSeq, otSeq, len(pamSeq))
                fullSeq = concatGuideAndPam(guideSeq, pamSeq)
                row = [guideDesc, fullSeq, otSeq, mismStr, editDist, mitScore, cfdScore, chrom, start, end, strand, gene]
                if seqId:
                    row.insert(0, seqId)
                row = [str(x) for x in row]
                otRows.append(row)

    if skipCount != 0:
        newRow = [""]*len(headers)
        newRow[0] = "# %d off-targets are not shown: guides with more than 5000 off-targets were considered too repetitive" % skipCount
        otRows.insert(0, newRow)

    return otRows
#--------------------------------------------
def main(inSeqFname,genomeDir,org,outGuideFname,offtargetFname,genomefn):
    
#     inSeqFname='../../../04_offtarget/data/04_specificity/in/sample.sacCer3.fa',
#     # faFname='tmp/in/x50sPGMoTvUagv3zWjGg.fa'
#     # genomeDir='tmp/in/genomes/'
#     genomeDir='../../../04_offtarget/pub/release-92/fasta',
#     org='saccharomyces_cerevisiae',
#     outGuideFname='../../../04_offtarget/data/04_specificity/out/sample.sacCer3.mine.out.tsv',
#     offtargetFname='../../../04_offtarget/data/04_specificity/out/sample.sacCer3.mine.offs.tsv',
#     genomefn='dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.chromosome.I.fa'
    
# def dguide2fltofftarget(cfg):
    # get dna and protein sequences 
    #--
    # junk
    # write debug output to stdout
    DEBUG = False
    #DEBUG = True
    # calculate the efficienc scores?
    doEffScoring = True

    # system-wide temporary directory
    TEMPDIR = dirname(dirname(offtargetFname))+'/tmp'
    if not exists(TEMPDIR):
        makedirs(TEMPDIR)

    # filename of this script, usually crispor.py
    # myName = basename(__file__)

    # the segments.bed files use abbreviated genomic region names
    segTypeConv = {"ex":"exon", "in":"intron", "ig":"intergenic"}

    # for some PAMs, there are alternative main PAMs. These are also shown on the main sequence panel
    multiPams = {
        "NGN" : ["GAW"]
    }

    offtargetPams = {
        "NGG" : ["NAG","NGA"],
        "NGN" : ["GAW"],
        "NGK" : ["GAW"],
        "NGA" : ["NGG"],
        "NNGRRT" : ["NNGRRN"]
    }

    # BWA: allow up to X mismatches
    maxMMs=4

    # maximum number of occurences in the genome to get flagged as repeats. 
    # This is used in bwa samse, when converting the same file
    # and for warnings in the table output.
    MAXOCC = 60000

    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000/MAXOCC

    # the length of the guide sequence, set by setupPamInfo
    GUIDELEN=None
    # length of the PAM sequence
    PAMLEN=None

    # input sequences are extended by X basepairs so we can calculate the efficiency scores
    # and can better design primers
    FLANKLEN=100

    # the name of the currently processed batch, assigned only once 
    # in readBatchParams and only for json-type batches
    batchName = ""

    # are we doing a Cpf1 run?
    # this variable changes almost all processing and
    # has to be set on program start, as soon as we know 
    # the PAM we're running on
    cpf1Mode=None

    ALTPAMMINSCORE = 1.0

    # how much shall we extend the guide after the PAM to match restriction enzymes?
    pamPlusLen = 5

    # global flag to indicate if we're run from command line or as a CGI

    # names/order of efficiency scores to show in UI
    scoreNames = ["fusi", "crisprScan"]
    # allScoreNames = ["fusi", "fusiOld", "chariRank", "ssc", "doench", "wang", "crisprScan", "aziInVitro"]

    cpf1ScoreNames = ["seqDeepCpf1"]

    saCas9ScoreNames = ["najm"]



    # labels and descriptions of eff. scores
    scoreDescs = {
        "doench" : ("Doench '14", "Range: 0-100. Linear regression model trained on 880 guides transfected into human MOLM13/NB4/TF1 cells (three genes) and mouse cells (six genes). Delivery: lentivirus. The Fusi score can be considered an updated version this score, as their training data overlaps a lot. See <a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench et al.</a>"),
        "ssc" : ("Xu", "Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on data from &gt;1000 genes in human KBM7/HL60 cells (Wang et al) and mouse (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2. See <a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu et al.</a>"),
        "crisprScan" : ["Moreno-Mateos", "Also called 'CrisprScan'. Range: mostly 0-100. Linear regression model, trained on data from 1000 guides on &gt;100 genes, from zebrafish 1-cell stage embryos injected with mRNA. See <a target=_blank href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos et al.</a>. Recommended for guides transcribed <i>in-vitro</i> (T7 promoter). Click to sort by this score."],
        "wang" : ("Wang", "Range: 0-100. SVM model trained on human cell culture data on guides from &gt;1000 genes. The Xu score can be considered an updated version of this score, as the training data overlaps a lot. Delivery: lentivirus. See <a target='_blank' href='http://http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang et al.</a>"),
        "chariRank" : ("Chari", "Range: 0-100. Support Vector Machine, converted to rank-percent, trained on data from 1235 guides targeting sequences that were also transfected with a lentivirus into human 293T cells. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari et al.</a>"),
        "fusi" : ("Doench '16", "Aka the 'Fusi-Score', since V4.4 using the version 'Azimuth', scores are slightly different than before April 2018 but very similar (click 'show all' to see the old scores). Range: 0-100. Boosted Regression Tree model, trained on data produced by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished additional data). Delivery: lentivirus. See <a target='_blank' href='http://biorxiv.org/content/early/2015/06/26/021568'>Fusi et al. 2015</a> and <a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench et al. 2016</a> and <a target=_blank href='https://crispr.ml/'>crispr.ml</a>. Recommended for guides expressed in cells (U6 promoter). Click to sort the table by this score."),
        "fusiOld" : ("Doench '16-Old", "The original implementation of the Doench 2016 score, as received from John Doench. The scores are similar, but not exactly identical to the 'Azimuth' version of the Doench 2016 model that is currently the default on this site, since Apr 2018."),
        "najm" : ("Najm 2018", "A modified version of the Doench 2016 score ('Azimuth'), by Mudra Hegde for S. aureus Cas9. Range 0-100. See <a target=_blank href='https://www.nature.com/articles/nbt.4048'>Najm et al 2018</a>."),
        "aziInVitro" : ("Azimuth in-vitro", "The Doench 2016 model trained on the Moreno-Mateos zebrafish data. Unpublished model, gratefully provided by J. Listgarden"),
        "housden" : ("Housden", "Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA injections. See <a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden et al.</a>"),
        "proxGc" : ("ProxGCCount", "Number of GCs in the last 4pb before the PAM"),
        "seqDeepCpf1" : ("DeepCpf1", "Range: ~ 0-100. Convolutional Neural Network trained on ~20k Cpf1 lentiviral guide results. This is the score without DNAse information, 'Seq-DeepCpf1' in the paper. See <a target='_blank' href='https://www.nature.com/articles/nbt.4061'>Kim et al. 2018</a>"),
        "oof" : ("Out-of-Frame", "Range: 0-100. Predicts the percentage of clones that will carry out-of-frame deletions, based on the micro-homology in the sequence flanking the target site. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al.</a>. Click the score to show the most likely deletions for this guide.")
    }

    # the headers for the guide and offtarget output files
    guideHeaders = ["guideId", "targetSeq", "mitSpecScore", "offtargetCount", "targetGenomeGeneLocus"]
    offtargetHeaders = ["guideId", "guideSeq", "offtargetSeq", "mismatchPos", "mismatchCount", "mitOfftargetScore", "cfdOfftargetScore", "chrom", "start", "end", "strand", "locusDesc"]

    # a file crispor.conf in the directory of the script allows to override any global variable
    myDir = dirname(__file__)
#     confPath =myDir+"crispor.conf"
#     if isfile(confPath):
#         exec(open(confPath))

    # cgiParams = None

    # ====== END GLOBALS ============
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M',
        "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n", "V" : "B", "v":"b", 
        "B" : "V", "b": "v", "W" : "W", "w" : "w"}
    pamIdRe = re.compile(r's([0-9]+)([+-])g?([0-9]*)')

    #
    DEBUG=True
    doEffScoring=False
    genome=org
    pam='NGG'
    genomep='{}/{}/{}'.format(genomeDir,org,genomefn)
    # genomep='crisporWebsite/genomes/sacCer3/sacCer3.fa'
    bedFname='test.bed'
    class options(object):
        def __init__(self, pam):
            self.pam = pam
            self.skipAlign = False
            self.genomeDir=genomeDir
            self.offtargetFname=offtargetFname
            self.debug=DEBUG

    options=options(pam=pam)
    class ConsQueue:
        """ a pseudo job queue that does nothing but report progress to the console """
        def startStep(self, batchId, desc, label):
            logging.info("Progress %s - %s - %s" % (batchId, desc, label))

    effScores='na'
    GUIDELEN,cpf1Mode,addGenePlasmids,PAMLEN,scoreNames=setupPamInfo(pam,setupPamInfo,scoreNames)
    skipAlign = False
    if options.skipAlign:
        skipAlign = True

    # different genomes directory?
    if options.genomeDir != None:
        genomesDir = options.genomeDir

    params_findofftargetbwa={'GUIDELEN':GUIDELEN,
                    'MFAC':MFAC,
                    'MAXOCC':MAXOCC,
                    'offtargetPams':offtargetPams,
                    'ALTPAMMINSCORE':ALTPAMMINSCORE,
                    'DEBUG':DEBUG}

    # get sequence
    seqs = parseFasta(open(inSeqFname))
    guideFh = None
    offtargetFh = None
    batchId='batchId'
    for seqId, seq in seqs.items():
        seq = seq.upper()
#         print(seqId, GUIDELEN, len(seq))
        logging.info(" * running on sequence '%s', guideLen=%d, seqLen=%d" % (seqId, GUIDELEN, len(seq)))
        # get the other parameters and write to a new batch
        seq = seq.upper()
#         pamPat = options.pam
#         batchId, position, extSeq = newBatch(seqId, seq, org, pamPat, genome,TEMPDIR,skipAlign,
#                                              flank=FLANKLEN,genomep=genomep)
#         logging.debug("Temporary output directory: %s/%s" % (TEMPDIR, batchId))

#         if position=="?":
#             logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

#         startDict, endSet = findAllPams(seq, pamPat,GUIDELEN,multiPams)
#        position=get_position(genome, seq, batchId,TEMPDIR=TEMPDIR,genomep=genomep)
        position='chr:1:10' #FIXME
        startDict={33: '+', 40: '+', 57: '+', 103: '+', 136: '+', 79: '-', 139: '-'}
        endSet={36, 106, 43, 139, 142, 82, 60}
#         print('>>>>>')
#         print(position)
#         print(endSet)
#         print('>>>>>')
        otBedFname = getOfftargets(seq, org, options.pam, batchId, startDict,
                                   ConsQueue(),TEMPDIR,GUIDELEN,pamPlusLen,maxMMs,genomep,
#                                   GUIDELEN,
                                    MFAC,
                                    MAXOCC,
                                    offtargetPams,
                                    ALTPAMMINSCORE,
                                    DEBUG)
        otMatches = parseOfftargets(otBedFname)
        guideData, guideScores, hasNotFound, pamIdToSeq = \
            mergeGuideInfo(seq, startDict, options.pam, otMatches, position, GUIDELEN, pamPlusLen, maxMMs,cpf1Mode,params_findofftargetbwa,effScores)

        # write guide headers
        if guideFh is None:
            guideFh = open(join(TEMPDIR, "guideInfo.tab"), "w")
            guideHeaders, _ = makeGuideHeaders(guideHeaders,scoreNames,
                     cpf1Mode,scoreDescs)
            guideHeaders.insert(0, "#seqId")
            guideFh.write("\t".join(guideHeaders)+"\n")

    #     write offtarget headers
        if options.offtargetFname and (offtargetFh is None):
            print('# write offtarget headers')
            offtargetFh = open(join(TEMPDIR, "offtargetInfo.tab"), "w")
            offtargetHeaders.insert(0, "seqId")
            offtargetFh.write("\t".join(offtargetHeaders)+"\n")
        print(offtargetHeaders)
        for row in iterGuideRows(guideData, guideHeaders,scoreNames,
                     cpf1Mode,scoreDescs,pamIdRe,seqId=seqId):
            guideFh.write("\t".join(row))
            guideFh.write("\n")

        if options.offtargetFname:
            for row in iterOfftargetRows(guideData,offtargetHeaders,pamIdRe, seqId=seqId, skipRepetitive=False):
                offtargetFh.write("\t".join(row))
                offtargetFh.write("\n")

    guideFh.close()
    shutil.move(guideFh.name, outGuideFname)
    logging.info("guide info written to %s" % outGuideFname)

    if options.offtargetFname:
        offtargetFh.close()
        shutil.move(offtargetFh.name, options.offtargetFname)
        logging.info("off-target info written to %s" % options.offtargetFname)

    if not options.debug and not options.tempDir:
        shutil.rmtree(TEMPDIR)

# if __name__ == '__main__':
#     main()
