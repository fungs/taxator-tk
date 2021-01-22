# This file was part of LAST v454 and released under the GNU GPL3 or later
# This is a partial extract and included as a python module to avoid rewriting
# functions for parsing MAF files produced by LAST

from itertools import *
import sys, os, fileinput, math, operator, optparse, signal, string

def quantify(iterable, pred=bool):
    """Count how many times the predicate is true."""
    return sum(map(pred, iterable))

class Maf:
    def __init__(self, aLine, sLines, qLines, pLines):
        self.namesAndValues = dict(i.split("=") for i in aLine[1:])

        if not sLines: raise Exception("empty alignment")
        cols = zip(*sLines)
        self.seqNames = cols[1]
        self.alnStarts = map(int, cols[2])
        self.alnSizes = map(int, cols[3])
        self.strands = cols[4]
        self.seqSizes = map(int, cols[5])
        self.alnStrings = cols[6]

        self.qLines = qLines
        self.pLines = pLines

def dieUnlessPairwise(maf):
    if len(maf.alnStrings) != 2:
        raise Exception("pairwise alignments only, please")

def isMatch(alignmentColumn):
    # No special treatment of ambiguous bases/residues: same as NCBI BLAST.
    first = alignmentColumn[0].upper()
    for i in alignmentColumn[1:]:
        if i.upper() != first: return False
    return True

def mafInput(lines, isKeepComments):
    a = []
    s = []
    q = []
    p = []
    for i in lines:
        w = i.split()
        if not w:
            if a: yield a, s, q, p
            a = []
            s = []
            q = []
            p = []
        elif w[0] == "a":
            a = w
        elif w[0] == "s":
            s.append(w)
        elif w[0] == "q":
            q.append(w)
        elif w[0] == "p":
            p.append(w)
        elif i[0] == "#" and isKeepComments:
            print i,
    if a: yield a, s, q, p
