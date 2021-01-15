# -*- coding: utf-8 -*-
# taxator-tk predicts the taxon for DNA sequences based on sequence alignment.

#Copyright (C) 2010 Johannes Dröge

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""This module acts as a python bridge to work with taxator-tk-style
data formats"""

# common packages etc.
from itertools import izip
from math import isnan


def dot2nan(x):
    if x == ".":
        return float("nan")
    return float(x)


def nan2dot(x):
    if isnan(x):
        return "."
    return str(x)  #TODO: precision?


class TaxonFeatureRecord(object):
    """class to represent taxonomic features
    read and write to GFF3 taxator-tk format"""

    def __init__(self, line):
        self.__parse_gff3_line(line)

    def __get_gff3_line(self):
        tax = "-".join(["%s:%i" % (t, s) for t, s in self.taxon_support])  #TODO: avoid writing support if equals length
        rest = ";".join(["%s=%s" % (k, v) for k, v in self.attributes.iteritems()])
        if len(rest) > 0:
            sep = ";"
        else:
            sep = ""
        return "%s\t%s\tsequence_feature\t%i\t%i\t%s\t.\t.\tseqlen=%i;tax=%s%s%s" % (
        self.sequence_identifier, self.generator, self.feature_begin, self.feature_end, nan2dot(self.feature_score),
        self.sequence_length, tax, sep, rest)

    def __parse_gff3_line(self, line):
        columns = line.split("\t")
        self.sequence_identifier = columns[0]
        self.generator = columns[1]
        self.feature_begin = int(columns[3])
        self.feature_end = int(columns[4])
        self.feature_score = dot2nan(columns[5])

        attributes = dict()
        for entry in columns[8].strip().split(";"):
            key, val = entry.split("=")
            attributes[key] = val

        # manditory attribute "seqlen"
        self.sequence_length = int(attributes.pop("seqlen"))

        # manditory attribute "tax"
        taxon_support = attributes.pop("tax").split("-")
        taxid = taxon_support[0].split(":")
        if len(taxid) > 1:
            support = int(taxid[1])
        else:
            support = self.sequence_length
        taxid = taxid[0]

        self.taxon_support = [(taxid, support)]
        for entry in taxon_support[1:]:
            taxid = entry.split(":")
            if len(taxid) > 1:
                support = int(taxid[1])
            taxid = taxid[0]
            self.taxon_support.append((taxid, support))

        self.attributes = attributes

    line = property(__get_gff3_line, __parse_gff3_line)

    __repr__ = lambda self: self.line


__regular_dna_lowercase = frozenset(["a", "c", "g", "t"])


class AlignmentCode(object):
    """class to generate AlignmentCode from complete pairwise alignment (like CIGAR)"""

    def __init__(self):
        self.symbol = ""
        self.duplication = ""
        self.code = []

    def append(self, s):
        if s == self.symbol_:
            self.duplication_ += 1
        else:
            if self.duplication == 1:
                self.code_.append(self.symbol_)
            else:
                self.code_.append(str(self.duplication_) + self.symbol_)
            self.symbol = s
            self.duplication = 1

    @property
    def code(self):
        return "".join(self.code_) + str(self.duplication_) + self.symbol_

    def clear(self):
        self.symbol = ""
        self.duplication = ""
        del self.code_[:]

    __repr__ = lambda self: self.code


def get_alignment_encoding_dna(seq, ref):
    """returns an alignment encoding for two aligned sequences"""
    aln = AlignmentCode()
    for s, r in izip(seq, ref):
        if s == "-":
            if r != "-":
                aln.append("D")
            else:
                continue
        elif r == "-":
            aln.append("I")
        elif r != s:
            aln.append("X")  # CIGAR extension for mismatch pair
        else:
            aln.append("=")  # CIGAR extension for match pair
    return aln.get()


# compatability function, returns empty alignment encoding
get_alignment_encoding_dummy = lambda seq, ref: ""


# this is a stub, it will be extended to generate alignment encodings for protein
get_alignment_encoding_prot = get_alignment_encoding_dna


class AlignmentScoreCalculatorDNA(object):
    """class for calculation of a pairwise alignment score with affine gap costs"""

    def __init__(self, match, mismatch, gapopen, gapext):
        self.__match = match
        self.__mismatch = mismatch
        self.__gapopen = gapopen
        self.__gapext = gapext

    # other symbols are a mismatch, always! this follows BLAST practice
    def __call__(self, line1, line2):
        score = 0.0
        is_gapopen = True
        for q, s in izip(line1.lower(), line2.lower()):
            if q == s and q in __regular_dna_lowercase and s in __regular_dna_lowercase:
                score += self.__match
                is_gapopen = True
            else:
                if q != "-" and s != "-":
                    score += self.__mismatch
                    is_gapopen = True
                else:
                    if is_gapopen:
                        score -= self.__gapopen
                        is_gapopen = False
                    score -= self.__apext
        return score


class AlignmentScoreCalculatorBLOSUM62(object):
    """class for calculation of pairwise alignments scores using BLOSUM62 (proteins)"""
    __general_mismatch = -4
    __matrix = {
        ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
        ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
        ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
        ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
        ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
        ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
        ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
        ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
        ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
        ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
        ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
        ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
        ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
        ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
        ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
        ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
        ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
        ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
        ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
        ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
        ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
        ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
        ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
        ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
        ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
        ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
        ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
        ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
        ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
        ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
        ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
        ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
        ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
        ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
        ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
        ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
        ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
        ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
        ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
        ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
        ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
        ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
        ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
        ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
        ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
        ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
        ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
        ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
        ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
        ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
        ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
        ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
        ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
        ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
        ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
        ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
        ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
        ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
        ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
        ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
        ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
        ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
        ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
        ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
        ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
        ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
        ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
        ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
        ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
    }

    def __init__(self, gapopen, gapext):
        self.__gapopen = gapopen
        self.__gapext = gapext

    def __call__(self, line1, line2):
        score = 0.0
        for q, s in izip(line1.upper(), line2.upper()):
            if q == "-":
                if s != "-":
                    if is_gapopen:
                        score -= self.__gapopen
                    else:
                        score -= self.__gapext
                    is_gapopen = False
            elif s == "-":
                if q != "-":
                    if is_gapopen:
                        score -= self.__gapopen
                    else:
                        score -= self.__gapext
                    is_gapopen = False
            else:
                try:
                    if q < s:
                        score += self.__matrix[(q, s)]
                    else:
                        score += self.__matrix[(s, q)]
                except KeyError:
                    score += self.__general_mismatch
                is_gapopen = True
        return score


def extract_gi_from_identifier(string):
    """extracts the GI number from a NCBI-style fasta header"""
    try:
        target_id = string.replace('|', ' ').split()
        for i in range(0, len(target_id)):
            if target_id[i][-2:] == 'gi':
                return target_id[i + 1]
    except IndexError:
        pass
    return None


class AlignmentRecord(object):
    """class to represent alignment records
    read and write in taxator-tk tab-separated format"""

    def __init__(self, line=None):
        if line:
            self.__parse_record_line(line)
            return

        self.query_identifier = None
        self.query_start = None
        self.query_stop = None
        self.query_length = None
        self.reference_identifier = None
        self.reference_start = None
        self.reference_stop = None
        self.score = None
        self.evalue = None
        self.identities = None
        self.alignment_length = None
        self.alignment_code = None
        self.blacklist_this = False

    def __parse_record_line(self, line):
        try:
            if line[0] == "*":
                self.blacklist_this = True
                fields = line[1:].strip().split(self.__field_separator)
            else:
                self.blacklist_this = False
                fields = line.strip().split(self.__field_separator)

            self.query_identifier = fields[0]
            self.query_start = int(fields[1])
            self.query_stop = int(fields[2])
            self.query_length = int(fields[3])
            self.reference_identifier = fields[4]
            self.reference_start = int(fields[5])
            self.reference_stop = int(fields[6])
            self.score = float(fields[7])
            self.evalue = float(fields[8])
            self.identities = int(fields[9])
            self.alignment_length = int(fields[10])
            if len(fields) > 11:
                self.alignment_code = fields[11]
            else:
                self.alignment_code = ""
            return True

        except (IndexError, ValueError):
            return False

    def __get_record_line(self):
        if self.blacklist_this:
            masksymbol = "*"
        else:
            masksymbol = ""
        return "%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        masksymbol, self.query_identifier, self.query_start, self.query_stop, self.query_length,
        self.reference_identifier, self.reference_start, self.reference_stop, self.score, self.evalue, self.identities,
        self.alignment_length, self.alignment_code )

    @property
    def file_header(self):
        return "#%s\n" % (self.__field_separator.join(self.__field_names))

    __field_names = ["query ID", "query begin", "query end", "query length", "reference ID", "reference begin",
                     "reference end", "score", "E-value", "identities", "alignment length", "alignment code"]
    __field_separator = "\t"

    line = property(__get_record_line, __parse_record_line)

    __repr__ = lambda self: self.line



class ProteinDnaMapper(object):
    """class takes a taxator-tk mapping file and creates an index
    to map positions from dna to protein"""

    def __init__(self, filename, ftype="CDS"):
        self.__table = {}
        for line in open(filename, "r"):
            if line[0] == "#":
                continue

            fields = line.strip().split("\t")

            if fields[2] != ftype:
                continue

            for key, val in [map(str.strip, f.split("=")) for f in fields[8].split(";")]:
                if key == "ID":
                    protein_id = val
                    break
            if fields[6] == "+":
                self.__table[protein_id] = ( fields[0], int(fields[3]), int(fields[4]) )
            elif fields[6] == "-":
                self.__table[protein_id] = ( fields[0], int(fields[4]), int(fields[3]) )
            else:
                raise SyntaxError("cannot handle other strand symbols than + and - for protein sequences")

    def map_position(self, pidentifier, pos1, pos2):
        didentifier, dstart, dstop = self.__table[pidentifier]  #TODO: catch KeyError
        if dstart < dstop:
            pos1 = dstart + 3 * (pos1 - 1)
            pos2 = dstart + 3 * (pos2 - 1) + 2
        else:
            pos1 = dstart - 3 * (pos1 - 1)
            pos2 = dstart - 3 * (pos2 - 1) - 2
        return didentifier, pos1, pos2
