#!/usr/bin/env python
# This script will produce NCBI Blast tab-like format from taxator-tk
# intermediate alignments format, expecting sorted input (same query in block)
#
# TODO: set database parameter
# TODO: eliminate undefined columns

from taxatortk import AlignmentRecord
from itertools import chain, takewhile

class simpleCache( object ):
	def __init__( self, generator ):
		self._gen = generator
		self._cache = None
		self._notempty = True
	
	def __iter__( self ):
		for i in self._gen:
			self._cache = i
			yield i
		self._notempty = False
	
	__nonzero__ = lambda self: self._notempty
	last = lambda self: self._cache

def parseAlignmentRecords( fh ):
		records = simpleCache( ( AlignmentRecord( line ) for line in fh if line and line[0] != "#" ) )
		current = iter( records ).next()
		while records:
			yield current, chain( [current], takewhile( lambda f: f.query_identifier == current.query_identifier, records ) )
			current = records.last()

header = "\
# taxator\n\
# Query: %s\n\
# Database: unknown\n\
# Fields: Query id, Subject id, %% identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score"

if __name__ == "__main__":
	from sys import stdin
	
	# handle broken pipes
	from signal import signal, SIGPIPE, SIG_DFL
	signal( SIGPIPE, SIG_DFL )
	
	# process and convert alignment records
	for rec, block in parseAlignmentRecords( stdin ):
		print header % rec.query_identifier
		for rec in block:
			mismatches = rec.alignment_length - rec.identities
			if rec.alignment_code: 
				gap_openings = rec.alignment_code.count("I") + rec.alignment_code.count("D") #TODO: check this simple formula
			else:
				gap_openings = 0 #set to zero as placeholder; this is not correct!
			print "%s\t%s\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%0.1f" % (rec.query_identifier,rec.reference_identifier,100*rec.identities/float( rec.alignment_length ),rec.alignment_length,mismatches,gap_openings,rec.query_begin,rec.query_end,rec.reference_begin,rec.reference_end,rec.evalue,rec.score)
