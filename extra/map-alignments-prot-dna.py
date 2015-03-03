#!/usr/bin/env python
# -*- coding: utf8 -*-
# Conversion tool for LAST MAF files
# Copyright (C) 2010 Johannes Dröge
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from sys import argv, stdin, stderr, stdout
from common import ProteinDnaMapper, AlignmentRecord

def Usage():
	print >> stderr, 'Program reads from standard input and writes to standard output.'
	print >> stderr, 'Usage: ', argv[0], '[--help | --mapping map.gff3 ]'

## standalone program reads from stdin
if __name__ == "__main__":
	# handle broken pipes
	import signal
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)
	
	# parse command line options
	import getopt
	try:
		opts, args = getopt.getopt( argv[1:], 'hm:', ['help', 'mapping='] )
	except getopt.GetoptError, err:
		# print help information and exit:
		print str( err ) # will print something like "option -a not recognized"
		Usage()
		exit(1)
		
	mf = None

	for o, a in opts:
		if o in ("-h", "--help"):
			Usage()
			exit( 1 )
		elif o in ("-m", "--mapping"):
			mf = a
		else:
			assert False, "unhandled option"

	if not mf:
		Usage()
		exit( 1 )

	mapper = ProteinDnaMapper( mf )
	record = AlignmentRecord()
	record.printHeader()
	for line in stdin:
		try:
			if line[0] != "#" and record.parseRecord( line ):
				if record.query_start < record.query_stop:
					record.reference_identifier, record.reference_start, record.reference_stop = mapper.mapPosition( record.reference_identifier, record.reference_start, record.reference_stop )
				else:
					# swap ranges such that query positions are ordered
					record.query_start, record.query_stop = record.query_stop, record.query_start
					record.reference_identifier, record.reference_stop, record.reference_start = mapper.mapPosition( record.reference_identifier, record.reference_start, record.reference_stop )
				record.identities *= 3
				record.alignment_length *= 3
				record.printRecord()
		except KeyError:
			print >> stderr, "Could not map protein GI %s" % (record.reference_identifier)
			pass
