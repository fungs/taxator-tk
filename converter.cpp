/*
taxator-tk predicts the taxon for DNA sequences based on sequence alignment.

Copyright (C) 2010 Johannes Dröge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*! \file converter.cpp
 * A program that provides conversion between different prediction formats.
 *
 * \author Johannes Dröge <johannes.droege@uni-duesseldorf.de>
 */

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/scoped_ptr.hpp>
#include <assert.h>
#include "src/taxonomyinterface.hh"
#include "src/ncbidata.hh"
#include "src/utils.hh"
#include "src/constants.hh"

using namespace std;

int main( int argc, char** argv ) {
// 	string oformat_args;
	vector< string > cami_headers, oformat_args;
	unsigned int ident_pos, taxon_pos;

	namespace po = boost::program_options;

	// general options
	po::options_description general("Allowed options");
	general.add_options()
	( "help,h", "show help message")
	( "identifier-field,i", po::value< unsigned int >( &ident_pos )->default_value( 1 ), "input column of sequence identifier" )
	( "taxon-field,t", po::value< unsigned int >( &taxon_pos )->default_value( 2 ), "input column of taxon identifier" )
	( "output-format,o", po::value< vector<string> >( &oformat_args)->multitoken()->required(), "output format specification");
	
	
	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( general ).run(), vm);
	po::notify(vm);

	// global sanity checks
	if( vm.count("help")) { cout << general << endl; return EXIT_SUCCESS; }
	if( ident_pos < 1 ) { cerr << "Field number index is 1-based" << endl; return EXIT_FAILURE; }
	if( taxon_pos < 1 ) { cerr << "Field number index is 1-based" << endl; return EXIT_FAILURE; }
	
	// option parser for CAMI format //TODO: I'd like to pass following arguments to this option parser
// 	po::options_description cami("Options for CAMI output");
// 	cami.add_options()
// 	( "header,h", po::value< vector <string> >( &cami_headers )->multitoken(), "Additional CAMI header in key=value syntax\n" );
	
	if( oformat_args[0] == "CAMI" ) { 
		// add CAMI key=value headers
		for( uint i = 1; i <= oformat_args.size(); ++i ) {
			cout << '@' << oformat_args[i] << endl;
		}
		
		cout << "@@SEQUENCEID\tTAXID" << endl;
		
		// parse line by line
		uint max_pos = max( ident_pos, taxon_pos );
		--ident_pos;
		--taxon_pos;
		const string cami_field_separator = "\t";
		vector< string > fields;
		string line;
		while( getline( cin, line ) ) {
			if ( ignoreLine( line ) ) continue;
			else {
				tokenizeSingleCharDelim( line, fields, cami_field_separator, max_pos );
				cout << fields[ident_pos] << cami_field_separator << fields[taxon_pos] << endl;
			}
		}
	} else {
		cerr << "Output format " << oformat_args[0] << " is unrecognized" << endl;
	}

	return EXIT_SUCCESS;
}
