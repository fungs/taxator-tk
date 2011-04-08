PROJECT(ttp)

cmake_minimum_required( VERSION 2.6 )

# find_package(sqlite3)
set(CMAKE_CXX_FLAGS "-g -Wall")
# set(CMAKE_CXX_FLAGS "-m32")

# unittest: constructs the taxonomy from NCBI dump files and tests the structure thoroughly
add_executable( unittest_ncbitaxonomy unittest_ncbitaxonomy.cpp src/ncbidata.cpp src/accessconv.cpp src/taxontree.cpp src/sqlite3pp.cpp )
target_link_libraries( unittest_ncbitaxonomy sqlite3 )

# traverse nodes up to a given rank
add_executable( rank-filter rank-filter.cpp src/taxontree.cpp src/ncbidata.cpp )
target_link_libraries( rank-filter boost_program_options )

# looks up the name or rank name of a node
add_executable( name-filter name-filter.cpp src/taxontree.cpp src/ncbidata.cpp )
target_link_libraries( name-filter boost_program_options )

# general converter for string or id mappings
add_executable( acc2taxid acc2taxid.cpp src/accessconv.cpp src/sqlite3pp )
target_link_libraries( acc2taxid sqlite3 boost_program_options )

# parses NCBI-style identifiers in the FASTA comment header
add_executable( extract-fastacomment-ncbifield extract-fastacomment-ncbifield.cpp src/ncbidata.cpp src/taxontree.cpp )
target_link_libraries( extract-fastacomment-ncbifield boost_program_options )

# apply filtering to alignments file
add_executable( alignments-prefilter alignments-prefilter.cpp src/alignmentrecord.cpp )
target_link_libraries( alignments-prefilter boost_program_options )

# takes input alignments and predicts a taxon for each query id using various methods and parameters
add_executable( predictor predictor.cpp src/taxontree.cpp src/ncbidata.cpp src/accessconv.cpp src/sqlite3pp.cpp src/alignmentrecord.cpp src/taxonprediction.cpp )
target_link_libraries( predictor sqlite3 boost_program_options )

# measures the prediction error in terms of tree path distances
add_executable( prediction2distances prediction2distances.cpp src/taxontree.cpp src/ncbidata.cpp )
target_link_libraries( prediction2distances boost_program_options )

# various formulas for error calculation
add_executable( prediction-error prediction-error.cpp )
target_link_libraries( prediction-error boost_program_options )

# filter alignments (supervised mode, needed for benchmarks only)
add_executable( alignments-prefilter-training alignments-prefilter-training.cpp src/alignmentrecord.cpp src/accessconv.cpp src/sqlite3pp src/taxontree.cpp src/ncbidata.cpp )
target_link_libraries( alignments-prefilter-training boost_program_options sqlite3 )
