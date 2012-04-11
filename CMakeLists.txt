PROJECT(taxator-tk)
cmake_minimum_required( VERSION 2.6 )

# internal variables
set( boost_version_required 1.34.1 )
set( boost_version_required_enc 103401 )
set( project_name taxator-tk )

# cmake requirements detection
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)
find_package( Sqlite3 )
find_package( Boost "${boost_version_required}" REQUIRED COMPONENTS program_options filesystem )
if( Boost_FOUND )
	if( Boost_VERSION LESS "${boost_version_required_enc}" )
    message( FATAL_ERROR "FOUND INCORRECT BOOST VERSION: ${Boost_LIB_VERSION},\n${project_name} requires at least version ${boost_version_required} ")
  else()
# 		#SET(Boost_USE_STATIC_LIBS OFF)
# 		#SET(Boost_USE_MULTITHREADED ON)
# 		#SET(Boost_USE_STATIC_RUNTIME OFF)
	endif()
else()
  message( FATAL_ERROR "Unable to find a Boost version. Did you set BOOST_ROOT?" )
endif()
include_directories( "includes-external" )
set(CMAKE_CXX_FLAGS "-g -Wall -pedantic -Wno-long-long -Wno-variadic-macros") #add -fopenmp for OpenMP

# set(CMAKE_CXX_FLAGS "-m32")

# unittest: constructs the taxonomy from NCBI dump files and tests the structure thoroughly
add_executable( unittest_ncbitaxonomy unittest_ncbitaxonomy.cpp src/ncbidata.cpp src/accessconv.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/sqlite3pp.cpp )
target_link_libraries( unittest_ncbitaxonomy sqlite3 boost_filesystem )

# traverse nodes up to a given rank
add_executable( rank-filter rank-filter.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp )
target_link_libraries( rank-filter boost_program_options boost_filesystem )

# looks up the name or rank name of a node
add_executable( name-filter name-filter.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp )
target_link_libraries( name-filter boost_program_options boost_filesystem )

# general converter for string or id mappings
add_executable( acc2taxid acc2taxid.cpp src/accessconv.cpp src/sqlite3pp )
target_link_libraries( acc2taxid sqlite3 boost_program_options boost_filesystem )

# parses NCBI-style identifiers in the FASTA comment header
add_executable( extract-fastacomment-ncbifield extract-fastacomment-ncbifield.cpp src/ncbidata.cpp src/taxontree.cpp )
target_link_libraries( extract-fastacomment-ncbifield boost_program_options boost_filesystem )

# apply filtering to alignments file
add_executable( alignments-prefilter alignments-prefilter.cpp src/alignmentrecord.cpp )
target_link_libraries( alignments-prefilter boost_program_options boost_system )

# takes input alignments and predicts a taxon for each query id using various methods and parameters
add_executable( predictor predictor.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp src/accessconv.cpp src/sqlite3pp.cpp src/taxonpredictionmodel.cpp src/predictionrecord.cpp )
target_link_libraries( predictor sqlite3 boost_program_options boost_filesystem boost_thread )

# apply filtering to predictions file
add_executable( predictions-filter predictions-filter.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp src/predictionrecord.cpp )
target_link_libraries( predictions-filter boost_program_options boost_filesystem )

# measures the prediction error in terms of tree path distances
add_executable( prediction2distances prediction2distances.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp )
target_link_libraries( prediction2distances boost_program_options boost_filesystem )

# various formulas for error calculation
add_executable( prediction-error prediction-error.cpp )
target_link_libraries( prediction-error boost_program_options )

# filter alignments (supervised mode, needed for benchmarks only)
add_executable( alignments-prefilter-training alignments-prefilter-training.cpp src/accessconv.cpp src/sqlite3pp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp src/alignmentrecord.cpp )
target_link_libraries( alignments-prefilter-training boost_program_options sqlite3 boost_filesystem boost_system boost_regex )

#debugging
add_executable( testpredictionrecord testpredictionrecord.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp src/predictionrecord.cpp )
target_link_libraries( testpredictionrecord boost_filesystem )