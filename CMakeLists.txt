PROJECT(taxator-tk)
cmake_minimum_required( VERSION 2.6 )

# internal variables
set( boost_version_required 1.34.1 )
set( boost_version_required_enc 103401 )
set( project_name taxator-tk )

# cmake requirements detection
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)
# find_package( Sqlite3 )
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost "${boost_version_required}" REQUIRED COMPONENTS program_options system filesystem thread )
if( Boost_FOUND )
	if( Boost_VERSION LESS "${boost_version_required_enc}" )
    message( FATAL_ERROR "FOUND INCORRECT BOOST VERSION: ${Boost_LIB_VERSION},\n${project_name} requires at least version ${boost_version_required} ")
  else()
		include_directories(${Boost_INCLUDE_DIRS})
		link_directories(${Boost_LIBRARY_DIRS})
	endif()
else()
  message( FATAL_ERROR "Unable to find a Boost version. Did you set BOOST_ROOT?" )
endif()

include_directories( "includes-external" )
set(CMAKE_CXX_FLAGS "-Wall -pedantic -Wno-long-long -Wno-variadic-macros -fpermissive -O2 -march=native") #-g for debuggin, -m32 for x32

# traverse nodes up to a given rank
#add_executable( rank-filter rank-filter.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp )
#target_link_libraries( rank-filter ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# looks up the name or rank name of a node
#add_executable( name-filter name-filter.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp )
#target_link_libraries( name-filter ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# apply filtering to alignments file
add_executable( alignments-filter alignments-filter.cpp src/alignmentrecord.cpp src/accessconv.cpp )
target_link_libraries( alignments-filter ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# takes input alignments and predicts a taxon for each query id using various methods and parameters
add_executable( taxator taxator.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp src/accessconv.cpp src/predictionrecord.cpp )
#target_link_libraries( taxator boost_program_options boost_system boost_filesystem boost_thread )
target_link_libraries( taxator ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} )

# apply filtering to predictions file
add_executable( binner binner.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp src/predictionrecord.cpp )
target_link_libraries( binner ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# taxknife 
add_executable( taxknife taxknife.cpp src/taxontree.cpp src/taxonomyinterface.cpp src/ncbidata.cpp )
target_link_libraries( taxknife ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# The following targets are usually not needed or for development and testing: uncomment if needed

# general converter for string or id mappings
#add_executable( acc2taxid acc2taxid.cpp src/accessconv.cpp )
#target_link_libraries( acc2taxid ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# parses NCBI-style identifiers in the FASTA comment header
#add_executable( extract-fastacomment-ncbifield extract-fastacomment-ncbifield.cpp src/ncbidata.cpp src/taxontree.cpp )
#target_link_libraries( extract-fastacomment-ncbifield ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )

# unittest: constructs the taxonomy from NCBI dump files and tests the structure thoroughly
#add_executable( unittest_ncbitaxonomy unittest_ncbitaxonomy.cpp src/ncbidata.cpp src/accessconv.cpp src/taxontree.cpp src/taxonomyinterface.cpp )
#target_link_libraries( unittest_ncbitaxonomy ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )
