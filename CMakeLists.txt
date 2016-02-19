#########################################################################################
# Header Message
#
MESSAGE(STATUS "+--------------------------------------------------------------------+")
MESSAGE(STATUS "|                         RepeatDetecter                             |")
MESSAGE(STATUS "+--------------------------------------------------------------------+")
MESSAGE(STATUS "|     (C) Copyright SIB Swiss Institute of Bioinformatics            |")
MESSAGE(STATUS "|         Thierry Schuepbach (thierry.schuepbach@sib.swiss)          |")
MESSAGE(STATUS "|                                                                    |")
MESSAGE(STATUS "+--------------------------------------------------------------------+")

#########################################################################################
# CMAKE
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
PROJECT(RepeatDetecter C CXX)
IF(COMMAND cmake_policy)
  IF(POLICY CMP0003) 
	CMAKE_POLICY(SET CMP0003 NEW)
  ENDIF(POLICY CMP0003)
  IF(POLICY CMP0053) 
	CMAKE_POLICY(SET CMP0053 NEW)
  ENDIF(POLICY CMP0053) 
  IF(POLICY CMP0042)
	CMAKE_POLICY(SET CMP0042 NEW)
  ENDIF(POLICY CMP0042)
ENDIF(COMMAND cmake_policy)
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#########################################################################################
# Version with GIT
#
SET(MAJOR_VERSION "1")
SET(MINOR_VERSION "0")
SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

FIND_PACKAGE(Git REQUIRED)
cmake_host_system_information(RESULT buildhost QUERY FQDN)
FILE(WRITE ${CMAKE_BINARY_DIR}/Version.h.in
"\#define PROJECT_MAJOR_VERSION ${MAJOR_VERSION}\n"
"\#define PROJECT_MINOR_VERSION ${MINOR_VERSION}\n"
"\#define PROJECT_PATCH_VERSION @PATCH_VERSION@\n"
"\#define PROJECT_VERSION \"${MAJOR_VERSION}.${MINOR_VERSION}-@PATCH_VERSION@\"\n"
"\#define PROJECT_VERSION_FULL_STRING \"RepeatDecoder version ${MAJOR_VERSION}.${MINOR_VERSION}-@PATCH_VERSION@ built \" __DATE__ \", \" __TIME__ \" on ${buildhost}\"\n"
)
FILE(WRITE ${CMAKE_BINARY_DIR}/version.cmake
"EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} rev-parse --verify --short HEAD
     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
     OUTPUT_VARIABLE PATCH_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
 )
 CONFIGURE_FILE(\${SRC} \${DST} @ONLY)
")
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
ADD_CUSTOM_TARGET(
    git_version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_BINARY_DIR}/Version.h.in
                     -D DST=${CMAKE_BINARY_DIR}/Version.h
                     -P ${CMAKE_BINARY_DIR}/version.cmake
)

EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} rev-parse --verify --short HEAD
     OUTPUT_VARIABLE PATCH_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
)

MESSAGE(STATUS "This is patch version ${PATCH_VERSION}")

#########################################################################################
# Architectures
#
INCLUDE(CMakeDetermineSystem)
#INCLUDE(CMakePrintSystemInformation)
IF(CMAKE_SYSTEM_PROCESSOR)
    MESSAGE(STATUS "Compilation on architecture ${CMAKE_SYSTEM_PROCESSOR}.")
ENDIF(CMAKE_SYSTEM_PROCESSOR)

#########################################################################################
# What type of build
#
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#########################################################################################
# Determine C Compiler
#
IF(NOT CMAKE_COMPILER_IS_GNUCC)
    EXEC_PROGRAM(${CMAKE_C_COMPILER} ARGS ${CMAKE_BOOT_C_FLAGS} -E "\"${PFTOOLS_SOURCE_DIR}/cmake/CMakeTestIntel.c\"" OUTPUT_VARIABLE CMAKE_COMPILER_OUTPUT RETURN_VALUE CMAKE_COMPILER_RETURN)
    SET(CMAKE_COMPILER_IS_GNUCC_RUN 1)
    IF(NOT CMAKE_COMPILER_RETURN)
        IF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*THIS_IS_INTEL.*" )
             INCLUDE(${PFTOOLS_SOURCE_DIR}/cmake/Linux-icc.cmake)
        ENDIF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*THIS_IS_INTEL.*" )
    ENDIF(NOT CMAKE_COMPILER_RETURN)
ENDIF(NOT CMAKE_COMPILER_IS_GNUCC)

INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG(-msse HANDLE_SSE)
IF(HANDLE_SSE)
  SET(CMAKE_C_SSE_FLAGS "-msse")
ELSE(HANDLE_SSE)
  MESSAGE(FATAL_ERROR "C compiler does not support SSE")
ENDIF(HANDLE_SSE)
CHECK_C_COMPILER_FLAG(-msse2 HANDLE_SSE2)
IF(HANDLE_SSE2)
  SET(CMAKE_C_SSE2_FLAGS "-msse2")
ELSE(HANDLE_SSE2)
  MESSAGE(FATAL_ERROR "C compiler does not support SSE 2")
ENDIF(HANDLE_SSE2)
SET(SSE41_FLAGS "-msse4.1" "-msse41")
UNSET(CMAKE_C_SSE41_FLAGS)
FOREACH(flag IN LISTS SSE41_FLAGS)
  MESSAGE(STATUS "testing flag ${flag}...")
  CHECK_C_COMPILER_FLAG("${flag}" HANDLE_SSE41)
  IF(HANDLE_SSE41)
    SET(CMAKE_C_SSE41_FLAGS "${flag}")
    BREAK()
  ENDIF(HANDLE_SSE41)
ENDFOREACH(flag)
IF(NOT DEFINED CMAKE_C_SSE41_FLAGS)
  MESSAGE(FATAL_ERROR "C compiler does not support SSE 4.1")
ENDIF(NOT DEFINED CMAKE_C_SSE41_FLAGS)

CHECK_C_COMPILER_FLAG("-std=c99" HANDLE_C99)
IF(NOT HANDLE_C99)
  MESSAGE(FATAL_ERROR "PFTOOLS require a C compiler with C99 support")
ELSE(NOT HANDLE_C99)
  IF(CMAKE_C_FLAGS)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
  ELSE(CMAKE_C_FLAGS)
		SET(CMAKE_C_FLAGS "-std=c99")
  ENDIF(CMAKE_C_FLAGS)
  MESSAGE(STATUS "Add -std=c99 to C compiler options")
ENDIF(NOT HANDLE_C99)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_SSE2_FLAGS}")
MESSAGE(STATUS "Add SSE2 to C compiler options")

IF(CMAKE_COMPILER_IS_GNUCC)
  CHECK_C_COMPILER_FLAG("-fms-extensions" MS_EXTENSION)
  IF(MS_EXTENSION)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fms-extensions")
  ENDIF(MS_EXTENSION)
ENDIF(CMAKE_COMPILER_IS_GNUCC)

#########################################################################################
# Determine Fortran Compiler
#

#########################################################################################
# External Dependencies
#
SET(FIND_LIBRARY_USE_LIB64_PATHS 1)

# Check for pthread
FIND_PACKAGE(Threads REQUIRED)
IF( ${CMAKE_USE_WIN32_THREADS_INIT})
  MESSAGE(STATUS "Using Windows API for threads")
  SET(USE_WINAPI 1)
ENDIF(${CMAKE_USE_WIN32_THREADS_INIT})

# Check for presence of Perl for profile generation
FIND_PACKAGE(Perl)
IF(PERL_FOUND)
	MESSAGE(STATUS "Perl found - full test suite usable")
ELSE(PERL_FOUND)
	MESSAGE(STATUS "Perl not found - test suite only partially usable")
ENDIF(PERL_FOUND)
### insert here condition for which we do full test in case Perl is found
### or a slim version if Perl is absent

#########################################################################################
# Potentially Internal Dependencies
#
OPTION(STANDALONE "This will compile executable in static mode." OFF)

#########################################################################################
# EXTRA HEADER FILE
#
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (emmintrin.h HAVE_SSE2_H)
CHECK_INCLUDE_FILES (smmintrin.h HAVE_SSE_41_H)
IF (NOT HAVE_SSE2_H OR NOT HAVE_SSE_41_H)
#   MESSAGE(FATAL_ERROR "PFTOOLS require a compiler with both SSE 2 and 4.1 support")
ENDIF(NOT HAVE_SSE2_H OR NOT HAVE_SSE_41_H)

CHECK_INCLUDE_FILES (mm_malloc.h HAVE_MM_MALLOC_H)
CHECK_INCLUDE_FILES (alloca.h HAVE_ALLOCA_H)

SET(MMAP ON)
IF(MINGW)
  CHECK_INCLUDE_FILES (sys/mman.h HAVE_MMAN_H)
  IF(NOT HAVE_MMAN_H)
      MESSAGE(WARNING "MinGW package mman is missing, memory mapping cannot be used")
      SET(MMAP OFF)
  ENDIF(NOT HAVE_MMAN_H)
ENDIF(MINGW)

#########################################################################################
# Compiling options
#
# Use RPATH?
IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  # No rpath on Darwin. Setting it will only cause trouble.
  #SET(CMAKE_INSTALL_RPATH @executable_path/../lib)
ELSE(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  OPTION(USE_RPATH "Use -rpath when linking libraries, executables" ON)
  IF(USE_RPATH)
    SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
      SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
  ENDIF(USE_RPATH)
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")

OPTION(USE_AFFINITY "CPU core affinity will be used for threads." ON)
IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  SET(USE_AFFINITY OFF)
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")

IF(MMAP)
  OPTION(USE_MMAP "DB files wil be mapped to memory for speed." ON)
ENDIF(MMAP)

#########################################################################################
# Enabling TESTING 
IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/data")
  include(CTest)
  IF(BUILD_TESTING)
    SET(TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/Tests")
    IF (NOT EXISTS "${TESTS_DIRECTORY}" )
      FILE(MAKE_DIRECTORY "${TESTS_DIRECTORY}" )
    ENDIF (NOT EXISTS "${TESTS_DIRECTORY}" )
  ENDIF(BUILD_TESTING)
  SET(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
  SET(UTIL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/utils")
  MESSAGE(STATUS "Found data directory in ${DATA_DIR}")
  MESSAGE(STATUS "Found utils directory in ${UTIL_DIR}")
ELSE(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/data")
  MESSAGE(STATUS "No data directory found for tests")
ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/data")

#########################################################################################
# CONFIGURATION FILE
#
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#########################################################################################
# STATIC BUILD
#
IF(STANDALONE)
  ADD_DEFINITIONS("-DNO_USERNAME")
  SET(USE_RPATH OFF)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
#   SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  SET (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
  
	SET(CMAKE_SKIP_INSTALL_RPATH ON)
	SET(CMAKE_SKIP_RPATH ON)
	
ENDIF(STANDALONE)

#########################################################################################
# Compiling external tools
#
# FIND_PACKAGE(Git REQUIRED)
# INCLUDE(ExternalProject)
# SET(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/extern/)
# INCLUDE(cmake/build_libprf.cmake)

SET(LIBPRF_DIR CACHE PATH "Path to the lib PRF")
FIND_PATH(PRF_LIBRARY_DIR NAMES cmake/options.cmake
          PATHS /home/tschuepb/Projects/libprf
                ${CMAKE_SOURCE_DIR}/../libprf
                ${LIBPRF_DIR}
          DOC "option file for libPRF"
          NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH)
IF(NOT PRF_LIBRARY_DIR)
	MESSAGE(FATAL_ERROR "Unable to locate libprf source directory")
ELSE(NOT PRF_LIBRARY_DIR)
	MESSAGE(STATUS "libprf found in ${PRF_LIBRARY_DIR}")
ENDIF(NOT PRF_LIBRARY_DIR)

INCLUDE(${PRF_LIBRARY_DIR}/cmake/options.cmake)

SET(BUILD_SHARED_LIBS ON)
SET(BUILD_STATIC_LIBS OFF)
SET(BUILD_TOOLS OFF)
SET(PRF_USE_AFFINITY ${USE_AFFINITY})
SET(PRF_USE_MMAP ${USE_MMAP})
SET(PRF_USE_32BIT_INTEGER OFF)

IF (BUILD_STATIC_LIBS)
 	SET(BUILD_STATIC_DEPENDENCIES OFF)
ENDIF(BUILD_STATIC_LIBS)

SET(TOOLS_STATIC_BUILD OFF)

SET(PRF_INPUT_FASTA ON)
SET(PRF_INPUT_HDF5 ON)
SET(PRF_INPUT_PBBAM ON)

SET(PRF_CORE_STD ON)
SET(PRF_CORE_PCRE ON)
SET(PRF_CORE_DISPATCH ON)
SET(PRF_CORE_REPEAT ON)
SET(PRF_CORE_EXT_PROFILE OFF)
SET(PRF_CORE_ZONE OFF)
SET(PRF_CORE_FPGA OFF)
SET(PRF_WRAPPER_JAVA OFF)

SET(PRF_OUTPUT_PDF OFF)
SET(PRF_OUTPUT_GRAPHICS OFF)
SET(PRF_OUTPUT_DATA ON)

SET(PRF_OUTPUT_FORMAT_FASTA ON)
SET(PRF_OUTPUT_FORMAT_ONELINE ON)
SET(PRF_OUTPUT_FORMAT_XPSA ON)
SET(PRF_OUTPUT_FORMAT_TSV ON)
SET(PRF_OUTPUT_FORMAT_INTERPRO OFF)
SET(PRF_OUTPUT_FORMAT_FASEARCH OFF)
SET(PRF_OUTPUT_FORMAT_INCMATCH OFF)
SET(PRF_OUTPUT_FORMAT_PFSCAN OFF)
SET(PRF_OUTPUT_FORMAT_PSMAKER OFF)
SET(PRF_OUTPUT_FORMAT_SIMPLE OFF)
SET(PRF_OUTPUT_FORMAT_CLASSIFICATION OFF)
SET(PRF_OUTPUT_FORMAT_TEST ON)

ADD_SUBDIRECTORY(${PRF_LIBRARY_DIR} ${CMAKE_BINARY_DIR}/libprf)


#########################################################################################
## PATHs 
SET(SHARE_PATH "share/PFrepeat-${VERSION}" CACHE STRING "Installation share folder" )
SET(DOC_PATH "share/doc/PFrepeat-${VERSION}" CACHE STRING "Installation doc folder" )
SET(MAN_PATH "share/man/man1" CACHE STRING "Installation man folder" )

#########################################################################################
# Now treats subdirectories
#
ADD_SUBDIRECTORY(src)

#########################################################################################
## Data files
ADD_SUBDIRECTORY(data)

#########################################################################################
## Documentation 
#ADD_SUBDIRECTORY(doc)
#ADD_SUBDIRECTORY(man)

#########################################################################################
## Packaging
SET (CPACK_PACKAGE_NAME ${PROJECT_NAME})
SET (CPACK_PACKAGE_CONTACT "Thierry Schuepbach")
SET (CPACK_PACKAGE_VENDOR "SIB Swiss Institute of Bioinformatics")
SET (CPACK_PACKAGE_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")
SET (CPACK_PACKAGE_VERSION_MAJOR   ${MAJOR_VERSION})
SET (CPACK_PACKAGE_VERSION_MINOR  ${MINOR_VERSION})
SET (CPACK_PACKAGE_VERSION_PATCH  ${PATCH_VERSION})
SET (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_BUILD_TYPE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
SET (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}/${CPACK_PACKAGE_VERSION}")
SET (CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)

SET(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "${PROJECT_NAME}")
SET(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION  "Binary executables")

SET(CPACK_COMPONENT_DATAPRF_DISPLAY_NAME "Profiles files")
SET(CPACK_COMPONENT_DATAPRF_DESCRIPTION  "Profiles data files")
SET(CPACK_COMPONENT_DATAPRF_DISABLED TRUE)
SET(CPACK_COMPONENT_DATAPRF_DEPENDS Applications)
SET(CPACK_COMPONENT_DATAPRF_GROUP "data")

SET(CPACK_COMPONENT_DEV_DISPLAY_NAME "Development library")
SET(CPACK_COMPONENT_DEV_DESCRIPTION  "Development library")
SET(CPACK_COMPONENT_DEV_DISABLED TRUE)
SET(CPACK_COMPONENT_DEV_GROUP "dev")

SET (CPACK_COMPONENTS_ALL Applications DataPrf DataPrfSeq DataCalibration DataBlossum Documentation)


SET (CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/License.txt")

UNSET(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE)
IF(UNIX AND NOT WIN32)
    SET (CPACK_GENERATOR "TBZ2")
    SET (CPACK_DEB_COMPONENT_INSTALL ON)
    INSTALL(FILES ${PROJECT_SOURCE_DIR}/License.txt DESTINATION .
            PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
            COMPONENT Applications)
            
    IF(USE_PDF)
      SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libhpdf")
    ENDIF(USE_PDF)
    IF(USE_PCRE)
      SET(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS},pcre")
    ENDIF(USE_PCRE)
ELSE(UNIX AND NOT WIN32)
	#SET (CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE TRUE)
ENDIF(UNIX AND NOT WIN32)

SET(CPACK_SOURCE_GENERATOR "TBZ2")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "${RepeatDetecter_SOURCE_DIR}/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_IGNORE_FILES 
   "/CVS/"
   "/Build/*"
   "/build/*"
   "/\\\\.svn/*"
   "/doc/Article/pfsearch3_v0.8.odt"
   "/doc/poster/"
   "/data/Benchmarks/"
   ".*~"
   "/doc/Presentation/"
   "/src/C/experimental/"
   "/src/C/external/"
   "/src/C/generic/"
   "/data/.*\\\\.index"
   "/data/.*\\\\.fas"
   ".*\\\\.am"
   ".*\\\\.kdev4"
   "/\\\\.git/*"
   "/generic-cmake.config"
   "/output.c.diff" 
)
# SET(CPACK_SOURCE_PACKAGE_FILE_NAME "RepeatDecoder-${VERSION}")
SET(CPACK_SOURCE_STRIP_FILES "")
INCLUDE(CPack)


