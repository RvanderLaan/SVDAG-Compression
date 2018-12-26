#.rst:
# FindSL
# ------------------
#
# Try to find SpaceLib (SL) library and include path.
# Once done this will define
# SL_FOUND
# SL_INCLUDE_DIR
# SL_LIBRARIES
#
# Internally uses
# SL_LIBRARY
# SL_LIBRARY_DEBUG

SET(SL_ROOT_DIR "" CACHE PATH "SL (SpaceLib) directory")

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 ) # 64 bits
	SET( SL_WIN_LIB_DIRS
		${SL_ROOT_DIR}/build64/src/sl/Release
		${SL_ROOT_DIR}/build64/src/sl/Debug
		$ENV{PROGRAMFILES}/sl/lib)
		
	SET( SL_UNIX_LIB_DIRS
		/usr/lib64
		/usr/local/lib64
		/sw/lib
		/opt/local/lib)
ELSE()  # 32 bits
	SET( SL_WIN_LIB_DIRS
		${SL_ROOT_DIR}/build32/src/sl/Release
		${SL_ROOT_DIR}/build32/src/sl/Debug
		${SL_ROOT_DIR}/build/src/sl/Release
		${SL_ROOT_DIR}/build/src/sl/Debug
		$ENV{PROGRAMFILES}/sl/lib)
		
	SET( SL_UNIX_LIB_DIRS
		/usr/lib
		/usr/local/lib
		/sw/lib
		/opt/local/lib)
ENDIF()

FIND_PATH( SL_INCLUDE_DIR sl/any.hpp
	/usr/include
	/usr/local/include
	/sw/include
	/opt/local/include
	${SL_ROOT_DIR}/include
	${SL_ROOT_DIR}/src
	DOC "The directory where sl/any.hpp resides")

FIND_LIBRARY( SL_LIBRARY
    NAMES 
	libsl.a	# unix-based
	sl.lib	# windows
    PATHS
	${SL_WIN_LIB_DIRS}
	${SL_UNIX_LIB_DIRS}
	DOC "The SL library")

FIND_LIBRARY( SL_LIBRARY_DEBUG
    NAMES 
	libsl_d.a	# unix-based
	sl_d.lib	# windows
    PATHS
	${SL_WIN_LIB_DIRS}
	${SL_UNIX_LIB_DIRS}
	DOC "The SL library")

SET(SL_FOUND "NO")

IF (SL_INCLUDE_DIR AND SL_LIBRARY)
	SET(SL_FOUND "YES")
	IF(SL_LIBRARY_DEBUG)
		MESSAGE( "SL Library: Found! (both Release and Debug)" )
		SET(SL_LIBRARIES optimized ${SL_LIBRARY} debug ${SL_LIBRARY_DEBUG})
	ELSE()
		MESSAGE( "SL Library: Found! (only Release)" )
		SET(SL_LIBRARIES ${SL_LIBRARY})
	ENDIF()
ELSE (SL_INCLUDE_DIR AND SL_LIBRARY)
	MESSAGE( WARNING "SL Library: Not found!" )
ENDIF (SL_INCLUDE_DIR AND SL_LIBRARY)


INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SL DEFAULT_MSG SL_LIBRARY)

IF (SL_FOUND AND CMAKE_COMPILER_IS_GNUCXX)
    ##define _ISOC9X_SOURCE 1
    ADD_DEFINITIONS(-D_ISOC9X_SOURCE) 
    ADD_DEFINITIONS(-D_ISOC99_SOURCE)
    ADD_DEFINITIONS(-D__USE_ISOC9X)
    ADD_DEFINITIONS(-D__USE_ISOC99)

    # Compiler flags
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"  )
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC"  )
ENDIF (SL_FOUND AND CMAKE_COMPILER_IS_GNUCXX)
