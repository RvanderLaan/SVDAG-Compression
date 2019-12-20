
#=============================================================
#  This file is part of the SymVox (Symmetry Voxelization) software
#  Copyright (C) 2016 by CRS4 Visual Computing Group, Pula, Italy
#
#  For more information, visit the CRS4 Visual Computing Group 
#  web pages at http://vic.crs4.it
#
#  This file may be used under the terms of the GNU General Public
#  License as published by the Free Software Foundation and appearing
#  in the file LICENSE included in the packaging of this file.
#
#  CRS4 reserves all rights not expressly granted herein.
#  
#  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
#  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
#  FOR A PARTICULAR PURPOSE.
#=============================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(SymVox)

SET(SymVox_VERSION_MAJOR 0)
SET(SymVox_VERSION_MINOR 1)

############## Required Libraries

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake_modules" )

# Required for the whole system
find_package(SL REQUIRED) 

# Optional, for building the viewer
find_package(OpenGL)
find_package(GLFW)
find_package(GLEW)

# Optional, for parsing LAS files
find_package(LIBLAS)

############### Options
IF(OpenGL_FOUND AND GLFW_FOUND AND GLEW_FOUND)
	SET(BUILD_SVVIEWER TRUE CACHE BOOL "Build SymVox GL Viewer")
ELSE()
	SET(BUILD_SVVIEWER FALSE CACHE BOOL "Build SymVox GL Viewer")
    Message("OpenGL, GLFW or GLEW libraries missing. Not building SymVoX GL Viewer")
ENDIF(OpenGL_FOUND AND GLFW_FOUND AND GLEW_FOUND)

IF(LIBLAS_FOUND)
        SET(BUILD_LASPARSER TRUE CACHE BOOL "Build SymVox LAS Parser")
        MESSAGE( "LIBLAS library: Found!" )
ELSE()
        SET(BUILD_LASPARSER FALSE CACHE BOOL "Build SymVox LAS Parser")
        Message("LIBLAS library missing. Not building SymVoX LAS Parser")
ENDIF(LIBLAS_FOUND)

############## Default target to build

if(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

############## Output to the "bin" and "lib" folders

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
	STRING( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_SOURCE_DIR}/bin )
	SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_SOURCE_DIR}/lib )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

############## Removing some warnings on Visual Studio

if(MSVC)
	ADD_DEFINITIONS(/D_CRT_SECURE_NO_WARNINGS)
	ADD_DEFINITIONS(/D_SCL_SECURE_NO_WARNINGS)
endif()

############## C++ 11 compatibility on GCC and misc. stuff

# Update to C++ 17
set(CMAKE_CXX_STANDARD 14)

if(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	SET(CMAKE_COLOR_MAKEFILE ON)
endif()

############## Adding OpenMP support

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

######## Includes for everybody ############

INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/src" )
INCLUDE_DIRECTORIES( ${SL_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${GLFW_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${GLEW_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${LIBLAS_INCLUDE_DIRS} )

#############################################################
##################### SYMVOX LIBRARY ###########################
#############################################################

ADD_LIBRARY( symvox
	src/symvox/encoded_octree.hpp
	src/symvox/encoded_svo.hpp
	src/symvox/encoded_svo.cpp
	src/symvox/encoded_svdag.hpp
	src/symvox/encoded_svdag.cpp
	src/symvox/encoded_ussvdag.hpp
	src/symvox/encoded_ussvdag.cpp
	src/symvox/encoded_ssvdag.hpp
	src/symvox/encoded_ssvdag.cpp
	src/symvox/scene.cpp
	src/symvox/scene.hpp
	src/symvox/test_triangle_box.cpp
	src/symvox/test_triangle_box.hpp
	src/symvox/octree.cpp
	src/symvox/octree.hpp
	src/symvox/geom_octree.cpp
	src/symvox/geom_octree.hpp
	src/symvox/octree_node.cpp
	src/symvox/util.hpp
	shaders/octree_dda.frag.glsl
	src/symvox/geom_octree_extension.cpp
	src/symvox/cluster.cpp src/symvox/cluster.hpp
)

TARGET_LINK_LIBRARIES(symvox ${SL_LIBRARIES} ${LIBLAS_LIBRARIES})


#############################################################
####################### SYMVOX BUILDER #########################
#############################################################

ADD_EXECUTABLE( svbuilder
	src/svbuilder/main.cpp
)
TARGET_LINK_LIBRARIES(svbuilder symvox)

#############################################################
####################### SYMVOX VIEWER ##########################
#############################################################

IF(BUILD_SVVIEWER)

	# Setup imgui
	set(imgui_DIR "${PROJECT_SOURCE_DIR}/include/imgui" )
	include_directories(${imgui_DIR})

	FILE(GLOB SRC_IMGUI
		${imgui_DIR}/imgui_impl_opengl3.cpp
		${imgui_DIR}/imgui_impl_opengl3.h
		${imgui_DIR}/imgui_impl_glfw.cpp
		${imgui_DIR}/imgui_impl_glfw.h
		${imgui_DIR}/imconfig.h
		${imgui_DIR}/imgui.cpp
		${imgui_DIR}/imgui.h
		${imgui_DIR}/imgui_draw.cpp
		${imgui_DIR}/imgui_internal.h
		${imgui_DIR}/imgui_widgets.cpp
		${imgui_DIR}/imstb_rectpack.h
		${imgui_DIR}/imstb_textedit.h
		${imgui_DIR}/imstb_truetype.h
	)

	SOURCE_GROUP("imgui" FILES ${SRC_IMGUI})

	# Add source code
	ADD_EXECUTABLE( svviewer
		src/svviewer/main.cpp
		src/svviewer/camera.hpp
		src/svviewer/camera.cpp
		src/svviewer/glslprogram.hpp
		src/svviewer/glslprogram.cpp
		src/svviewer/screenquadrenderer.hpp
		src/svviewer/screenquadrenderer.cpp
		src/svviewer/renderer.hpp
		src/svviewer/renderer_monitor.hpp
		src/svviewer/renderer_monitor.cpp
		src/svviewer/octree_dda_renderer.hpp
		src/svviewer/octree_dda_renderer.cpp
		src/svviewer/halton.hpp
		src/svviewer/halton.cpp
		${SRC_IMGUI}
	)

	if(UNIX)
		set(ADDITIONAL_LIBRARIES X11 Xxf86vm m rt pthread) # Xinerama Xcursor
	else()
		set(ADDITIONAL_LIBRARIES "")
	endif()

	TARGET_LINK_LIBRARIES( svviewer 
		symvox
		${OPENGL_LIBRARIES}
		${GLFW_LIBRARIES}
		${GLEW_LIBRARIES}
		${CMAKE_DL_LIBS}
		${ADDITIONAL_LIBRARIES}
	)

ENDIF(BUILD_SVVIEWER)
