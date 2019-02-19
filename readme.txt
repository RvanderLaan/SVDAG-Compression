This repo is used for experimenting with SVDAG compression techniques for my MSc thesis.
The source code of Symmetry-aware Sparse Voxel DAGs is used as a foundation.

Original readme:

====================================================
	SymVox v0.1 (Symmetry Voxelator) Software
	(c)2016	ViC / CRS4	http://vic.crs4.it  
====================================================
	
:: Introduction
----------------
	This software is intended to show the techniques described on the papers
	on SSVDAGs (Symmetry-aware Sparse Voxel DAGs) by Visual Computing
	Group (CRS4 - Italy). It provides two tools:
		
		- svbuilder: to construct Sparse Voxel DAGs, with and without
		symmetries or encoding, as described in the paper.
		
		- svviewer: an OpenGL viewer for the encoded octrees, with an
		GLSL based traversal.
							
:: Requirments
----------------
	CMake for construction, C++ compiler, OpenMP for parallelization.
	The compilation and use as been tested in Windows (7, 8.1 and 10)
	and Linux (Ubuntu, ArchLinux) platforms, with Visual Studio
	(2013 and 2015) and GCC 4.2+.

:: Dependencies
------------------
	SpaceLand (SL): set of C++ classes for geometric computations.
	http://vic.crs4.it/vic/download/
	
	For building the OpenGL viewer, the following libraries are also required:
		
		GLFW: cross-platform OpenGL context management
		http://www.glfw.org/
		
		GLEW: OpenGL Extension Wrangler Library
		http://glew.sourceforge.net/

:: Build Instructions
---------------------
	The building is based on CMake tool. It will detect the
	required libraries and produce the platform target files
	for compiling.
	
	mkdir build; cd build; cmake ..; make
	
	The binaries will be produced in a "bin" folder. It is recommended
	to launch the viewer from this folder, as it will search for the
	required shader "octree_dda.glsl.frag" in the path "../shaders".
	
	Usage and examples will be shown when launching
	the built tools without parameters.
	
	More info on cmake: https://cmake.org/runningcmake/

:: License
------------
	The SymVox software has been developed by CRS4
	(Centro di Ricerca, Sviluppo e Studi Superiori in Sardegna), and is
	available under the terms of the GNU General Public License v3.
	https://www.gnu.org/licenses/gpl-3.0.txt

	The development of commercial software based on SymVox and/or
	software based on SymVox whose source code you wish to keep
	private is not allowed without an authorization granted by CRS4.
	This license is governed by the Laws of Italy. Disputes shall be
	settled by Cagliari City Court.
	
:: Changes
-----------
	[31/8/2016] v0.1: First release
	
:: Contact
-----------
	Visual Computing Group (CRS4 - Pula, Italy)
	http://vic.crs4.it
	
	For technical questions, you can contact directly 
	Alberto Jaspe Villanueva [ ajaspe@crs4.it ]

