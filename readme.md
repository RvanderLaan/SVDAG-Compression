# Sparse Voxel DAG compression
This repo is used for experimenting with SVDAG compression techniques for my MSc thesis.
The source code of Symmetry-aware Sparse Voxel DAGs is used as a foundation.

It includes 
* A lossy compression format (LSVDAG) that clusters similar nodes together.
* A cross-level merging approach (CSVDAG) that merges nodes across all levels of the DAG.
* An attribute compression enhancement for the bit-tree representation, by encoding attributes as Gray codes instead of a conventional binary encoding.
* A work-in-progress implementation of hidden-geometry exploitation, a type of visibly lossless compression, where nodes that are identical from the outside of water-tight meshes are merged.

Most of these implementations are contained in `src/symvox/geom_octree_extension.cpp`, though the rest of the project has been adjusted here and there as well.
The viewer application (svviewer) has been extended with a GUI, diffuse lighting, hard shadows and some debugging view modes.

The accompanying MSc thesis report can be found [here](https://repository.tudelft.nl/islandora/object/uuid:83057534-111d-43bc-84f3-67a6ffe1af3b?collection=education).
For additional information or questions, I would be happy to answer those through e-mail: `rrm.remi {AT} gmail.com`.


## Instructions
Follow the instructions of the original readme for the initial setup.
Two additional optional requirements can be installed:

* For lossy compression, the clustering program [MCL](https://micans.org/mcl/) needs to be installed in its default location (`$HOME/local/bin/mcl`). This path can be adjusted in `src/symvox/cluster.hpp`.
* For the voxelization of LAS/LAZ files (lidar point clouds), the [libLAS](https://liblas.org/) package needs to be installed before building the project.

## Original readme

====================================================\
	SymVox v0.1 (Symmetry Voxelator) Software\
	(c)2016	ViC / CRS4	http://vic.crs4.it  
\====================================================
	
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
	Alberto Jaspe Villanueva [ ajaspe {AT} crs4.it ]

