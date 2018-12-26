
//=============================================================
//  This file is part of the SymVox (Symmetry Voxelization) software
//  Copyright (C) 2016 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://vic.crs4.it
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//=============================================================

#include <symvox/octree.hpp>

Octree::id_t Octree::nullNode = std::numeric_limits<id_t>::max() - 1;

Octree::Octree(Scene * scene) : _scene(scene) {
	_nNodes = _nVoxels = _nLeaves = 0;
	_levels = 0;
	_rootSide = 0;
	_bbox.to_empty();
}

Octree::Octree(const Octree &other) {
	_bbox = other._bbox;
	_rootSide = other._rootSide;
	_levels = other._levels;
	_scene = other._scene;
	_nVoxels = other._nVoxels;
	_nNodes = other._nNodes;
	_nLeaves = other._nLeaves;
}

