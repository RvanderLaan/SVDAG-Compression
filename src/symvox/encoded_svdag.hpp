
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


#pragma once

#include <sl/fixed_size_vector.hpp>
#include <sl/external_array.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/fixed_size_point.hpp>

#include <queue>
#include <limits>

#include <symvox/encoded_octree.hpp>


class EncodedSVDAG : public EncodedOctree {

public:
	EncodedSVDAG();

	virtual bool load(const std::string filename);
	virtual bool save(const std::string filename) const;
	virtual void encode(const GeomOctree & octree);
	virtual int getNodeIndex(sl::point3f p) const;
	virtual void * getDataPtr() const { return (void *)_data.data(); }
	virtual size_t getDataSize() const { return _data.size() * sizeof(sl::uint32_t); }
	virtual std::string getDescription() const { return "SVDAG. Encoded just like [Kampe et al.] paper."; }

	virtual TravNode getRootTravNode() const;
	virtual bool hasChild(const TravNode &node, const int c) const;
	virtual TravNode getChild(const TravNode &node, const int c, bool &mX, bool &mY, bool &mZ) const;
	virtual bool isLeaf(const TravNode &node) const;
	virtual int getLeafSize() const { return 2; }
	virtual bool hasVoxel(const TravNode &leaf, const int i, const int j, const int k) const;
	
	inline sl::uint32_t getFirstLeafPtr() { return _firstLeafPtr; }

protected:
	std::vector<sl::uint32_t> _data;
	sl::uint32_t _firstLeafPtr;
};