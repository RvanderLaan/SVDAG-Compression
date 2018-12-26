
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

#include <queue>
#include <limits>
#include <sl/fixed_size_vector.hpp>

#include <symvox/encoded_octree.hpp>

class EncodedSVO : public EncodedOctree {
private:
	struct Node {
		sl::uint32_t ptr;
		sl::uint32_t data;
	};
public:
	EncodedSVO();

	// EncodedOctree interface
	virtual bool load(const std::string filename);
	virtual bool save(const std::string filename) const;
	virtual void encode(const GeomOctree & octree);
	virtual void * getDataPtr() const { return (void *)_data.data(); }
	virtual size_t getDataSize() const { return _data.size() * sizeof(Node); }
	virtual std::string getDescription() const { return "SVO. Simple SVO encoding, fixed nodes and 32b pointers."; }
	
	virtual TravNode getRootTravNode() const;
	virtual bool hasChild(const TravNode &node, const int c) const;
	virtual TravNode getChild(const TravNode &node, const int c, bool &mX, bool &mY, bool &mZ) const;
	virtual bool isLeaf(const TravNode &node) const;
	virtual int getLeafSize() const { return 2; }
	virtual bool hasVoxel(const TravNode &leaf, const int i, const int j, const int k) const;


	// others
	inline sl::uint32_t * getCompactPointer() { return _compactData.data(); }
	inline std::size_t getCompactSize() { return sizeof(sl::uint32_t) * _compactData.size(); }

	// tools
	unsigned int getPointMaxLevel(sl::point3f p);
	void compactData();

	// debug stuff
	void print();



private:
	std::vector<Node> _data;
	std::vector<sl::uint32_t> _compactData;
};