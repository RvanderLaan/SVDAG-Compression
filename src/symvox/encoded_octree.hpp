
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
#include <sl/axis_aligned_box.hpp>
#include <sl/clock.hpp>

#include <symvox/geom_octree.hpp>

class EncodedOctree {


public:
	EncodedOctree() : _rootSide(0), _levels(0), _nVoxels(0), _nNodes(0), _nLeaves(0) { _sceneBBox.to_empty(); }

	virtual ~EncodedOctree() { }

	// Encoding and storaging interface---------------------------
	virtual void encode(const GeomOctree & octree) = 0;
	virtual std::string getDescription() const = 0;
	// optional
	virtual bool load(const std::string filename) {
		printf("WARNING:EncodedOctree:load() Not implemented\n");
		return false;
	};
	virtual bool save(const std::string filename) const {
		printf("WARNING:EncodedOctree:save() Not implemented\n");
		return false;
	}
	virtual int traverse(sl::point3f p) const {
		printf("WARNING:EncodedOctree:traverse() Not implemented\n");
		return 0;
	}

	virtual int getNodeIndex(sl::point3f p, int level) const {
		printf("WARNING:EncodedOctree:getNodeIndex() Not implemented\n");
		return 0;
	}

	// Traversal interface ---------------------------------------

public:
	struct TravNode {
		TravNode() : idx(0), level(0), mX(0), mY(0), mZ(0) { center = sl::point3f(0,0,0); };
		sl::uint32_t idx;
		unsigned int level;
		sl::point3f center;
		bool mX, mY, mZ;
	};
	inline float getHalfSide(const unsigned int level = 0) const {
		return _rootSide / float(1U << (level + 1));
	}
	virtual TravNode getRootTravNode() const = 0;
	virtual bool hasChild(const TravNode &node, const int c) const = 0;
	virtual TravNode getChild(const TravNode &node, const int c, bool &mX, bool &mY, bool &mZ) const = 0;
	virtual bool isLeaf(const TravNode &node) const = 0;
	virtual int getLeafSize() const = 0;
	virtual bool hasVoxel(const TravNode &leaf, const int i, const int j, const int k) const = 0;
	
	// Getters ---------------------

public:
	inline sl::aabox3f getSceneBBox() const { return _sceneBBox; }
	inline float getRootSide() const { return _rootSide; }
	inline sl::point3f getCenter() const { return _sceneBBox.center(); }
	inline unsigned int getNLevels() const { return _levels; }
	
	inline size_t getNVoxels() const { return _nVoxels; }
	inline size_t getNLeaves() const { return _nLeaves; }
	inline size_t getNNodes() const { return _nNodes; }


protected:
	sl::aabox3f _sceneBBox;
	float _rootSide;
	unsigned int _levels;
	size_t _nVoxels, _nNodes, _nLeaves;
	sl::real_time_clock _clock;
	sl::time_duration _encodingTime;
};
