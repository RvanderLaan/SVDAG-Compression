
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

#include <symvox/scene.hpp>

class Octree {

public:  //////// Definitions
	typedef sl::uint32_t id_t;
	typedef sl::fixed_size_point<3, sl::uint32_t> VoxelCoord;

	// null node definition
	static id_t nullNode;

	enum ChildrenIdx {
		PXPYPZ = 7,
		PXPYNZ = 6,
		PXNYPZ = 5,
		PXNYNZ = 4,
		NXPYPZ = 3,
		NXPYNZ = 2,
		NXNYPZ = 1,
		NXNYNZ = 0
	};

	enum NeighbourIdx {
	    PZ = 5,
	    PY = 4,
	    PX = 3,
	    NZ = 2,
	    NY = 1,
	    NX = 0
	};

public:  //////// Octree Node
	struct Node {

		// Node's attributes
		id_t children[8];
		sl::uint8_t childrenBitmask;
        /** Indicates for every child in which level it is located */
        unsigned int childLevels[8];

		// Constructor
        Node() {
			childrenBitmask = 0;
			children[0] = children[1] = children[2] = children[3] =
				children[4] = children[5] = children[6] = children[7] = nullNode;
            // By default, a child is located in the level below that of its parent
            childLevels[0] = childLevels[1] = childLevels[2] = childLevels[3] =
                childLevels[4] = childLevels[5] = childLevels[6] = childLevels[7] = -1;
		};

		// Methods
		inline bool existsChild(int childId) const { return (childrenBitmask & (1U << childId)) != 0; }
		inline bool existsChildPointer(int childId) const { return (children[childId] != nullNode); }
		inline bool hasChildren() const { return (childrenBitmask != 0); }
		inline unsigned int getNChildren() const { return (childrenBitmask * 01001001001ULL & 042104210421ULL) % 017; }

        inline void setChildBit(int childId) { childrenBitmask |= (1U << childId); }
		inline void unsetChildBit(int childId) { childrenBitmask &= ~(1U << childId); }

		virtual void print() const = 0;
	};

protected: ///////////// Private attributes
	sl::aabox3f _bbox;
	double _rootSide;
	unsigned int _levels;
	Scene * _scene;
	size_t _nVoxels, _nNodes, _nLeaves;
public:
	// Constructrs
	Octree(Scene * scene = NULL);
	Octree(const Octree &other);
	virtual ~Octree() {};

	// Gettters & Setters
	inline void setScene(Scene * scene) { _scene = scene; }
	inline size_t getNNodes() const { return _nNodes; }
	inline size_t getNVoxels() const { return _nVoxels; }
	inline sl::aabox3f getSceneBBox() const { return _bbox; }
	inline double getRootSideD() const { return _rootSide; }
	inline float getRootSide() const { return float(_rootSide); }
	inline unsigned int getLevels() const { return _levels; }

	inline void resizeSceneBbox(sl::aabox3f bbox) {
		_bbox = bbox;
		sl::vector3f sides = _bbox.half_side_lengths() * 2.0f;
		_rootSide = sl::max(sl::max(sides[0], sides[1]), sides[2]);
	}

	inline float getHalfSide(const unsigned int level) const { return float(_rootSide / double(1U << (level + 1))); }
	inline double getHalfSideD(const unsigned int level) const { return double(_rootSide) / double(1U << (level + 1)); }

	inline VoxelCoord getVoxelCoord(sl::point3f point) {
		sl::vector3f k = float(1<<_levels) * ((point - _bbox[0]) / float(_rootSide));
		return VoxelCoord(sl::uint32_t(k[0]), sl::uint32_t(k[1]), sl::uint32_t(k[2]));
	}

};
