
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

#include <sl/cstdint.hpp>
#include <sl/axis_aligned_box.hpp>
#include <sl/clock.hpp>
#include <sl/external_array.hpp>

#include <symvox/scene.hpp>
#include <symvox/octree.hpp>

class GeomOctree : public Octree {

public:
	enum State {
		S_EMPTY,
		S_SVO,
		S_DAG,
		S_SDAG
	};

	struct Stats {
		Stats() :
			nTotalVoxels(0),
			nNodesSVO(0), nNodesDAG(0), nNodesSDAG(0),
			nNodesLastLevSVO(0), nNodesLastLevDAG(0), nNodesLastLevSDAG(0),
			simulatedEncodedSVOSize(0), memFootprint(0) {}
		size_t nTotalVoxels;
		size_t nNodesSVO, nNodesDAG, nNodesSDAG;
		size_t nNodesLastLevSVO, nNodesLastLevDAG, nNodesLastLevSDAG;
		size_t simulatedEncodedSVOSize;
		size_t memFootprint;
		sl::time_duration buildSVOTime, buildDAGTime, toDAGTime, toSDAGTime;

	};

public:  //////// Octree Node
	struct Node : public Octree::Node {

		sl::uint8_t childrenMirroredBitmask[3];
		sl::uint8_t invariantBitmask;

		// Constructor
		Node();

		// Operators
		virtual bool operator==(const Node& other) const;
		virtual bool operator<(const Node& other) const;

		inline bool getInvariantBit(int axis) const { return (invariantBitmask & (1U << axis)) != 0; }
		inline void setInvariantBit(int axis) { invariantBitmask |= (1U << axis); }
		inline void unsetInvariantBit(int axis) { invariantBitmask &= ~(1U << axis); }
		inline void invInvariantBit(int axis) { invariantBitmask ^= (1U << axis); }

		inline bool getChildMirrrorBit(int axis, int childId) const { return (childrenMirroredBitmask[axis] & (1U << childId)) != 0; }
		inline void setChildMirrrorBit(int axis, int childId) { childrenMirroredBitmask[axis] |= (1U << childId); }
		inline void unsetChildMirrrorBit(int axis, int childId) { childrenMirroredBitmask[axis] &= ~(1U << childId); }
		inline void invChildMirrrorBit(int axis, int childId) { childrenMirroredBitmask[axis] ^= (1U << childId); }

		Node mirror(bool x, bool y, bool z, bool applyToChildren = true) const;
		Node getCanonical(bool &x, bool &y, bool &z) const;
		virtual void print() const;
	};

	typedef std::vector< std::vector< Node > > NodeData;

private: ///////////// Private attributes
	NodeData _data;
	State _state;
	Stats _stats;
	sl::real_time_clock _clock;

public:
	GeomOctree(Scene * scene);
	GeomOctree(const GeomOctree &other);
	
	// Gettters & Setters
	inline void setScene(Scene * scene) { _scene = scene; }
	inline const NodeData& getNodeData() const { return _data; }
	inline Stats getStats() const { return _stats; }
	inline State getState()  const { return _state; }

	inline virtual void * getData() { return (void *)&_data[0][0]; }

	// Main methods
	void buildDAG(unsigned int levels, unsigned int stepLevel, sl::aabox3d bbox, bool verbose = false);
	void buildSVO(unsigned int levels, sl::aabox3d bbox, bool internalCall = false, std::vector< sl::point3d > * leavesCenters = NULL, bool putMaterialIdInLeaves = false);
    void buildSVOFromPoints(std::string fileName, unsigned int levels, sl::aabox3d bbox, bool internalCall = false, std::vector< sl::point3d > * leavesCenters = NULL);
    void toDAG(bool internalCall = false);
	void toSDAG(bool internalCall = false);
	void toSDAGCanonical();

	// External tools
	bool checkIntegrity();
	//std::map<size_t, size_t> getParentLinksHistogram(int level);
	size_t getMemFootprint();
	int traverse(sl::point3f p) const;

private: // Internal tools
	unsigned int cleanEmptyNodes();
	void invertInvs(Node &n, int lev, bool inX, bool inY, bool inZ);

};
