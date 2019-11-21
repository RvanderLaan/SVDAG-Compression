
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

#include <set>

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
			simulatedEncodedSVOSize(0), memFootprint(0),
			totalLossyVoxelDifference(0), nClusteredNodes(0), nClusters(0), nEdges(0), nCrossLevelMerged(0) {}
		size_t nTotalVoxels;
		size_t nNodesSVO, nNodesDAG, nNodesSDAG;
		size_t nNodesLastLevSVO, nNodesLastLevDAG, nNodesLastLevSDAG;
		size_t simulatedEncodedSVOSize;
		size_t memFootprint;
		sl::time_duration buildSVOTime, buildDAGTime, toDAGTime, toLSVDAGTime, toSDAGTime, lHashing, lSimNodes, lClustering, crossMergeTime;
		size_t totalLossyVoxelDifference, nClusteredNodes, nClusters, nEdges, nCrossLevelMerged;
	};

public:  //////// Octree Node
	struct Node : public Octree::Node {

		sl::uint8_t childrenMirroredBitmask[3];
		sl::uint8_t invariantBitmask;
		// Indicates for every child whether it is inside or outside a surface
		sl::uint8_t outsideMask;

		// Constructor
        /**
         * @brief Node
         * @param level The level in which this Node is (initially) located. Used to initialize the level of this node's children.
         */
        Node();
//        using Octree::Node::Node;

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

        inline bool getChildOutsideBit(int childId) const { return (outsideMask & (1U << childId)) != 0; }
        inline void setChildOutsideBit(int childId) { outsideMask |= (1U << childId); }
        inline bool isInside() { return outsideMask == 0; }
        inline sl::uint8_t getOutsideMask() const { return outsideMask; }

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
    void buildDAG(unsigned int levels, unsigned int stepLevel, sl::aabox3d bbox, bool verbose = false, bool attributes = false);
	void buildSVO(unsigned int levels, sl::aabox3d bbox, bool internalCall = false, std::vector< sl::point3d > * leavesCenters = NULL, bool putMaterialIdInLeaves = false);
    void buildSVOFromPoints(std::string fileName, unsigned int levels, sl::aabox3d bbox, bool internalCall = false, std::vector< sl::point3d > * leavesCenters = NULL);
    void toDAG(bool internalCall = false);
    void toSDAG(bool internalCall = false, bool skipSymmetry = false);
	void toSDAGCanonical();

	// Extensions
    void toLossyDag(float lossyInflation, float allowedLossyDiffFactor, int includedNodeRefCount);
    unsigned int mergeAcrossAllLevels();
    void symMergeAcrossAllLevels();
    void toHiddenGeometryDAG();
		void toAttributeSVO();

    // External helpers
    void initChildLevels();
    void hiddenGeometryFloodfill();
    std::vector<std::vector<unsigned int>> sortByRefCount();
    std::vector<std::vector<unsigned int>> sortByEffectiveRefCount();

	// External tools
	bool checkIntegrity();
	//std::map<size_t, size_t> getParentLinksHistogram(int level);
	size_t getMemFootprint();
	int traverse(sl::point3f p) const;

    void DebugHashing();

private: // Internal tools
	unsigned int cleanEmptyNodes();
	void invertInvs(Node &n, int lev, bool inX, bool inY, bool inZ);

    bool compareSubtrees(unsigned int levA, unsigned int levB, Node &nA, Node &nB, std::vector<std::map<id_t, std::pair<unsigned int, id_t>>> &nodesInSubtree);
    bool compareSymSubtrees(unsigned int levA, unsigned int levB, Node &nA, Node &nB, bool sX, bool sY, bool sZ);

    unsigned int findAllSymDuplicateSubtrees();

    void diffSubtrees(unsigned int levA, unsigned int levB, const Node &nA, const Node &nB,
                 const unsigned int abortThreshold, unsigned int &accumulator, unsigned int &numLeaves);

    uint64_t computeNodeHash(const Node &node, unsigned int depth, std::vector<uint64_t> &childHashes);

    void buildMultiMap(
            unsigned int depth,
            std::vector<std::multimap<uint64_t, id_t>> &matchMaps,
            std::vector<std::vector<uint64_t>> &hashes,
            unsigned int levStart,
            unsigned int levEnd);

    inline void buildMultiMap(unsigned int depth, std::vector<std::multimap<uint64_t, id_t>> &matchMaps, std::vector<std::vector<uint64_t>> &hashes) {
        this->buildMultiMap(depth, matchMaps, hashes, 1, _levels);
    };

    uint64_t computeNodeHashBotUp(const Node &node, const unsigned int lev, const std::vector<std::vector<uint64_t>> &hashes,
                         uint8_t childMaskMask);

    void buildHiddenGeomMultiMap(std::vector<std::multimap<uint64_t, id_t>> &matchMaps,
                       std::vector<std::vector<uint64_t>> &hashes);

};
