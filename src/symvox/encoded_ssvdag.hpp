
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


class EncodedSSVDAG : public EncodedOctree {
public:
	EncodedSSVDAG();

public:
	virtual bool load(const std::string filename);
	virtual bool save(const std::string filename) const;
	virtual void encode(const GeomOctree & octree);
	virtual size_t getDataSize() const {
		return _dataInner.size() * sizeof(sl::uint16_t) + _dataLeaves.size() * sizeof(sl::uint8_t) + _levelOffsets.size() * sizeof(sl::uint32_t);
	}
	virtual std::string getDescription() const { return "SSVDAG. Symmetry-aware DAG, with compressed encoding."; }

	virtual TravNode getRootTravNode() const;
	virtual bool hasChild(const TravNode &node, const int c) const;
	virtual TravNode getChild(const TravNode &node, const int c, bool &mX, bool &mY, bool &mZ) const;
	virtual bool isLeaf(const TravNode &node) const;
	virtual int getLeafSize() const { return 4; }
	virtual bool hasVoxel(const TravNode &leaf, const int i, const int j, const int k) const;


	virtual int traverse(sl::point3f p) const;

	inline void * getInnerNodesDataPtr() const { return (void *)&_dataInner[0]; }
	inline unsigned int getInnerNodesSize() const { return (unsigned int)_dataInner.size() * sizeof(sl::uint16_t); }
	inline unsigned int getNInnerNodes() const { return (unsigned int)_dataInner.size(); }
	inline void * getLeafNodesDataPtr() const { return (void *)&_dataLeaves[0]; }
	inline unsigned int getLeafNodesSize() const { return (unsigned int)_dataLeaves.size() * sizeof(sl::uint8_t); }
	inline void * getLevelOffsetsPtr() const { return (void *)&_levelOffsets[0]; }
	inline unsigned int getNLevelOffsets() const { return (unsigned int)_levelOffsets.size(); }
	inline void expandInnerBuffer(unsigned int upTo) { _dataInner.resize(upTo); }
	inline void compressInnerBuffer(unsigned int upTo) { _dataInner.resize(upTo); }
protected:
	std::vector<sl::uint8_t> _dataLeaves;
	std::vector<sl::uint16_t> _dataInner;
	std::vector<sl::uint32_t> _levelOffsets;
	bool compute_compact_histogram_enabled_;
	bool voxel_matrix_order_enabled_;
	bool check_decode_enabled_;

public: // stats
	struct LevelStats {
		LevelStats() : nNodes(0), avgNumPtrsPerNode(0), avgSizeNode(0), nPtr16b(0), nPtr28b(0), nPtr29b(0) {}
		size_t nNodes;
		float avgNumPtrsPerNode;
		float avgSizeNode;
		size_t nPtr16b, nPtr28b, nPtr29b;
		std::vector<sl::uint32_t> compact_histogram;
	};
	std::vector<LevelStats> levelsStats;

	void set_compute_compact_histogram(bool x);
	void set_voxel_matrix_order(bool x);

protected:
	static const int sizesMask[];
	bool existsChild(const sl::uint16_t * nodePtr, int c) const;
	bool getChildMirrorBit(int axis, const sl::uint16_t * nodePtr, int c) const;
	int getChildPtrSize(const sl::uint16_t * nodePtr, int childId) const;
	sl::uint32_t getChild(const sl::uint16_t * nodePtr, int childId, bool &mX, bool &mY, bool &mZ) const;
	int getChildMask(const sl::uint16_t * nodePtr, int childId) const;


	void build_compact_histogram_in(std::vector<sl::uint32_t>& compact_histogram,
		const std::vector< std::pair< Octree::id_t, Octree::id_t > >& hist) const;

	bool check_decode(const sl::uint16_t* encoded_node,
		sl::uint32_t child_id,
		const sl::uint16_t child_offset,
		sl::uint32_t verify_mask,
		sl::uint32_t verify_addr,
		bool verify_mx,
		bool verify_my,
		bool verify_mz) const;

	void calculateChildrenOffsetsTable(std::vector<sl::uint8_t> & table) const;
};
