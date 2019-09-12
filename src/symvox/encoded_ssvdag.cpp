
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


#include <iostream>
#include <fstream>
#include <cstring>

#include <symvox/encoded_ssvdag.hpp>
#include <symvox/util.hpp>

const int EncodedSSVDAG::sizesMask[] = { 0,1,2,2 };

EncodedSSVDAG::EncodedSSVDAG() {
	compute_compact_histogram_enabled_ = false;
	voxel_matrix_order_enabled_ = true;
	check_decode_enabled_ = true;
}

EncodedSSVDAG::~EncodedSSVDAG() {
	printf("Cleaning up EncodedSSVDAG...\n");
    _dataLeaves.clear();
    _dataLeaves.shrink_to_fit();
    _dataInner.clear();
    _dataInner.shrink_to_fit();
    _levelOffsets.clear();
    _levelOffsets.shrink_to_fit();
}

bool EncodedSSVDAG::load(const std::string filename)
{
	printf("* Loading SSVDAG '%s'... ", filename.c_str()); fflush(stdout);
	std::ifstream ifs;
	ifs.open(filename, std::ios::in | std::ios::binary);
	if (!ifs.is_open()) {
		printf("FAILED!!!\n");
		return false;
	}

	_dataInner.clear();

	ifs.read((char *)_sceneBBox[0].to_pointer(), 12);
	ifs.read((char *)_sceneBBox[1].to_pointer(), 12);
	ifs.read((char *)&_rootSide, 4);
	ifs.read((char *)&_levels, 4);
	ifs.read((char *)&_nNodes, 4);
	unsigned int count;
	// inner nodes
	count = (unsigned int)_dataInner.size();
	ifs.read((char *)&count, 4);
	_dataInner.resize(count);
	ifs.read((char *)_dataInner.data(), count * sizeof(sl::uint16_t));
	// leaves
	count = (unsigned int)_dataLeaves.size();
	ifs.read((char *)&count, 4);
	_dataLeaves.resize(count);
	ifs.read((char *)_dataLeaves.data(), count * sizeof(sl::uint8_t));
	// level offset table
	count = (unsigned int)_levelOffsets.size();
	ifs.read((char *)&count, 4);
	_levelOffsets.resize(count);
	ifs.read((char *)_levelOffsets.data(), count * sizeof(sl::uint32_t));

	printf("OK!\n");
	return true;
}

bool EncodedSSVDAG::save(const std::string filename) const
{
	printf("* Saving SSVDAG '%s'... ", filename.c_str()); fflush(stdout);

	std::ofstream ofs;
	ofs.open(filename, std::ios::out | std::ios::binary);
	if (!ofs.is_open()) {
		printf("FAILED!\n");
		return false;
	}

	ofs.write((char *)_sceneBBox[0].to_pointer(), 12);
	ofs.write((char *)_sceneBBox[1].to_pointer(), 12);
	ofs.write((char *)&_rootSide, 4);
	ofs.write((char *)&_levels, 4);
	ofs.write((char *)&_nNodes, 4);
	unsigned int count;
	// inner nodes
	count = (unsigned int)_dataInner.size();
	ofs.write((char *)&count, 4);
	ofs.write((char *)_dataInner.data(), count * sizeof(sl::uint16_t));
	// leaves
	count = (unsigned int)_dataLeaves.size();
	ofs.write((char *)&count, 4);
	ofs.write((char *)_dataLeaves.data(), count * sizeof(sl::uint8_t));
	// level offset table
	count = (unsigned int)_levelOffsets.size();
	ofs.write((char *)&count, 4);
	ofs.write((char *)_levelOffsets.data(), count * sizeof(sl::uint32_t));

	ofs.close();
	printf("OK!\n");
	return true;
}

static std::pair<sl::uint8_t, sl::uint8_t> xyz_to_byte_voxel(sl::uint8_t x, sl::uint8_t y, sl::uint8_t z) {
	// Convert xyz matrix voxel position to  byte/index of the hierarchical pos
	sl::uint8_t byte_x = x / 2;
	sl::uint8_t byte_y = y / 2;
	sl::uint8_t byte_z = z / 2;
	sl::uint8_t bit_x = x % 2;
	sl::uint8_t bit_y = y % 2;
	sl::uint8_t bit_z = z % 2;
#if 0
	sl::uint8_t byte_id = byte_x * 4 + byte_y * 2 + byte_z;
	sl::uint8_t bit_id = bit_x * 4 + bit_y * 2 + bit_z;
#else
	sl::uint8_t byte_id = byte_x + byte_y * 2 + byte_z * 4;
	sl::uint8_t bit_id = bit_x + bit_y * 2 + bit_z * 4;
#endif
	return std::make_pair(byte_id, bit_id);
}

bool EncodedSSVDAG::check_decode(const sl::uint16_t* encoded_node,
	sl::uint32_t child_id,
	const sl::uint16_t child_offset,
	sl::uint32_t verify_mask,
	sl::uint32_t verify_addr,
	bool verify_mx,
	bool verify_my,
	bool verify_mz) const {
	// CHECK MASK
	int mask = getChildMask(encoded_node, child_id);
	if (mask != verify_mask) {
		std::cerr << "EncodedSSVDAG::check_decode(): WRONG CHILD MASK " << mask << " vs " << verify_mask << std::endl;
		return false;
	}
	if (mask == 0) {
		std::cerr << "EncodedSSVDAG::check_decode(): NULL CHILD MASK" << std::endl;
		return false;
	}

	// CHECK SIMMETRY
	sl::uint16_t ptr_1 = encoded_node[child_offset];
	sl::uint32_t result = 0;
	bool mx = ((ptr_1 >> 13) & 1) == 1;
	bool my = ((ptr_1 >> 14) & 1) == 1;
	bool mz = ((ptr_1 >> 15) & 1) == 1;

	if (mx != verify_mx) { std::cerr << "EncodedSSVDAG::check_decode(): WRONG MX" << std::endl; return false; }
	if (my != verify_my) { std::cerr << "EncodedSSVDAG::check_decode(): WRONG MY" << std::endl; return false; }
	if (mz != verify_mz) { std::cerr << "EncodedSSVDAG::check_decode(): WRONG MZ" << std::endl; return false; }

	// CHECK PTR
	ptr_1 &= ~(7U << 13);

	if (mask == 1) {
		if (ptr_1 != verify_addr) {
			std::cerr << "EncodedSSVDAG::check_decode(): WRONG 1US PTR LO " << ptr_1 << " !=  " << verify_addr << std::endl;
			return false;
		}
		else {
			result = ptr_1;
		}
	}
	else if (mask > 1) {
		sl::uint16_t ptr_2 = encoded_node[child_offset + 1];
		result = (sl::uint32_t(ptr_1) << 16) | ptr_2;

		sl::uint32_t x_bit = mask & 1;
		result |= (x_bit << 29);
		if (result != verify_addr) {
			std::cerr << "EncodedSSVDAG::check_decode(): WRONG 2US PTR " << result << " !=  " << verify_addr << " x bit " << x_bit << std::endl;
			return false;
		}
	}
	//  std::cerr << "test passed, c-addr " << result << std::endl;
	return true;
}

void EncodedSSVDAG::encode(const GeomOctree & octree) {
	// ENCCODING
	// ---------
	//	leaves:	raw 4^3 voxels
	//          NO MORE! -> adressing ChildIdx[maxlev-1] * 8 + ChildIdx[maxlev]
	//
	//	nodes:	variable size
	//			header (2 bytes):	[0..1] childrenMask:
	//											00 : no child
	//											01 : 16b pointer
	//											1X : 32b pointer
	//											1X : 32b pointer
	//											(X is the 30th bit of 32bits pointers)
	//			body (2..32 bytes): children pointers
	//						if(16bits) -> [][][][.... 2^13 pointer ]
	//									  ^ ^ ^--- bits symmetry (x,y,z)
	//
	//											   __bit 30th -> bit X in header
	//						if(32bits) -> [][][][ []... 2^13 pointer ]
	//									  ^ ^ ^--- bits symmetry (x,y,z)
	// NO PADDING

	printf("* Encoding to SSVDAG... ");

	if (octree.getState() != GeomOctree::S_SDAG) {
		printf("FAILED! Octree is not in SDAG state\n");
		return;
	}

	_sceneBBox = octree.getSceneBBox();
	_rootSide = octree.getRootSide();
	_levels = octree.getLevels();
	_nVoxels = octree.getNVoxels();
	_nNodes = octree.getNNodes();

	levelsStats.resize(_levels);

	const GeomOctree::NodeData &octData = octree.getNodeData();

	std::vector< std::vector<sl::uint16_t> > encInners(_levels - 2);

	// histogram vector with pairs of idx, number of pointers to the nodes
	std::vector< std::pair< GeomOctree::id_t, GeomOctree::id_t > > hist;

	// indirection vector to map between one level indices (position in the vector)
	// and their children offsets
	std::vector<GeomOctree::id_t> indirection;

	// init for stats
	//SDAG2LevelsStats.resize(_levels - 1);

	_clock.restart();

	for (int lev = _levels - 2; lev >= 0; --lev) {

		LevelStats ls;

		if (octData[lev].size() > (1U << 30)) {
			printf("EncodedSSVDAG:encode(): ERROR!!! Level %i size %zu.\nToo big for encoding in 30 bits!\n", lev, octData[lev].size());
			exit(1);
		}

		// clean and fit the refs vector
		hist.resize(octData[lev].size());
		hist.shrink_to_fit();

		// refs initialization to ordered indices and zeros
		for (GeomOctree::id_t i = 0; i < octData[lev].size(); ++i) {
			hist[i].first = i; // idx
			hist[i].second = 0; // num of refs
		}

		if (lev > 0) { // Don't do it for the root node
			// count num of references in the superior level
			for (GeomOctree::Node n : octData[lev - 1])
				for (int c = 0; c < 8; ++c)
					if (n.existsChildPointer(c)) hist[n.children[c]].second++;

			// sort by number of references (c++11's lambdas r00lez) ;D
			std::sort(hist.begin(), hist.end(),
				[](std::pair< GeomOctree::id_t, GeomOctree::id_t > a, std::pair< GeomOctree::id_t, GeomOctree::id_t> b) { return a.second > b.second; }
			);
		}

		std::vector<GeomOctree::id_t> newIndirection(octData[lev].size());

		if (lev == _levels - 2) { // LEAVES (two last levels, 4^3)

		    // Why last 2 levels instead of 1?
		    // - such a level (2^3) of granularity leads to a high structure overhead.

			_dataLeaves.resize(hist.size() * 8);
			sl::uint8_t buf[8];
			sl::uint8_t new_buf[8];
			for (GeomOctree::id_t i = 0; i < hist.size(); ++i) {
				newIndirection[hist[i].first] = i;
				const GeomOctree::Node &n = octData[lev][hist[i].first];
				if (!voxel_matrix_order_enabled_) {
					for (int c = 7; c >= 0; --c) {
						if (n.existsChildPointer(c)) {
							const GeomOctree::Node &cn = octData[lev + 1][n.children[c]];
							GeomOctree::Node cnm = cn.mirror(n.getChildMirrrorBit(0, c),
								n.getChildMirrrorBit(1, c),
								n.getChildMirrrorBit(2, c));
							_dataLeaves[i * 8 + c] = cnm.childrenBitmask;
						}
						else _dataLeaves[i * 8 + c] = 0;
					}
				}
				else {
					// Fetch the 8 leaf bytes 
					for (int c = 7; c >= 0; --c) {
						if (n.existsChildPointer(c)) {
							const GeomOctree::Node &cn = octData[lev + 1][n.children[c]];
							GeomOctree::Node cnm = cn.mirror(n.getChildMirrrorBit(0, c),
								n.getChildMirrrorBit(1, c),
								n.getChildMirrrorBit(2, c));
							buf[c] = cnm.childrenBitmask;
						}
						else {
							buf[c] = 0;
						}
					}

					memset(new_buf, 0, 8);

					// Shuffle bytes and bits to store bits in matrix row order
					// Todo: ???
					sl::uint8_t byte_pos = 0;
					sl::uint8_t bit_pos = 0;
					sl::uint8_t output_offset = 0;
					for (sl::uint8_t z = 0; z < 4; ++z) {
						for (sl::uint8_t y = 0; y < 4; ++y) {
							for (sl::uint8_t x = 0; x < 4; ++x) {
								sl::tie(byte_pos, bit_pos) = xyz_to_byte_voxel(x, y, z);

								sl::uint8_t voxel = (buf[byte_pos] >> bit_pos) & 1;

								sl::uint8_t output_byte_index = output_offset / 8;
								sl::uint8_t output_bit_index = output_offset % 8;

								new_buf[output_byte_index] |= (voxel << output_bit_index);

								++output_offset;
							}
						}
					}

					// Write to output
					for (sl::uint8_t c = 0; c < 8; ++c) {
						_dataLeaves[i * 8 + c] = new_buf[c];
					}
				}

				ls.nNodes++;
				ls.avgSizeNode += 8;
			}
		}
		else { // INNER NODES

			for (GeomOctree::id_t i = 0; i < hist.size(); ++i) {
				const GeomOctree::Node &n = octData[lev][hist[i].first];

				sl::uint16_t tmpEncNode[100];

				// tmpEncNode[0..1] <- initialize to zero
				tmpEncNode[0] = 0;

				int tmpNodeSize = 1;

				for (int c = 7; c >= 0; --c) {
					if (n.existsChild(c)) {
						sl::uint32_t child_offset = tmpNodeSize;
						sl::uint32_t verify_mask = 0;
						GeomOctree::id_t childAddress = indirection[n.children[c]];
						if (childAddress < (1U << 13)) { // < 2^13
							sl::uint16_t encChildPtr = sl::uint16_t(childAddress);
							tmpEncNode[0] |= (1U << (2 * c)); // mask 01
							if (n.getChildMirrrorBit(0, c)) encChildPtr |= 1U << 13U; //mX
							if (n.getChildMirrrorBit(1, c)) encChildPtr |= 1U << 14U; //mY
							if (n.getChildMirrrorBit(2, c)) encChildPtr |= 1U << 15U; //mZ
							tmpEncNode[tmpNodeSize] = encChildPtr;
							tmpNodeSize += 1;
							ls.nPtr16b++;
							verify_mask = 1;
						}
						else if (childAddress < (1U << 30)) { // < 2^30
							sl::uint32_t encChildPtr = sl::uint32_t(childAddress);
							if (encChildPtr & (1U << 29U)) {
								tmpEncNode[0] |= (3U << (2 * c)); // mask 11
								encChildPtr &= (~(1U << 29));       // Set to 0 29 bit
								ls.nPtr29b++;
								verify_mask = 3;
							}
							else {
								tmpEncNode[0] |= (2U << (2 * c)); // mask 10
								ls.nPtr28b++;
								verify_mask = 2;
							}

							if (n.getChildMirrrorBit(0, c)) encChildPtr |= 1U << 29U;
							if (n.getChildMirrrorBit(1, c)) encChildPtr |= 1U << 30U;
							if (n.getChildMirrrorBit(2, c)) encChildPtr |= 1U << 31U;
							//							*(sl::uint32_t *)(&tmpEncNode[tmpNodeSize]) = encChildPtr;
							// Encode in the first ushort the most significant bits 
							tmpEncNode[tmpNodeSize] = encChildPtr >> 16;
							++tmpNodeSize;
							// Encode in the first ushort the least significant bits 
							tmpEncNode[tmpNodeSize] = encChildPtr & 0xFFFF;
							++tmpNodeSize;
							//							ls.nPtr32b++;
						}
						else {
							std::cerr << "ERRROR Unsupported pointer sz " << childAddress << " > 2^30 " << (1U << 30) << std::endl;
						}

						if (check_decode_enabled_ &&
							!check_decode(tmpEncNode, c, child_offset, verify_mask, childAddress,
								n.getChildMirrrorBit(0, c),
								n.getChildMirrrorBit(1, c),
								n.getChildMirrrorBit(2, c))) {
							exit(1);
						}

						ls.avgNumPtrsPerNode++;
					}
				}

				// update the new indirection table
				newIndirection[hist[i].first] = sl::uint32_t(encInners[lev].size());

				encInners[lev].insert(encInners[lev].end(), tmpEncNode, tmpEncNode + tmpNodeSize);

				ls.nNodes++;
				ls.avgSizeNode += tmpNodeSize;
			}
		}

		indirection = newIndirection;

		if (compute_compact_histogram_enabled_) {
			//std::cerr << "Compact histogram level " << lev << std::endl;
			build_compact_histogram_in(ls.compact_histogram, hist);
		}
		else {
			//std::cerr << "Disabled compact histogram" << std::endl;
		}

		ls.avgNumPtrsPerNode /= ls.nNodes;
		ls.avgSizeNode /= ls.nNodes;
		levelsStats[lev] = ls;
	}

	// create a index table per level
	_levelOffsets.resize(_levels - 2);
	_levelOffsets[0] = 0;
	for (int i = 1; i < _levelOffsets.size(); ++i)
		_levelOffsets[i] = sl::uint32_t(encInners[i - 1].size()) + _levelOffsets[i - 1];

	// join all the encoded inner levels into a the single encoded data array
	_dataInner.clear();
	for (int i = 0; i < encInners.size(); ++i)
		_dataInner.insert(_dataInner.end(), encInners[i].begin(), encInners[i].end());

	_encodingTime = _clock.elapsed();
	printf("OK! [%s]\n", sl::human_readable_duration(_encodingTime).c_str());

#if 1 // Debug output, with type of pointers for level
	for (int i = 0; i < levelsStats.size(); ++i) {
		printf("\t-Level %2i:  Ptr16:%12zu\tPtr28:%12zu\tPtr29:%12zu\n", i, levelsStats[i].nPtr16b, levelsStats[i].nPtr28b, levelsStats[i].nPtr29b);
	}
#endif
}

EncodedOctree::TravNode EncodedSSVDAG::getRootTravNode() const
{
	TravNode tn;
	tn.idx = 0;
	tn.level = 0;
	return tn;
}

bool EncodedSSVDAG::hasChild(const TravNode & node, const int c) const
{
	return ((_dataInner[node.idx] >> (2U * c)) & 3U) != 0;
}

EncodedOctree::TravNode EncodedSSVDAG::getChild(const TravNode & node, const int c, bool & mX, bool & mY, bool & mZ) const
{
	TravNode tn;
	if (node.level >= _levels - 2) {
		printf("EncodedSSVDAG::getChild: WARNING! Asking for a node child, but node is in its max level\n");
		return tn;
	}

	unsigned int offset = 1;
	for (int i = 7; i > c; --i) {
		offset += sl::min((_dataInner[node.idx] >> (2U * i)) & 3U, 2U);
	}
	tn.idx = _dataInner[node.idx + offset];
	sl::uint32_t childMask = (_dataInner[node.idx] >> (2U * c)) & 3U;
	mX = ((tn.idx >> 13U) & 1U) == 1;
	mY = ((tn.idx >> 14U) & 1U) == 1;
	mZ = ((tn.idx >> 15U) & 1U) == 1;
	tn.idx &= ~(7U << 13U);

	if (childMask > 1) {
		sl::uint16_t tmp = _dataInner[node.idx + offset + 1];
		tn.idx = (tn.idx << 16) | tmp;
		tn.idx |= (childMask & 1U) << 29U;
	}
	tn.level = node.level + 1;
	if (tn.level < _levels - 2)
		tn.idx += _levelOffsets[tn.level];
	else
		tn.idx *= 8;

	return tn;
}

bool EncodedSSVDAG::isLeaf(const TravNode & node) const
{
	return node.level == (_levels - 2);
}

bool EncodedSSVDAG::hasVoxel(const TravNode & leaf, const int i, const int j, const int k) const
{
	int c = i + 4 * j + 16 * k;
	int idxByte = leaf.idx + (c / 8);
	int idxBit = c % 8;
	return (_dataLeaves[idxByte] & (1U << idxBit)) != 0;
}

int EncodedSSVDAG::traverse(sl::point3f p) const
{
	if (!_sceneBBox.contains(p)) return -1;

	size_t nodeIdx = 0;
	unsigned int nodeLevel = 0;
	sl::point3f nodeCenter = _sceneBBox.center();
	bool mX = false, mY = false, mZ = false;

	while (nodeLevel < _levels) {
		const sl::uint16_t * nodePtr = &_dataInner[nodeIdx];
		if (mX) p[0] = 2.f * nodeCenter[0] - p[0];
		if (mY) p[1] = 2.f * nodeCenter[1] - p[1];
		if (mZ) p[2] = 2.f * nodeCenter[2] - p[2];
		int selectedChild = 4 * int(p[0] > nodeCenter[0]) + 2 * int(p[1] > nodeCenter[1]) + int(p[2] > nodeCenter[2]);
		if (!existsChild(nodePtr, selectedChild)) break;
		nodeLevel++;
		float hs = getHalfSide(nodeLevel);
		nodeCenter[0] += (p[0] > nodeCenter[0]) ? hs : -1.f * hs;
		nodeCenter[1] += (p[1] > nodeCenter[1]) ? hs : -1.f * hs;
		nodeCenter[2] += (p[2] > nodeCenter[2]) ? hs : -1.f * hs;
		sl::uint32_t childIdx = getChild(nodePtr, selectedChild, mX, mY, mZ);
		if (nodeLevel < _levels - 2) { // INNER
			nodeIdx = _levelOffsets[nodeLevel] + childIdx;
		}
		else if (nodeLevel == _levels - 2) { // 4^3 LEAVES
			size_t leafIdx = childIdx * 8;
			const sl::uint8_t * leafPtr = &_dataLeaves[leafIdx];
			if (mX) p[0] = 2.f * nodeCenter[0] - p[0];
			if (mY) p[1] = 2.f * nodeCenter[1] - p[1];
			if (mZ) p[2] = 2.f * nodeCenter[2] - p[2];
			selectedChild = 4 * int(p[0] > nodeCenter[0]) + 2 * int(p[1] > nodeCenter[1]) + int(p[2] > nodeCenter[2]);
			if (*(leafPtr + selectedChild) == 0) break;

			nodeLevel++;
			hs = getHalfSide(nodeLevel);
			nodeCenter[0] += (p[0] > nodeCenter[0]) ? hs : -1.f * hs;
			nodeCenter[1] += (p[1] > nodeCenter[1]) ? hs : -1.f * hs;
			nodeCenter[2] += (p[2] > nodeCenter[2]) ? hs : -1.f * hs;
			sl::uint8_t selectedVoxel = 4 * int(p[0] > nodeCenter[0]) + 2 * int(p[1] > nodeCenter[1]) + int(p[2] > nodeCenter[2]);
			if (_dataLeaves[leafIdx + selectedChild] & (1U << selectedVoxel)) nodeLevel++;

			break;
		}
	}

	return nodeLevel;
}


bool EncodedSSVDAG::existsChild(const sl::uint16_t * nodePtr, int c) const {
	return getChildMask(nodePtr, c) != 0;
}

bool EncodedSSVDAG::getChildMirrorBit(int axis, const sl::uint16_t * nodePtr, int c) const {
	bool m[3];
	getChild(nodePtr, c, m[0], m[1], m[2]);
	return m[axis];
}

int EncodedSSVDAG::getChildPtrSize(const sl::uint16_t * nodePtr, int childId) const {
	return sizesMask[getChildMask(nodePtr, childId)];
}

int EncodedSSVDAG::getChildMask(const sl::uint16_t * nodePtr, int childId) const {
	return (nodePtr[0] >> (2 * childId)) & 3U;
}


sl::uint32_t EncodedSSVDAG::getChild(const sl::uint16_t * nodePtr, int childId, bool &mX, bool &mY, bool &mZ) const {
	int offset = 1; // init to header size
	for (int i = 7; i > childId; --i) {
		offset += getChildPtrSize(nodePtr, i);
	}
	int mask = getChildMask(nodePtr, childId);
	const sl::uint16_t * childPtr = (nodePtr + offset);

	sl::uint32_t childIdx = 0;
	if (mask == 1) {
		childIdx = sl::uint32_t(*childPtr);
		mX = (childIdx & (1U << 13U)) != 0;
		mY = (childIdx & (1U << 14U)) != 0;
		mZ = (childIdx & (1U << 15U)) != 0;
		childIdx &= ~(7U << 13U);
	}
	else if (mask == 2 || mask == 3) {
		childIdx = sl::uint32_t(*((sl::uint32_t *)childPtr));
		mX = (childIdx & (1U << 29U)) != 0;
		mY = (childIdx & (1U << 30U)) != 0;
		mZ = (childIdx & (1U << 31U)) != 0;
		childIdx &= ~(7U << 29U);
		childIdx |= ((mask & 1U) << 29U);
	}
	else {
		// shouldnt arrive here
		printf("EncodedSSVDAG::getChild(): WARNING. Asking for a inexistent child\n");
	}

	return childIdx;
}

void EncodedSSVDAG::set_compute_compact_histogram(bool x) {
	compute_compact_histogram_enabled_ = x;
}

void EncodedSSVDAG::set_voxel_matrix_order(bool x) {
	voxel_matrix_order_enabled_ = x;
}

void EncodedSSVDAG::build_compact_histogram_in(std::vector<sl::uint32_t>& compact_histogram,
	const std::vector< std::pair< GeomOctree::id_t, GeomOctree::id_t > >& hist) const {
	// Make an histogram, with % values on x and y.
	// Use 10 bins
	std::size_t N = hist.size();
	std::size_t bin_count = 10;
	std::size_t compact_histogram_samples_per_bin = std::max(N / bin_count, std::size_t(1));
	compact_histogram.resize(bin_count);
	std::size_t total_count = 0;
	std::cerr << "build_compact_histogram from " << N << " samples, to " << compact_histogram_samples_per_bin << " samples per bin" << std::endl;
	for (std::size_t ci = 0; ci < bin_count; ++ci) {
		compact_histogram[ci] = 0;
		std::size_t offset = ci*compact_histogram_samples_per_bin;
		for (std::size_t s = 0; s < compact_histogram_samples_per_bin && offset < N; ++s) {
			compact_histogram[ci] += hist[offset].second;
			++offset;
		}
		total_count += compact_histogram[ci];
	}
	if (total_count) {
		for (std::size_t ci = 0; ci < bin_count; ++ci) {
			compact_histogram[ci] = (unsigned int)(std::size_t(compact_histogram[ci]) * 100 / total_count);
		}
	}
}

void EncodedSSVDAG::calculateChildrenOffsetsTable(std::vector<sl::uint8_t> & table) const {
	table.resize(65536 * 8);
	for (unsigned int childId = 0; childId < 8; ++childId) {
		for (int i = 0; i < 65536; ++i) {
			sl::uint16_t childMask(i);
			sl::uint8_t offset = 1; // init to header size
			for (unsigned int i = 7; i > childId; i--) {
				offset += sizesMask[(childMask >> (2 * childId)) & 3U];
			}
			table[childId * 65536 + childMask] = offset;
		}
	}
}
