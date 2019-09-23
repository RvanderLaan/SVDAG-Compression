
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

#include <symvox/encoded_ussvdag.hpp>
#include <symvox/util.hpp>


EncodedUSSVDAG::~EncodedUSSVDAG() {
	printf("Cleaning up EncodedUSSVDAG...\n");
	_data.clear();
	_data.shrink_to_fit();
}

bool EncodedUSSVDAG::load(const std::string filename)
{
	printf("* Loading USSVDAG '%s'... ", filename.c_str()); fflush(stdout);
	std::ifstream ifs;
	ifs.open(filename, std::ios::in | std::ios::binary);
	if (!ifs.is_open()) {
		printf("FAILED!!!\n");
		return false;
	}

	_data.clear();

	ifs.read((char *)_sceneBBox[0].to_pointer(), 12);
	ifs.read((char *)_sceneBBox[1].to_pointer(), 12);
	ifs.read((char *)&_rootSide, 4);
	ifs.read((char *)&_levels, 4);
	ifs.read((char *)&_nNodes, 4);
	ifs.read((char *)&_firstLeafPtr, 4);
	unsigned int count = (unsigned int)_data.size();
	ifs.read((char *)&count, 4);
	_data.resize(count);
	ifs.read((char *)_data.data(), count * sizeof(sl::uint32_t));

	printf("OK!\n");
	return true;
}

bool EncodedUSSVDAG::save(const std::string filename) const
{
	printf("* Saving USSVDAG '%s'... ", filename.c_str()); fflush(stdout);

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
	ofs.write((char *)&_firstLeafPtr, 4);
	unsigned int count = (unsigned int)_data.size();
	ofs.write((char *)&count, 4);
	ofs.write((char *)_data.data(), count * sizeof(sl::uint32_t));

	ofs.close();
	printf("OK!\n");
	return true;
}

void EncodedUSSVDAG::encode(const GeomOctree & octree)
{
	printf("* Encoding to USSVDAG... ");

	if (octree.getState() != GeomOctree::S_SDAG) {
		printf("FAILED! Octree is not in SDAG state\n");
		return;
	}

	_data.reserve(octree.getNNodes());
	_sceneBBox = octree.getSceneBBox();
	_rootSide = octree.getRootSide();
	_levels = octree.getLevels();
	_nVoxels = octree.getNVoxels();
	_nNodes = octree.getNNodes();

	const GeomOctree::NodeData &octData = octree.getNodeData();

	std::vector<sl::uint32_t> truePtrs;
	sl::uint32_t counter = 0;

	_clock.restart();

	for (unsigned int lev = 0; lev < _levels; ++lev) {
		for (unsigned int i = 0; i < octData[lev].size(); ++i) {
			truePtrs.push_back(counter);
			//printf("%i\t%i\n", truePtrs.size() - 1, counter);
			counter += (lev < _levels - 1) ? octData[lev][i].getNChildren() + 1 : 1;
		}
	}

	_firstLeafPtr = counter;

	unsigned int levSizeAcum = 0;

	for (unsigned int lev = 0; lev < _levels; ++lev) {
		unsigned int levSize = (unsigned int)octData[lev].size();
		levSizeAcum += levSize;
		for (unsigned int i = 0; i < levSize; ++i) {
			const GeomOctree::Node &n = octData[lev][i];
			sl::uint32_t nodeHeader =
				n.childrenMirroredBitmask[2] << 24U |
				n.childrenMirroredBitmask[1] << 16U |
				n.childrenMirroredBitmask[0] << 8U  |
				n.childrenBitmask;
				_data.push_back(nodeHeader);
			//printf("[%2d](%3d)[%i]", lev, dagEncoded.size() - 1, std::bitset<8>(bn.childrenBitmask).count());
			if (_data.size() < _firstLeafPtr) {
				counter = 0;
				for (int k = 7; k >= 0; --k) {
					if (n.children[k] != GeomOctree::nullNode) {
						_data.push_back(truePtrs[n.children[k] + levSizeAcum]);
						//printf("\t%i", dagEncoded[dagEncoded.size() - 1]);
						counter++;
					}
					//else printf("\t");
				}
#if 0
				if (check) {
					if (counter != n.getNChildren() && (lev < _levels - 1)) {
						printf("\nWRONG! bitcount != not null sons!!!\n");
						printf("Lev %u  NoID: %u  bitcount: %d  Bitmaks: \n", lev, i, n.getNChildren());
						for (int k = 0; k < 8; k++) {
							printf("  * [%i] ptr: %u\n", int((n.childrenBitmask&(1 << k))), n.children[k]);

						}
						return false;
					}
				}
#endif
			}
			else {
				//printf("\t");
				//for (int k = 7; k >= 0; --k)
				//bn.childrenBitmask & (1<<k) ? printf("1") : printf("0");
			}
			//	printf("\n");
		}
	}


	_encodingTime = _clock.elapsed();

	printf("OK! [%s]\n", sl::human_readable_duration(_encodingTime).c_str());
}

EncodedOctree::TravNode EncodedUSSVDAG::getRootTravNode() const
{
	TravNode tn;
	tn.idx = 0;
	tn.level = 0;
	return tn;
}

bool EncodedUSSVDAG::hasChild(const TravNode & node, const int c) const
{
	return (_data[node.idx] & (1U << c)) != 0;
}

EncodedOctree::TravNode EncodedUSSVDAG::getChild(const TravNode & node, const int c, bool & mX, bool & mY, bool & mZ) const
{
	TravNode tn;
	if (node.level >= _levels-1) {
		printf("EncodedUSSVDAG::getChild: WARNING! Asking for a node child, but node is in its max level\n");
		return tn;
	}
	sl::uint8_t mask = sl::uint8_t(_data[node.idx] & 0xFF);
	tn.idx = _data[node.idx + bitCount(mask >> c)];
	tn.level = node.level + 1;
	mX = (_data[node.idx] & (1U << (c + 8))) != 0;
	mY = (_data[node.idx] & (1U << (c + 16))) != 0;
	mZ = (_data[node.idx] & (1U << (c + 24))) != 0;
	return tn;
}

bool EncodedUSSVDAG::isLeaf(const TravNode & node) const
{
	return node.level == (_levels - 1);
}

bool EncodedUSSVDAG::hasVoxel(const TravNode & leaf, const int i, const int j, const int k) const
{
	int c = 4 * i + 2 * j + k;
	return hasChild(leaf, c);
}