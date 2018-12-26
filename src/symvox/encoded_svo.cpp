
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


#include <iterator>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <bitset>
#include <stack>

#include <symvox/encoded_svo.hpp>
#include <symvox/util.hpp>

EncodedSVO::EncodedSVO() : EncodedOctree()
{
	_data.clear();
}

bool EncodedSVO::load(const std::string filename)
{

	printf("* Loading SVO '%s'... ", filename.c_str()); fflush(stdout);

	std::ifstream ifs;
	ifs.open(filename, std::ios::in | std::ios::binary);
	if (!ifs.is_open()) {
		printf(" FAILED!!!\n");
		return false;
	}

	_data.clear();

	ifs.read((char *)_sceneBBox[0].to_pointer(), 12);
	ifs.read((char *)_sceneBBox[1].to_pointer(), 12);
	ifs.read((char *)&_rootSide, 4);
	ifs.read((char *)&_levels, 4);
	ifs.read((char *)&_nNodes, 4);
	unsigned int count = (unsigned int)_data.size();
	ifs.read((char *)&count, 4);
	_data.resize(count);
	ifs.read((char *)_data.data(), count * sizeof(Node));

	printf("OK!\n");

	return true;
}

bool EncodedSVO::save(const std::string filename) const
{
	printf("* Saving SVO '%s'... ", filename.c_str()); fflush(stdout);

	std::ofstream ofs;
	ofs.open(filename, std::ios::out | std::ios::binary);
	if (!ofs.is_open()) {
		printf("FAILED!!!\n");
		return false;
	}

	ofs.write((char *)_sceneBBox[0].to_pointer(), 12);
	ofs.write((char *)_sceneBBox[1].to_pointer(), 12);
	ofs.write((char *)&_rootSide, 4);
	ofs.write((char *)&_levels, 4);
	ofs.write((char *)&_nNodes, 4);
	unsigned int count = (unsigned int)_data.size();
	ofs.write((char *)&count, 4);
	ofs.write((char *)_data.data(), count * sizeof(Node));

	ofs.close();

	printf("OK!\n");

	return true;
}

void EncodedSVO::encode(const GeomOctree & octree) {
	
	printf("* Encoding SVO...      "); fflush(stdout);

	if (octree.getState() != GeomOctree::S_SVO) {
		printf("FAILED! Octree is not in SVO state\n");
		return;
	}

	struct EncodeTravNode {
		EncodeTravNode(sl::uint8_t level, GeomOctree::id_t buildNodeId, sl::uint32_t encodedNodeId) :
			level(level), buildNodeId(buildNodeId), encodedNodeId(encodedNodeId) {}
		sl::uint8_t level;
		GeomOctree::id_t buildNodeId;
		sl::uint32_t encodedNodeId;
	};
	
	typedef std::stack<EncodeTravNode> EncodeNodeStack;

	_data.reserve(octree.getNNodes());
	_sceneBBox = octree.getSceneBBox();
	_rootSide = octree.getRootSide();
	_levels = octree.getLevels();
	_nVoxels = octree.getNVoxels();
	_nNodes = octree.getNNodes();

	const GeomOctree::NodeData octData = octree.getNodeData();

	EncodeNodeStack stack;

	stack.push(EncodeTravNode(0, 0, 0));
	_data.push_back(Node());

	size_t nProcessed = 0;

	_clock.restart();

	while (!stack.empty()) {

		// Logging
		++nProcessed;
		//if (nProcessed % 10000 == 0) {
		//	printf("\b\b\b\b\b%3d %%", (100 * nProcessed) / (_nNodes - _nLeaves));
		//	fflush(stdout);
		//}
		// --------

		EncodeTravNode travNode = stack.top(); stack.pop();

		const GeomOctree::Node &node = octData[travNode.level][travNode.buildNodeId];
		_data[travNode.encodedNodeId].ptr = (travNode.level < (_levels - 1)) ? sl::uint32_t(_data.size()) : 0;
		for (int c = 7; c >= 0; --c) {
			if (node.children[c] != GeomOctree::nullNode && travNode.level < (_levels - 1)) {
				stack.push(EncodeTravNode(travNode.level + 1, node.children[c], sl::uint32_t(_data.size())));
				_data.push_back(Node());
			}
		}
		_data[travNode.encodedNodeId].data = (node.getNChildren() << 16) | node.childrenBitmask;
	}


	printf("\b\b\b\b\b OK! [%s]    \n", sl::human_readable_duration(_clock.elapsed()).c_str());
	printf("\t- Compact SVO elements: %zu (%s)\n", _data.size(), sl::human_readable_size(getDataSize()).c_str());
}

void EncodedSVO::compactData()
{
	if (_data.empty()) return;
	_compactData.clear();

	for (unsigned int i = 0; i < _data.size() / 3; i++)
	{
		_compactData.push_back(_data[3 * i + 0].ptr);
		_compactData.push_back(_data[3 * i + 1].ptr);
		_compactData.push_back(_data[3 * i + 2].ptr);

		sl::uint32_t masks =
			(_data[3 * i + 0].data & 0x000000FF) << 0 |
			(_data[3 * i + 1].data & 0x000000FF) << 8 |
			(_data[3 * i + 2].data & 0x000000FF) << 16;

		_compactData.push_back(masks);
	}

	if (_data.size() % 3 == 1) {
		_compactData.push_back(_data[_data.size() - 1].ptr);
		_compactData.push_back(0);
		_compactData.push_back(0);
		sl::uint32_t masks = (_data[_data.size() - 1].data & 0x000000FF);
		_compactData.push_back(masks);
	}
	else if (_data.size() % 3 == 2) {
		_compactData.push_back(_data[_data.size() - 2].ptr);
		_compactData.push_back(_data[_data.size() - 1].ptr);
		_compactData.push_back(0);
		sl::uint32_t masks =
			(_data[_data.size() - 2].data & 0x000000FF) << 0 |
			(_data[_data.size() - 1].data & 0x000000FF) << 8;
		_compactData.push_back(masks);

	}
}

unsigned int EncodedSVO::getPointMaxLevel(sl::point3f p)
{
	if (!_sceneBBox.contains(p)) return 0;

	sl::vector3f cv[8];
	cv[GeomOctree::PXPYPZ] = sl::vector3f(+1, +1, +1); // PXPYPZ
	cv[GeomOctree::PXPYNZ] = sl::vector3f(+1, +1, -1); // PXPYNZ
	cv[GeomOctree::PXNYPZ] = sl::vector3f(+1, -1, +1); // PXNYPZ
	cv[GeomOctree::PXNYNZ] = sl::vector3f(+1, -1, -1); // PXNYNZ
	cv[GeomOctree::NXPYPZ] = sl::vector3f(-1, +1, +1); // NXPYPZ
	cv[GeomOctree::NXPYNZ] = sl::vector3f(-1, +1, -1); // NXPYNZ
	cv[GeomOctree::NXNYPZ] = sl::vector3f(-1, -1, +1); // NXNYPZ
	cv[GeomOctree::NXNYNZ] = sl::vector3f(-1, -1, -1); // NXNYNZ

	bool found = false;
	unsigned int level = 0;
	Node node = _data[0];
	sl::point3f nodeCenter = _sceneBBox.center();
	while (1) {
		sl::uint8_t gt[3];
		gt[0] = (p[0] > nodeCenter[0]) ? 1 : 0;
		gt[1] = (p[1] > nodeCenter[1]) ? 1 : 0;
		gt[2] = (p[2] > nodeCenter[2]) ? 1 : 0;
		sl::uint32_t selectedChildID = 4 * gt[0] + 2 * gt[1] + gt[2];
		//printf("SELECTED: %i", selectedChildID);
		sl::uint32_t selectedChildMask = 1 << selectedChildID;
		if ((node.data & 0x000000FF) & selectedChildMask) {
			level++;
			if (level > (_levels - 1)) return level;
			std::bitset<8> bitsetChildren(sl::uint8_t((node.data & 0x000000FF) >> selectedChildID));
			node = _data[node.ptr + bitsetChildren.count() - 1];
			nodeCenter += cv[selectedChildID] * getHalfSide(level);
		}
		else return level;
	}
}


EncodedSVO::TravNode EncodedSVO::getRootTravNode() const
{
	TravNode tn;
	tn.idx = 0;
	tn.level = 0;
	return tn;
}

bool EncodedSVO::hasChild(const TravNode & node, const int c) const
{
	return (_data[node.idx].data & (1U << c)) != 0;
}

EncodedSVO::TravNode EncodedSVO::getChild(const TravNode & node, const int c, bool & mX, bool & mY, bool & mZ) const
{
	TravNode tn;
	if (node.level >= _levels-1) {
		printf("EncodedSVDAG::getChild: WARNING! Asking for a node child, but node is in its max level\n");
	}

	sl::uint8_t mask = sl::uint8_t(_data[node.idx].data & 0xFF);
	tn.idx = _data[node.idx].ptr + bitCount(mask >> c) - 1;
	tn.level = node.level + 1;
	mX = mY = mZ = false;
	return tn;
}

bool EncodedSVO::isLeaf(const TravNode & node) const
{
	return node.level == (_levels-1);
}

bool EncodedSVO::hasVoxel(const TravNode & leaf, const int i, const int j, const int k) const
{
	int c = 4 * i + 2 * j + k;
	return hasChild(leaf, c);
}



void EncodedSVO::print()
{
	printf("\n");
	for (std::size_t i = 0; i < _data.size(); i++) {
		if (_data[i].ptr != 0) {
			sl::uint32_t childrenCount = _data[i].data >> 16;
			char bitmask[9];
			bitmask[8] = '\0';
			for (int j = 0; j < 8; j++) (_data[i].data & (1 << j)) ? bitmask[j] = '1' : bitmask[j] = '0';
			printf("%zu\t%u\t%u\t%s\n", i, _data[i].ptr, childrenCount, bitmask);
		}
		else {
			sl::uint8_t nx = (_data[i].data & 0x00FF0000) >> 16;
			sl::uint8_t ny = (_data[i].data & 0x0000FF00) >> 8;
			sl::uint8_t nz = (_data[i].data & 0x000000FF);
			printf("%zu\tLEAF\t%u\t%u\t%u\n", i, nx, ny, nz);
		}
	}
}

