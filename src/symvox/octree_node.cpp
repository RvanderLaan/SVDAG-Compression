
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

#include <symvox/geom_octree.hpp>

GeomOctree::Node::Node()
{
	children[0] = nullNode;
	children[1] = nullNode;
	children[2] = nullNode;
	children[3] = nullNode;
	children[4] = nullNode;
	children[5] = nullNode;
	children[6] = nullNode;
	children[7] = nullNode;
	childrenBitmask = 0;
	childrenMirroredBitmask[0] = 0;
	childrenMirroredBitmask[1] = 0;
	childrenMirroredBitmask[2] = 0;
	invariantBitmask = 0;
}

bool GeomOctree::Node::operator==(const GeomOctree::Node& other) const {
	if (childrenBitmask != other.childrenBitmask) return false;
	if (children[0] != other.children[0]) return false;
	if (children[1] != other.children[1]) return false;
	if (children[2] != other.children[2]) return false;
	if (children[3] != other.children[3]) return false;
	if (children[4] != other.children[4]) return false;
	if (children[5] != other.children[5]) return false;
	if (children[6] != other.children[6]) return false;
	if (children[7] != other.children[7]) return false;
	if (childrenMirroredBitmask[0] != other.childrenMirroredBitmask[0]) return false;
	if (childrenMirroredBitmask[1] != other.childrenMirroredBitmask[1]) return false;
	if (childrenMirroredBitmask[2] != other.childrenMirroredBitmask[2]) return false;
	//if (invariantBitmask != other.invariantBitmask) return false;
	return true;
}

bool GeomOctree::Node::operator<(const GeomOctree::Node& other) const {
	
	if (childrenBitmask < other.childrenBitmask) return true;
	else if (childrenBitmask > other.childrenBitmask) return false;
	
	if (children[0] < other.children[0]) return true;
	else if (children[0] > other.children[0]) return false;
	if (children[1] < other.children[1]) return true;
	else if (children[1] > other.children[1]) return false;
	if (children[2] < other.children[2]) return true;
	else if (children[2] > other.children[2]) return false;
	if (children[3] < other.children[3]) return true;
	else if (children[3] > other.children[3]) return false;
	if (children[4] < other.children[4]) return true;
	else if (children[4] > other.children[4]) return false;
	if (children[5] < other.children[5]) return true;
	else if (children[5] > other.children[5]) return false;
	if (children[6] < other.children[6]) return true;
	else if (children[6] > other.children[6]) return false;
	if (children[7] < other.children[7]) return true;
	else if (children[7] > other.children[7]) return false;

	if (childrenMirroredBitmask[0] < other.childrenMirroredBitmask[0]) return true;
	else if (childrenMirroredBitmask[0] > other.childrenMirroredBitmask[0]) return false;
	if (childrenMirroredBitmask[1] < other.childrenMirroredBitmask[1]) return true;
	else if (childrenMirroredBitmask[1] > other.childrenMirroredBitmask[1]) return false;
	if (childrenMirroredBitmask[2] < other.childrenMirroredBitmask[2]) return true;
	else if (childrenMirroredBitmask[2] > other.childrenMirroredBitmask[2]) return false;

	if (yuv[0] < other.yuv[0]) return true;
	else if (yuv[0] > other.yuv[0]) return false;
	if (yuv[1] < other.yuv[1]) return true;
	else if (yuv[1] > other.yuv[1]) return false;
	if (yuv[2] < other.yuv[2]) return true;
	else if (yuv[2] > other.yuv[2]) return false;

	return false;
}

GeomOctree::Node GeomOctree::Node::mirror(bool x, bool y, bool z, bool applyToChildren) const {
	
	Node n = *this;
	if (x) {
		Node n2;
		if (n.childrenBitmask & (1 << PXPYPZ)) n2.childrenBitmask |= 1U << NXPYPZ;
		if (n.childrenBitmask & (1 << PXPYNZ)) n2.childrenBitmask |= 1U << NXPYNZ;
		if (n.childrenBitmask & (1 << PXNYPZ)) n2.childrenBitmask |= 1U << NXNYPZ;
		if (n.childrenBitmask & (1 << PXNYNZ)) n2.childrenBitmask |= 1U << NXNYNZ;
		if (n.childrenBitmask & (1 << NXPYPZ)) n2.childrenBitmask |= 1U << PXPYPZ;
		if (n.childrenBitmask & (1 << NXPYNZ)) n2.childrenBitmask |= 1U << PXPYNZ;
		if (n.childrenBitmask & (1 << NXNYPZ)) n2.childrenBitmask |= 1U << PXNYPZ;
		if (n.childrenBitmask & (1 << NXNYNZ)) n2.childrenBitmask |= 1U << PXNYNZ;
		for (unsigned int i = 0; i < 3; ++i) {
			if (n.childrenMirroredBitmask[i] & (1U << PXPYPZ)) n2.childrenMirroredBitmask[i] |= 1U << NXPYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXPYNZ)) n2.childrenMirroredBitmask[i] |= 1U << NXPYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXNYPZ)) n2.childrenMirroredBitmask[i] |= 1U << NXNYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXNYNZ)) n2.childrenMirroredBitmask[i] |= 1U << NXNYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXPYPZ)) n2.childrenMirroredBitmask[i] |= 1U << PXPYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXPYNZ)) n2.childrenMirroredBitmask[i] |= 1U << PXPYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXNYPZ)) n2.childrenMirroredBitmask[i] |= 1U << PXNYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXNYNZ)) n2.childrenMirroredBitmask[i] |= 1U << PXNYNZ;
		}
		n2.children[PXPYPZ] = n.children[NXPYPZ];
		n2.children[PXPYNZ] = n.children[NXPYNZ];
		n2.children[PXNYPZ] = n.children[NXNYPZ];
		n2.children[PXNYNZ] = n.children[NXNYNZ];
		n2.children[NXPYPZ] = n.children[PXPYPZ];
		n2.children[NXPYNZ] = n.children[PXPYNZ];
		n2.children[NXNYPZ] = n.children[PXNYPZ];
		n2.children[NXNYNZ] = n.children[PXNYNZ];
		n = n2;
	}

	if (y) {
		Node n2;
		if (n.childrenBitmask & (1 << PXPYPZ)) n2.childrenBitmask |= 1U << PXNYPZ;
		if (n.childrenBitmask & (1 << PXPYNZ)) n2.childrenBitmask |= 1U << PXNYNZ;
		if (n.childrenBitmask & (1 << PXNYPZ)) n2.childrenBitmask |= 1U << PXPYPZ;
		if (n.childrenBitmask & (1 << PXNYNZ)) n2.childrenBitmask |= 1U << PXPYNZ;
		if (n.childrenBitmask & (1 << NXPYPZ)) n2.childrenBitmask |= 1U << NXNYPZ;
		if (n.childrenBitmask & (1 << NXPYNZ)) n2.childrenBitmask |= 1U << NXNYNZ;
		if (n.childrenBitmask & (1 << NXNYPZ)) n2.childrenBitmask |= 1U << NXPYPZ;
		if (n.childrenBitmask & (1 << NXNYNZ)) n2.childrenBitmask |= 1U << NXPYNZ;
		for (unsigned int i = 0; i < 3; ++i) {
			if (n.childrenMirroredBitmask[i] & (1U << PXPYPZ)) n2.childrenMirroredBitmask[i] |= 1U << PXNYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXPYNZ)) n2.childrenMirroredBitmask[i] |= 1U << PXNYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXNYPZ)) n2.childrenMirroredBitmask[i] |= 1U << PXPYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXNYNZ)) n2.childrenMirroredBitmask[i] |= 1U << PXPYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXPYPZ)) n2.childrenMirroredBitmask[i] |= 1U << NXNYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXPYNZ)) n2.childrenMirroredBitmask[i] |= 1U << NXNYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXNYPZ)) n2.childrenMirroredBitmask[i] |= 1U << NXPYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXNYNZ)) n2.childrenMirroredBitmask[i] |= 1U << NXPYNZ;
		}
		n2.children[PXPYPZ] = n.children[PXNYPZ];
		n2.children[PXPYNZ] = n.children[PXNYNZ];
		n2.children[PXNYPZ] = n.children[PXPYPZ];
		n2.children[PXNYNZ] = n.children[PXPYNZ];
		n2.children[NXPYPZ] = n.children[NXNYPZ];
		n2.children[NXPYNZ] = n.children[NXNYNZ];
		n2.children[NXNYPZ] = n.children[NXPYPZ];
		n2.children[NXNYNZ] = n.children[NXPYNZ];
		n = n2;
	}

	if (z) {
		Node n2;
		if (n.childrenBitmask & (1 << PXPYPZ)) n2.childrenBitmask |= 1U << PXPYNZ;
		if (n.childrenBitmask & (1 << PXPYNZ)) n2.childrenBitmask |= 1U << PXPYPZ;
		if (n.childrenBitmask & (1 << PXNYPZ)) n2.childrenBitmask |= 1U << PXNYNZ;
		if (n.childrenBitmask & (1 << PXNYNZ)) n2.childrenBitmask |= 1U << PXNYPZ;
		if (n.childrenBitmask & (1 << NXPYPZ)) n2.childrenBitmask |= 1U << NXPYNZ;
		if (n.childrenBitmask & (1 << NXPYNZ)) n2.childrenBitmask |= 1U << NXPYPZ;
		if (n.childrenBitmask & (1 << NXNYPZ)) n2.childrenBitmask |= 1U << NXNYNZ;
		if (n.childrenBitmask & (1 << NXNYNZ)) n2.childrenBitmask |= 1U << NXNYPZ;
		for (unsigned int i = 0; i < 3; ++i) {
			if (n.childrenMirroredBitmask[i] & (1U << PXPYPZ)) n2.childrenMirroredBitmask[i] |= 1U << PXPYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXPYNZ)) n2.childrenMirroredBitmask[i] |= 1U << PXPYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXNYPZ)) n2.childrenMirroredBitmask[i] |= 1U << PXNYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << PXNYNZ)) n2.childrenMirroredBitmask[i] |= 1U << PXNYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXPYPZ)) n2.childrenMirroredBitmask[i] |= 1U << NXPYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXPYNZ)) n2.childrenMirroredBitmask[i] |= 1U << NXPYPZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXNYPZ)) n2.childrenMirroredBitmask[i] |= 1U << NXNYNZ;
			if (n.childrenMirroredBitmask[i] & (1U << NXNYNZ)) n2.childrenMirroredBitmask[i] |= 1U << NXNYPZ;
		}
		n2.children[PXPYPZ] = n.children[PXPYNZ];
		n2.children[PXPYNZ] = n.children[PXPYPZ];
		n2.children[PXNYPZ] = n.children[PXNYNZ];
		n2.children[PXNYNZ] = n.children[PXNYPZ];
		n2.children[NXPYPZ] = n.children[NXPYNZ];
		n2.children[NXPYNZ] = n.children[NXPYPZ];
		n2.children[NXNYPZ] = n.children[NXNYNZ];
		n2.children[NXNYNZ] = n.children[NXNYPZ];
		n = n2;
	}

	if (applyToChildren) {
		for (unsigned int i = 0; i < 8; ++i) {
			if (n.existsChildPointer(i)) {
				if (x) n.invChildMirrrorBit(0,i);
				if (y) n.invChildMirrrorBit(1,i);
				if (z) n.invChildMirrrorBit(2,i);
			}
		}
	}

	return n;
}

GeomOctree::Node GeomOctree::Node::getCanonical(bool &x, bool &y, bool &z) const {

	Node nc = *this, ntmp;
	x = false; y = false; z = false;

	// X
	ntmp = mirror(true, false, false);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = true; y = false; z = false;
	}

	// Y
	ntmp = mirror(false, true, false);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = false; y = true; z = false;
	}

	// Z
	ntmp = mirror(false, false, true);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = false; y = false; z = true;
	}

	// XY
	ntmp = mirror(true, true, false);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = true; y = true; z = false;
	}

	// XZ
	ntmp = mirror(true, false, true);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = true; y = false; z = true;
	}

	// YZ
	ntmp = mirror(false, true, true);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = false; y = true; z = true;
	}

	// XYZ
	ntmp = mirror(true, true, true);
	if (ntmp.childrenBitmask < nc.childrenBitmask) {
		nc = ntmp;
		x = true; y = true; z = true;
	}
	return nc;
}

#define BYTETOBINARY(byte)  \
  (byte & 0x80 ? 1 : 0), \
  (byte & 0x40 ? 1 : 0), \
  (byte & 0x20 ? 1 : 0), \
  (byte & 0x10 ? 1 : 0), \
  (byte & 0x08 ? 1 : 0), \
  (byte & 0x04 ? 1 : 0), \
  (byte & 0x02 ? 1 : 0), \
  (byte & 0x01 ? 1 : 0) 

void GeomOctree::Node::print() const {
	printf("%d%d%d%d%d%d%d%d\n", BYTETOBINARY(childrenBitmask));
}