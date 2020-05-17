
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

#include <queue>
#include <stack>
#include <set>
#include <utility>
#include <omp.h>

#if BUILD_LASPARSER
#include <liblas/liblas.hpp>
#endif

#include<sl/external_array.hpp>

#include <symvox/geom_octree.hpp>
#include <symvox/test_triangle_box.hpp>
#include <symvox/util.hpp>

GeomOctree::GeomOctree(Scene * scene) : Octree(scene) {
	_state = S_EMPTY;
}


GeomOctree::GeomOctree(const GeomOctree &other) : Octree (other) {
	_state = other._state;
	_stats = other._stats;
	_data.reserve(other._data.size());
	_data = other._data;
}

void GeomOctree::buildSVOFromPoints(std::string fileName, unsigned int levels, sl::aabox3d bbox, bool internalCall, std::vector< sl::point3d > * leavesCenters) {
    if (!internalCall) printf("* Building SVO... "); fflush(stdout);
#if BUILD_LASPARSER
	//_bbox = sl::conv_to<sl::aabox3f>::from(bbox);
	auto minD = bbox[0], maxD = bbox[1];
	sl::point3f minF = sl::point3f(minD[0], minD[1], minD[2]);
	sl::point3f maxF = sl::point3f(maxD[0], maxD[1], maxD[2]);
	_bbox = sl::aabox3f(minF, maxF);

    _levels = levels;
    sl::vector3f sides =_bbox.half_side_lengths() * 2.0f;
    _rootSide = sl::max(sl::max(sides[0], sides[1]), sides[2]);

    _data.resize(_levels);

    Node root;
    _data[0].push_back(root);

    struct QueueItem {
        QueueItem(id_t id, sl::uint8_t l, sl::point3d c) : nodeID(id), level(l), center(c) {}
        id_t nodeID;
        sl::uint8_t level;
        sl::point3d center;
    };
    sl::point3d childrenCenters[8];
    std::stack<QueueItem> queue;

    if(!internalCall) _clock.restart();

    // Prepare for reading the file
    std::ifstream ifs;
    ifs.open(fileName, std::ios::in | std::ios::binary);

    liblas::ReaderFactory f;
    liblas::Reader reader = f.CreateWithStream(ifs);

    int numPoints = reader.GetHeader().GetPointRecordsCount();
    int stepLogger = (int)round(numPoints / 10.f);

    // For every point
    unsigned int iPoint = 0;
    while (reader.ReadNextPoint()) {
        liblas::Point const& p = reader.GetPoint();

        if (!internalCall && (iPoint % stepLogger == 0)) {
            printf("%.0f%%..", round(100.f * (iPoint / (float)numPoints)));
            fflush(stdout);
        }

        // Push the root node to a queue
        queue.push(QueueItem(0, 0, bbox.center()));

        // Traverse through all nodes in the tree that intersect with the triangle, starting at the root
        while (!queue.empty()) {
            const QueueItem qi = queue.top(); queue.pop();

            // Store of the position (center) of all child nodes of the current node
            double k = getHalfSideD(qi.level + 1);
            childrenCenters[PXPYPZ] = qi.center + sl::vector3d(+k, +k, +k);
            childrenCenters[PXPYNZ] = qi.center + sl::vector3d(+k, +k, -k);
            childrenCenters[PXNYPZ] = qi.center + sl::vector3d(+k, -k, +k);
            childrenCenters[PXNYNZ] = qi.center + sl::vector3d(+k, -k, -k);
            childrenCenters[NXPYPZ] = qi.center + sl::vector3d(-k, +k, +k);
            childrenCenters[NXPYNZ] = qi.center + sl::vector3d(-k, +k, -k);
            childrenCenters[NXNYPZ] = qi.center + sl::vector3d(-k, -k, +k);
            childrenCenters[NXNYNZ] = qi.center + sl::vector3d(-k, -k, -k);

            // Traverse through all child nodes
            Node & node = _data[qi.level][qi.nodeID];
            for (int i = 7; i >= 0; --i) {
                // If the point lies inside this child...
                if ((p.GetX() >= childrenCenters[i][0] -k && p.GetX() < childrenCenters[i][0] +k)
                 && (p.GetY() >= childrenCenters[i][1] -k && p.GetY() < childrenCenters[i][1] +k)
                 && (p.GetZ() >= childrenCenters[i][2] -k && p.GetZ() < childrenCenters[i][2] +k)) {
                    // Mark the child mask
                    node.setChildBit(i);
                    // If there is no child node inserted yet, and it's not a leaf, insert a child node
                    if (!node.existsChildPointer(i) && (qi.level < (_levels - 1))) {
                        node.children[i] = (id_t)_data[qi.level + 1].size();
                        _data[qi.level + 1].emplace_back();
                        _nNodes++;
                        if (leavesCenters!=NULL && (qi.level == (_levels - 2))) leavesCenters->push_back(childrenCenters[i]);
                    }
                    // If this child is not a leaf, continue intersecting its children in a future iteration
                    if((qi.level+1U) < _levels) queue.push(QueueItem(node.children[i], qi.level + 1, childrenCenters[i]));
                }
            }
        }
        iPoint++;
    }

    cleanEmptyNodes();
    if (!internalCall) _stats.buildSVOTime = _clock.elapsed();

    // compute NVoxels
    for (id_t i = 0; i < _data[_levels - 1].size(); ++i) {
        _nVoxels += _data[_levels - 1][i].getNChildren();
    }

    _state = S_SVO;
    _stats.nNodesSVO = _nNodes;
    _stats.nNodesLastLevSVO = _data[_levels - 1].size();
    _stats.simulatedEncodedSVOSize = (_stats.nNodesSVO - _stats.nNodesLastLevSVO) * 4;
    _stats.nTotalVoxels = _nVoxels;

    if(!internalCall) printf("OK! [%s]\n",sl::human_readable_duration(_stats.buildSVOTime).c_str());
#endif
}

/*
 * This function converts an unsigned binary
 * number to reflected binary Gray code.
 *
 * The operator >> is shift right. The operator ^ is exclusive or.
 *
 * From https://en.wikipedia.org/wiki/Gray_code#Converting_to_and_from_Gray_code
 */
unsigned int BinaryToGray(unsigned int num) {
    return num ^ (num >> 1);
}

// https://stackoverflow.com/questions/1737726/how-to-perform-rgb-yuv-conversion-in-c-c
#define CLIP(X) ( (X) > 255 ? 255 : (X) < 0 ? 0 : X)

// RGB -> YCbCr
#define CRGB2Y(R, G, B) CLIP((19595 * R + 38470 * G + 7471 * B ) >> 16)
#define CRGB2Cb(R, G, B) CLIP((36962 * (B - CLIP((19595 * R + 38470 * G + 7471 * B ) >> 16) ) >> 16) + 128)
#define CRGB2Cr(R, G, B) CLIP((46727 * (R - CLIP((19595 * R + 38470 * G + 7471 * B ) >> 16) ) >> 16) + 128)

// YCbCr -> RGB
#define CYCbCr2R(Y, Cb, Cr) CLIP( Y + ( 91881 * Cr >> 16 ) - 179 )
#define CYCbCr2G(Y, Cb, Cr) CLIP( Y - (( 22544 * Cb + 46793 * Cr ) >> 16) + 135)
#define CYCbCr2B(Y, Cb, Cr) CLIP( Y + (116129 * Cb >> 16 ) - 226 )

void GeomOctree::buildSVO(unsigned int levels, sl::aabox3d bbox, bool internalCall, std::vector< sl::point3d > * leavesCenters, bool putMaterialIdInLeaves) {

    if (!internalCall) printf("* Building SVO... ");
    fflush(stdout);

    //_bbox = sl::conv_to<sl::aabox3f>::from(bbox);
    auto minD = bbox[0], maxD = bbox[1];
    sl::point3f minF = sl::point3f(minD[0], minD[1], minD[2]);
    sl::point3f maxF = sl::point3f(maxD[0], maxD[1], maxD[2]);
    _bbox = sl::aabox3f(minF, maxF);

    _levels = levels;
    sl::vector3f sides = _bbox.half_side_lengths() * 2.0f;
    _rootSide = sl::max(sl::max(sides[0], sides[1]), sides[2]);

    _data.resize(_levels);

    Node root;
    _data[0].push_back(root);

    struct QueueItem {
        QueueItem(id_t id, sl::uint8_t l, sl::point3d c) : nodeID(id), level(l), center(c) {}

        id_t nodeID;
        sl::uint8_t level;
        sl::point3d center;
    };
    sl::point3d childrenCenters[8];
    std::stack<QueueItem> queue;

		// Attribute parsing utils
		auto &materials = *_scene->getMaterials();
    float u, v, w; // barycentric coords
    sl::vector3f color; // temp color container
    sl::vector2f t0, t1, t2; // temp tex coord containers

    if (!internalCall) _clock.restart();

    int stepLogger = (int) round(_scene->getNRawTriangles() / 10.f);
    // For every triangle...
    for (std::size_t iTri = 0; iTri < _scene->getNRawTriangles(); iTri++) {
        if (!internalCall && (iTri % stepLogger == 0)) {
            printf("%.0f%%..", round(100.f * (iTri / (float) _scene->getNRawTriangles())));
            fflush(stdout);
        }
        unsigned int triMatId;
        if (putMaterialIdInLeaves) triMatId = (unsigned int) _scene->getTriangleMaterialId(iTri);

        // Push the root node to a queue
        queue.push(QueueItem(0, 0, bbox.center()));

		sl::vector3f yuvContainer;

        // Traverse through all nodes in the tree that intersect with the triangle, starting at the root
        while (!queue.empty()) {
            const QueueItem qi = queue.top();
            queue.pop();

            // Store of the position (center) of all child nodes of the current node
            double k = getHalfSideD(qi.level + 1);
            childrenCenters[PXPYPZ] = qi.center + sl::vector3d(+k, +k, +k);
            childrenCenters[PXPYNZ] = qi.center + sl::vector3d(+k, +k, -k);
            childrenCenters[PXNYPZ] = qi.center + sl::vector3d(+k, -k, +k);
            childrenCenters[PXNYNZ] = qi.center + sl::vector3d(+k, -k, -k);
            childrenCenters[NXPYPZ] = qi.center + sl::vector3d(-k, +k, +k);
            childrenCenters[NXPYNZ] = qi.center + sl::vector3d(-k, +k, -k);
            childrenCenters[NXNYPZ] = qi.center + sl::vector3d(-k, -k, +k);
            childrenCenters[NXNYNZ] = qi.center + sl::vector3d(-k, -k, -k);

            // Traverse through all child nodes
            Node &node = _data[qi.level][qi.nodeID];
            for (int i = 7; i >= 0; --i) {
                // If the triangle intersects with this child...
                if (_scene->getTrianglePtr(iTri) != NULL &&
                    testTriBox(childrenCenters[i], k, _scene->getTrianglePtr(iTri))) {
                    // Mark the child mask
                    node.setChildBit(i);
                    // If there is no child node inserted yet, and it's not a leaf, insert a child node
                    if (!node.existsChildPointer(i) && (qi.level < (_levels - 1))) {
                        node.children[i] = (id_t) _data[qi.level + 1].size();
                        _data[qi.level + 1].emplace_back();
                        _nNodes++;
                        if (leavesCenters != NULL && (qi.level == (_levels - 2)))
                            leavesCenters->push_back(childrenCenters[i]);
                    }
                    // If this child is not a leaf, continue intersecting its children in a future iteration
                    if ((qi.level + 1U) < _levels)
                        queue.push(QueueItem(node.children[i], qi.level + 1, childrenCenters[i]));
                    else {
											// If it's a leaf, do nothing unless we want attribute data here
                        if (putMaterialIdInLeaves && node.children[i] == nullNode) {
                            sl::uint8_t attr; // todo: only 1 uint8 for now (gray scale), later we can add more channels
														// sl::uint64_t attr; // like this for example

                            const auto& material = materials[triMatId];
                            std::string texName(material.texture);

                            // Check if triangle has texture
                            if (this->_scene->isTriangleTextured(iTri) && !texName.empty()) {
                                // Compute barycentric coordinates of the center of this voxel to this triangle
                                barycentric(qi.center, _scene->getTrianglePtr(iTri), u, v, w);

                                // Use that to find the texture coordinates of that point on the texture
                                _scene->getTriangleTexCoords(iTri, t0, t1, t2);
                                sl::vector2f voxTexCoords = u * t0 + v * t1 + w * t2;

								float x = voxTexCoords[0];
								// Fix for tiled textures, mod between 0 and 1
								voxTexCoords[0] = std::fmod(std::fmod(voxTexCoords[0], 1.f) + 1.f, 1.f);
								voxTexCoords[1] = std::fmod(std::fmod(voxTexCoords[1], 1.f) + 1.f, 1.f);
								
                                // Look up texture color at those coordinates
                                // TODO: Use texture LOD or bigger sample size depending on size ratio of triangle to voxel
                                _scene->getTexColor(texName, voxTexCoords, color);
//                                printf("Tri# %i: \tuv: %.2f \t%.2f \t - Attr: %.2f\n", iTri, voxTexCoords[0], voxTexCoords[1], f);
                            } else {
                                color = materials[triMatId].diffuseColor;
                            }

							// New approach: Set YUV in a separate field, can be encoded in node header (replaces padding)
							// TODO: This sets it for the entire leaf node, not its 8 individual voxels
							// meh, good enough for now. Maybe encode as 64-bit leaf nodes?
							// might work with 32 bit when only storing chroma difference from parent
							// (4 bits per voxel available, should be fine in most cases)
							//RGBToYUV(color, yuvContainer);
							int r = int(color[0] * 255.),
								g = int(color[1] * 255.),
								b = int(color[2] * 255.);
							node.yuv[0] = float(CRGB2Y(r, g, b)); // yuvContainer[0];
							node.yuv[1] = float(CRGB2Cb(r, g, b));; // yuvContainer[1];
							node.yuv[2] = float(CRGB2Cr(r, g, b));; // yuvContainer[2];

#if 0 // OLD APPROACH: Put grayscale color in the node.children (pointers)
                            // Average of RGB: Gray scale
                            float f = (color[0] + color[1] + color[2]) / 3.;
                            // Todo: try both for binary and for gray code
                            attr = sl::uint8_t(floor(f >= 1.0 ? 255 : f * 255.0));
                            //attr = BinaryToGray(attr);

                            // Todo: For more than 1 attribute, use a vector of attributes
                            // this->attributes.push_back(attr);


                            // For now, just put it in children
                            node.children[i] = (id_t) attr;
#endif
                        }
											
#if 0 // outputs a debug obj of voxels as points with their colours
                        if (putMaterialIdInLeaves) node.children[i] = (id_t) triMatId;
                        sl::color3f c;
                        _scene->getTriangleColor(iTri, c);
                        printf("v %f %f %f %f %f %f\n", childrenCenters[i][0], childrenCenters[i][1], childrenCenters[i][2], c[0], c[1], c[2]);
#endif


                    }
                }
            }
        }
    }

    cleanEmptyNodes();
    if (!internalCall) _stats.buildSVOTime = _clock.elapsed();

    // compute NVoxels
    for (id_t i = 0; i < _data[_levels - 1].size(); ++i) {
        _nVoxels += _data[_levels - 1][i].getNChildren();
    }

    _state = S_SVO;
    _stats.nNodesSVO = _nNodes;
    _stats.nNodesLastLevSVO = _data[_levels - 1].size();
    _stats.simulatedEncodedSVOSize = (_stats.nNodesSVO - _stats.nNodesLastLevSVO) * 4;
    _stats.nTotalVoxels = _nVoxels;

    if(!internalCall) printf("OK! [%s]\n",sl::human_readable_duration(_stats.buildSVOTime).c_str());

}

/**
 * @brief GeomOctree::buildDAG The main function that builds sub-SVOs, individually reduces them to DAGs and merges them
 * @param levels
 * @param stepLevel
 * @param bbox
 * @param verbose
 */
void GeomOctree::buildDAG(unsigned int levels, unsigned int stepLevel, sl::aabox3d bbox, bool verbose, bool attributes) {
	
	printf("* Building DAG [stepLevel: %i]\n", stepLevel); fflush(stdout);
	
	std::vector< sl::point3d > leavesCenters;
	
	unsigned int stepLevels = stepLevel + 1;

	_clock.restart();
	sl::time_duration timeStamp = _clock.elapsed();

	printf("\t- Building root SVO subtree... ");
	
	buildSVO(stepLevels, bbox, true, &leavesCenters);
	
	printf("OK! [%s]\n", sl::human_readable_duration(_clock.elapsed() - timeStamp).c_str()); fflush(stdout);

	double lhs = getHalfSideD(stepLevels - 1);
	sl::vector3d corners[8];
	corners[PXPYPZ] = sl::vector3d(+lhs, +lhs, +lhs); // PXPYPZ
	corners[PXPYNZ] = sl::vector3d(+lhs, +lhs, -lhs); // PXPYNZ
	corners[PXNYPZ] = sl::vector3d(+lhs, -lhs, +lhs); // PXNYPZ
	corners[PXNYNZ] = sl::vector3d(+lhs, -lhs, -lhs); // PXNYNZ
	corners[NXPYPZ] = sl::vector3d(-lhs, +lhs, +lhs); // NXPYPZ
	corners[NXPYNZ] = sl::vector3d(-lhs, +lhs, -lhs); // NXPYNZ
	corners[NXNYPZ] = sl::vector3d(-lhs, -lhs, +lhs); // NXNYPZ
	corners[NXNYNZ] = sl::vector3d(-lhs, -lhs, -lhs); // NXNYNZ

	_stats.nNodesSVO = _nNodes;
	_stats.nNodesLastLevSVO = 0;
	std::map<unsigned int, GeomOctree *> leavesOctrees;

	float nodesProcessed = 0;
	const size_t nodesToProcess = _nVoxels;
	size_t memConsumed = 0;
	size_t acumSubtreesDAGNNodes = 0;

	printf("\t- Building %zu subtress (SVO->DAG) [%i threads]... ", _nVoxels, (int)omp_get_max_threads());
	if (verbose) printf("\n");
	timeStamp = _clock.elapsed();
	float acumDAGTimeFactor = 0;
#if 1 // More memory consumtion, but ensures CPUs occupancy
#pragma omp parallel for schedule(dynamic,2)
#else
#pragma omp parallel for schedule(static)
#endif
	for (int i = 0; i < _data[_levels - 1].size(); ++i) {
		for (int j = 7; j >= 0; --j) {
			if (_data[_levels - 1][i].existsChild(j)) {
				sl::time_point threatInitTime = _clock.now();
				GeomOctree * leafOctree = new GeomOctree(_scene);
				sl::aabox3d lbbox;
				lbbox.to_empty();
				lbbox.merge(leavesCenters[i]);
				lbbox.merge(leavesCenters[i] + corners[j]);
				leafOctree->buildSVO(levels - stepLevels, lbbox, true, nullptr, attributes);

				if (attributes) {
					leafOctree->propagateYUV();
				}

//				leafOctree->hiddenGeometryFloodfill(); // todo: fix this
				size_t nNodesLeafSVO = leafOctree->getNNodes();
				size_t nNodesLeafLastLevSVO = leafOctree->_data[leafOctree->_data.size() - 1].size();
#pragma omp atomic
				_stats.nNodesSVO += nNodesLeafSVO;
#pragma omp atomic
				_stats.nNodesLastLevSVO += nNodesLeafLastLevSVO;
				size_t incVoxels = leafOctree->getNVoxels() - 1;
#pragma omp atomic
				_nVoxels += incVoxels; // root doesn't count

				leafOctree->toDAG(true);

#pragma omp atomic		
				nodesProcessed += 1.f;
				size_t incNNodes = leafOctree->getNNodes();
#pragma omp atomic		
				acumSubtreesDAGNNodes += incNNodes;
				
#pragma omp atomic
				acumDAGTimeFactor += leafOctree->getStats().toDAGTime.as_microseconds() / (float) (_clock.now()- threatInitTime).as_microseconds();
#pragma omp critical
				{
					leavesOctrees[i * 8 + j] = leafOctree;
					memConsumed += leafOctree->getMemFootprint();
					_stats.memFootprint = sl::max(_stats.memFootprint, leafOctree->getMemFootprint());
					if (verbose) {
						printf("\t\t- [%3d%%] [%3u-%u / %zu (%d)]\t%zu -> %zu\n", int(100.f * (nodesProcessed / (float)nodesToProcess)), i, j, _data[_levels - 1].size(), omp_get_thread_num(), nNodesLeafSVO, leafOctree->getNNodes());
						fflush(stdout);
					}
				}
			}
		}
	}
	acumDAGTimeFactor /= nodesProcessed;
	sl::time_duration acumDAGTime = _clock.elapsed() * acumDAGTimeFactor;

	if(verbose)
		printf("\t\t- [100%%] All subtrees done! [%s]\n", sl::human_readable_duration(_clock.elapsed() - timeStamp).c_str());
	else
		printf("OK! [%s]\n", sl::human_readable_duration(_clock.elapsed() - timeStamp).c_str());

	_stats.simulatedEncodedSVOSize = (_stats.nNodesSVO - _stats.nNodesLastLevSVO) * 4;
	_stats.nTotalVoxels = _nVoxels;

	_data.resize(levels);
	_levels = levels;

	timeStamp = _clock.elapsed();

	printf("\t- Joining subtrees... ");

	for (int i = 0; i < _data[stepLevels - 1].size(); ++i) {
		for (int j = 7; j >= 0; --j) {
			if (_data[stepLevels - 1][i].existsChild(j)) {
				GeomOctree * oct = leavesOctrees[i * 8 + j];
				_data[stepLevels - 1][i].children[j] = (id_t)_data[stepLevels].size();
				for (unsigned int k = stepLevels; k < levels; ++k) {
					if (k < (levels - 1)) { // update pointers
						id_t offset = (id_t)_data[k + 1].size();
						for (unsigned int m = 0; m < oct->_data[k - stepLevels].size(); ++m) {
							for (unsigned int n = 0; n < 8; ++n) {
								if (oct->_data[k - stepLevels][m].existsChildPointer(n))
									oct->_data[k - stepLevels][m].children[n] += offset;
							}
						}
					}
					_data[k].insert(_data[k].end(), oct->_data[k - stepLevels].begin(), oct->_data[k - stepLevels].end());
				}
				delete oct;
			}
		}
	}

	printf("OK! [%s]\n", sl::human_readable_duration(_clock.elapsed() - timeStamp).c_str());
	
	printf("\t- Last DAG pass to whole subtree... ");
	acumDAGTime += (_clock.elapsed() - timeStamp);
	timeStamp = _clock.elapsed();
	
	toDAG(true);

	acumDAGTime += (_clock.elapsed() - timeStamp);
	_stats.toDAGTime = acumDAGTime;
	printf("OK! [%s]\t(%zu -> %zu)\n", sl::human_readable_duration(_clock.elapsed() - timeStamp).c_str(), acumSubtreesDAGNNodes, _nNodes);

	_stats.buildDAGTime = _clock.elapsed();
	_stats.nNodesLastLevDAG = _data[_levels - 1].size();

	printf("\t- Finished! Total time [%s]\n", sl::human_readable_duration(_stats.buildDAGTime).c_str());
}

unsigned int GeomOctree::cleanEmptyNodes() {
	std::set<id_t> emptyNodes;
	unsigned int nDelPtrs = 0;
	for (int lev = _levels - 1; lev > 0; --lev) {
		emptyNodes.clear();
		for (id_t i = 0; i < _data[lev].size(); ++i) {
			if (!_data[lev][i].hasChildren()) emptyNodes.insert(i);
		}
		for (id_t i = 0; i < _data[lev-1].size(); ++i) {
			for (int j = 0; j < 8; ++j) {
				if (emptyNodes.find(_data[lev - 1][i].children[j]) != emptyNodes.end()) {
					_data[lev-1][i].unsetChildBit(j);
					_data[lev-1][i].children[j] = nullNode;
					nDelPtrs++;
				}
			}
		}
	}
	return nDelPtrs;
}

/**
 * @brief GeomOctree::toDAG Reduces an SVO to a DAG
 * @param iternalCall Whether the function is called internally (only adds extra print)
 */
void GeomOctree::toDAG(bool iternalCall) {

	if (_state == S_SVO || _state == S_DAG) {
		if (!iternalCall) printf("* Transforming SVO -> DAG ... "); fflush(stdout);
	} else {
		printf("ERROR! This is not a SVO!\n");
		return;
	}

	_nNodes = 1;
    /** Every index denotes the index of the first duplicate of that node in uniqueNodes. Reset for each level. */
	std::vector<id_t> correspondences;
	std::map<Node, id_t> uniqueNodesChecker;
	std::vector<Node> uniqueNodes;

	 sl::time_point ts = _clock.now();

	 if (!iternalCall)
        printf("\n");

    // For every level, starting at the leaves...
	for (unsigned int lev = _levels - 1; lev > 0; --lev) {
        // Clear the lists used to keep track of correspondences etc
		size_t oldLevSize = _data[lev].size();
		uniqueNodes.clear();
		uniqueNodes.shrink_to_fit();
		uniqueNodesChecker.clear();
		correspondences.clear();
		correspondences.resize(oldLevSize);

        // For all nodes in this level...
		for (id_t i = 0; i < _data[lev].size(); i++) {
			Node n = _data[lev][i];
			if (!n.hasChildren()) continue; // skip empty nodes
            auto k = uniqueNodesChecker.find(n); // find if duplicate

			if (k != uniqueNodesChecker.end()) { // found
                correspondences[i] = (*k).second; // store duplicate node
				auto target = &(uniqueNodes[k->second]);


				// experiment: store average color value
				// New total amount that this node is referenced
				target->numRefs += n.numRefs;

				// weight; contribution of this node to avg color, depends on nr. of dupes
				// with multiple build steps, the duplicate node can also be referenced more than once
				float w = n.numRefs / (float) target->numRefs; 
				float wInv = 1. - w;

				// Weighted average
				target->yuv[0] = target->yuv[0] * wInv + w * n.yuv[0];
				target->yuv[1] = target->yuv[1] * wInv + w * n.yuv[1];
				target->yuv[2] = target->yuv[2] * wInv + w * n.yuv[2];
			}
			else { // !found
                uniqueNodesChecker[n] = (id_t)uniqueNodes.size(); // store it as unique node
                correspondences[i] = (id_t)uniqueNodes.size(); // the correspondence is this node itself
				uniqueNodes.push_back(n);
			}
		}

		if (!iternalCall)
            printf("Reduced level %u from %lu to %lu nodes\n", lev, _data[lev].size(), uniqueNodes.size());

		_data[lev].clear();
		_data[lev].shrink_to_fit();
		uniqueNodes.shrink_to_fit();
        _data[lev] = uniqueNodes; // Replace all SVO nodes with the unique DAG nodes in this level
		_data[lev].shrink_to_fit();
		_nNodes += _data[lev].size();
		
        // Update all pointers in the level above
		for (id_t i = 0; i < _data[lev-1].size(); i++) {
			Node * bn = &_data[lev-1][i];
            // For all children...
			for (int j = 0; j < 8; j++) {
                // If this child exists...
                if (bn->existsChild(j)) {
                    // Set the child pointer to the unique node that replaced this child
                    bn->children[j] = correspondences[bn->children[j]];

//                    if (lev >= _levels - 1) {
//
//                        // test to see if recursion works. post-test: YES IT DOES. Need to set LEVELS uniform to *2 though
//                        if ((std::rand() % 100) / 100.0 > 0.5) {
//                            bn->children[j] = 0;
//                            bn->childLevels[j] = 1;
//                        }
//                    }
                }
			}
		}
	}

	_stats.toDAGTime = _clock.now() - ts;
	
	_state = S_DAG;
	_stats.nNodesDAG = _nNodes;
	_stats.nNodesLastLevDAG = _data[_levels - 1].size();

	if (!iternalCall) printf("OK! [%s]\n", sl::human_readable_duration(_stats.toDAGTime).c_str());
}


void GeomOctree::toSDAG(bool internalCall, bool skipSymmetry) {

	if (!internalCall) {
		if (_state == S_SVO) {
			printf("* Transforming SVO -> SDAG (not best way)... "); fflush(stdout);
		}
		else if (_state == S_DAG) {
			printf("* Transforming DAG -> SDAG (normal)... "); fflush(stdout);
		}
		else {
			printf("ERROR! This is not a DAG or SDAG!\n");
			return;
		}
	}

	struct MirroredNode{
		MirroredNode() : mirrorX(false), mirrorY(false), mirrorZ(false), id(0) {}
		MirroredNode(bool mX, bool mY, bool mZ, id_t id) :
			mirrorX(mX), mirrorY(mY), mirrorZ(mZ), id(id) {}
		bool mirrorX, mirrorY, mirrorZ;
		id_t id;
	};

	std::vector<MirroredNode> correspondences;
	std::map<Node, id_t> uniqueNodesChecker;
	std::vector<Node> uniqueNodes;

	_nNodes = 0;

	if (!internalCall) _clock.restart();

    // For every level, starting at the leaves...
	for (unsigned int lev = _levels - 1; lev > 0; lev--) {
        // Clear the lists used to keep track of correspondeces etc
		unsigned int oldLevSize = (unsigned int)_data[lev].size();
		uniqueNodes.clear();
		uniqueNodesChecker.clear();
		correspondences.resize(oldLevSize);

        // For each node in this level...
		for (id_t i = 0; i < _data[lev].size(); i++) {

			Node n = _data[lev][i];
			if (!n.hasChildren()) continue; // skip empty nodes

			if (skipSymmetry) { // Used when we only want pointer compression
                uniqueNodesChecker[n] = (id_t)uniqueNodes.size();
                correspondences[i] = MirroredNode(false, false, false, id_t(uniqueNodes.size()));
                uniqueNodes.push_back(n);
                continue;
			}

            // Check invariances
            if (n == n.mirror(true, false, false)) n.setInvariantBit(0);
            if (n == n.mirror(false, true, false)) n.setInvariantBit(1);
            if (n == n.mirror(false, false, true)) n.setInvariantBit(2);

            // Create mirrored versions of this node
            Node niX = n.mirror(true, false, false);
            Node niY = n.mirror(false, true, false);
            Node niZ = n.mirror(false, false, true);
            Node niXY = niX.mirror(false, true, false);
            Node niXZ = niX.mirror(false, false, true);
            Node niYZ = niY.mirror(false, false, true);
            Node niXYZ = niXY.mirror(false, false, true);

            // If a child is invariant on an axis, this removes the mirror bit there
			if (lev < (_levels - 1)) {
				invertInvs(niX, lev, true, false, false);
				invertInvs(niY, lev, false, true, false);
				invertInvs(niZ, lev, false, false, true);
				invertInvs(niXY, lev, true, true, false);
				invertInvs(niXZ, lev, true, false, true);
				invertInvs(niYZ, lev, false, true, true);
				invertInvs(niXYZ, lev, true, true, true);
			}

            // Find a correspondence of this node or any of its mirrored copies
			std::map<GeomOctree::Node, sl::uint32_t>::iterator it;
			if ((it = uniqueNodesChecker.find(n)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(false, false, false, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niX)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(true, false, false, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niY)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(false, true, false, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niZ)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(false, false, true, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niXY)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(true, true, false, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niXZ)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(true, false, true, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niYZ)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(false, true, true, (*it).second);
			}
			else if ((it = uniqueNodesChecker.find(niXYZ)) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(true, true, true, (*it).second);
			}
            else { // no correspondence found, so add as unique node
				uniqueNodesChecker[n] = (id_t)uniqueNodes.size();
				correspondences[i] = MirroredNode(false, false, false, id_t(uniqueNodes.size()));
				uniqueNodes.push_back(n);
			}
		}
        // Replace previous nodes on this level as the reduced set of symmetrically unique nodes
		_data[lev].clear();
		_data[lev] = uniqueNodes;
		_data[lev].shrink_to_fit();
		_nNodes += _data[lev].size();

        // Update all pointers in the level above
		for (id_t i = 0; i < _data[lev - 1].size(); i++) {
			Node &n = _data[lev - 1][i];
			for (id_t j = 0; j < 8; j++) {
				if (n.existsChildPointer(j)) {
                    // Set the child pointer to the unique node that replaced this child
					MirroredNode mn = correspondences[n.children[j]];
					n.children[j] = mn.id;
					if (mn.mirrorX) n.setChildMirrrorBit(0, j);
					if (mn.mirrorY) n.setChildMirrrorBit(1, j);
					if (mn.mirrorZ) n.setChildMirrrorBit(2, j);
				}
			}
		}
	}

	if (!internalCall) _stats.toSDAGTime = _clock.elapsed();
	_stats.nNodesSDAG = _nNodes;
	_state = S_SDAG;
	
#if 0 // DEGUG HOW MUCH MIRRORED NODES THERE ARE
	for (int lev = 0; lev < _maxLevel; ++lev) {
		int sCount = 0;
		for (int i = 0; i < _svoLevelsHDD[lev]->size(); ++i) {
			BuildNode bn = _svoLevelsHDD[lev]->item(i);
			if (bn.mirroredX != 0 || bn.mirroredY != 0 || bn.mirroredZ != 0) sCount++;
		}
		printf("LEVEL %d :\t %d nodes with some symm flag\n", lev, sCount);
	}
#endif
	if(!internalCall) printf("OK! [%s]\n", sl::human_readable_duration(_stats.toSDAGTime).c_str());
}


void GeomOctree::toSDAGCanonical() {

	printf("* Transforming DAG -> SDAG (canonical)... "); fflush(stdout);

	struct MirroredNode {
		MirroredNode() : mirrorX(false), mirrorY(false), mirrorZ(false), id(0) {}
		MirroredNode(bool mX, bool mY, bool mZ, id_t id) :
			mirrorX(mX), mirrorY(mY), mirrorZ(mZ), id(id) {}
		bool mirrorX, mirrorY, mirrorZ;
		id_t id;
	};

	std::vector<MirroredNode> correspondences;
	std::map<Node, id_t> uniqueNodesChecker;
	std::vector<Node> uniqueNodes;

	_nNodes = 0;

	_clock.restart();

	for (unsigned int lev = _levels - 1; lev > 0; lev--) {
		unsigned int oldLevSize = (unsigned int)_data[lev].size();
		uniqueNodes.clear();
		uniqueNodesChecker.clear();
		correspondences.resize(oldLevSize);

		for (id_t i = 0; i < _data[lev].size(); i++) {

			Node n = _data[lev][i];
			if (!n.hasChildren()) continue; // skip empty nodes

			// Set invariances
			if (n == n.mirror(true, false, false)) n.setInvariantBit(0);
			if (n == n.mirror(false, true, false)) n.setInvariantBit(1);
			if (n == n.mirror(false, false, true)) n.setInvariantBit(2);

			bool cx, cy, cz;
			Node nC = n.getCanonical(cx, cy, cz);

			// Check children invariances
			if (lev < (_levels - 1)) invertInvs(nC, lev, cx, cy, cz);

			if (uniqueNodesChecker.find(nC) != uniqueNodesChecker.end()) {
				correspondences[i] = MirroredNode(cx, cy, cz, uniqueNodesChecker[nC]);
			} else { // !found
				uniqueNodesChecker[nC] = (id_t)uniqueNodes.size();
				correspondences[i] = MirroredNode(cx, cy, cz, id_t(uniqueNodes.size()));
				uniqueNodes.push_back(nC);
			}
		}
		_data[lev].clear();
		_data[lev] = uniqueNodes;
		_data[lev].shrink_to_fit();
		_nNodes += _data[lev].size();

		for (id_t i = 0; i < _data[lev - 1].size(); i++) {
			Node &n = _data[lev - 1][i];
			for (id_t j = 0; j < 8; j++) {
				if (n.existsChildPointer(j)) {
					MirroredNode mn = correspondences[n.children[j]];
					n.children[j] = mn.id;
					if (mn.mirrorX) n.setChildMirrrorBit(0, j);
					if (mn.mirrorY) n.setChildMirrrorBit(1, j);
					if (mn.mirrorZ) n.setChildMirrrorBit(2, j);
				}
			}
		}
	}

	_stats.toSDAGTime = _clock.elapsed();
	_stats.nNodesSDAG = _nNodes;
	_state = S_SDAG;

	printf("OK! [%s]\n", sl::human_readable_duration(_stats.toSDAGTime).c_str());
}


bool GeomOctree::checkIntegrity() {
	struct QueueItem {
		QueueItem(id_t id, sl::uint8_t l) : nodeID(id), level(l) {}
		id_t nodeID;
		sl::uint8_t level;
	};
	std::queue<QueueItem> queue;

	queue.push(QueueItem(0, 0));

	unsigned int nVisitedNodes = 0;

	while (!queue.empty()) {
		QueueItem qi = queue.front(); queue.pop();
		nVisitedNodes++;
		if (qi.level >= _levels) return false;
		if (qi.nodeID >= _data[qi.level].size()) return false;
		Node &n = _data[qi.level][qi.nodeID];
		if (!n.hasChildren()) return false;
		for (int i = 0; i < 8; ++i) {
			if (qi.level < _levels - 1) {
				if (n.existsChild(i) && !n.existsChildPointer(i)) return false;
				if (!n.existsChild(i) && n.existsChildPointer(i)) return false;
				if (n.existsChild(i)) queue.push(QueueItem(n.children[i], qi.level+1));
			}
		}
	}

	for (unsigned int i = 0; i < _levels; i++) {
		printf("%u\t%zu\n", i, _data[i].size());
	}

	printf("Visited %i nodes\n", nVisitedNodes);

	return true;
}

/**
 * @brief GeomOctree::invertInvs Unsets mirror bits for axis' that are invariant for child nodes
 * @param n
 * @param lev
 * @param inX
 * @param inY
 * @param inZ
 */
void GeomOctree::invertInvs(Node &n, int lev, bool inX, bool inY, bool inZ) {
	for (int i = 0; i < 8; ++i) {
		if (n.existsChildPointer(i)) {
			Node & nc = _data[lev + 1][n.children[i]];
			if (inX && nc.getInvariantBit(0)) n.unsetChildMirrrorBit(0, i);
			if (inY && nc.getInvariantBit(1)) n.unsetChildMirrrorBit(1, i);
			if (inZ && nc.getInvariantBit(2)) n.unsetChildMirrrorBit(2, i);
		}
	}
}

size_t GeomOctree::getMemFootprint() {
	size_t mem = sizeof(Octree);
	for (int i = 0; i < _data.size(); ++i) {
		mem += _data[i].capacity()*sizeof(Node);
		mem += sizeof(std::vector<Node>);
	}
	return mem;
}

int GeomOctree::traverse(sl::point3f p) const
{
	if (_state == S_EMPTY) return -1;

	if (!_bbox.contains(p)) return -1;

	if (_state == S_SVO || _state == S_DAG) {

		const Node *n = &_data[0][0];
		sl::point3f nodeCenter = _bbox.center();
		unsigned int nodeLevel = 0;

		while (nodeLevel < _levels) {
			int selectedChild = 4 * int(p[0] > nodeCenter[0]) + 2 * int(p[1] > nodeCenter[1]) + int(p[2] > nodeCenter[2]);
			if (!n->existsChild(selectedChild)) break;
			nodeLevel++;
			float hs = getHalfSide(nodeLevel);
			nodeCenter[0] += (p[0] > nodeCenter[0]) ? hs : -1.f * hs;
			nodeCenter[1] += (p[1] > nodeCenter[1]) ? hs : -1.f * hs;
			nodeCenter[2] += (p[2] > nodeCenter[2]) ? hs : -1.f * hs;
			if (nodeLevel < _levels) n = &_data[nodeLevel][ n->children[selectedChild] ];
		}
		return nodeLevel;
	}

	if (_state == S_SDAG) {

		const Node *n = &_data[0][0];
		sl::point3f nodeCenter = _bbox.center();
		unsigned int nodeLevel = 0;
		bool mX = false, mY = false, mZ = false;
		
		while (nodeLevel < _levels) {
			if (mX) p[0] = 2.f * nodeCenter[0] - p[0];
			if (mY) p[1] = 2.f * nodeCenter[1] - p[1];
			if (mZ) p[2] = 2.f * nodeCenter[2] - p[2];
			int selectedChild = 4 * int(p[0] > nodeCenter[0]) + 2 * int(p[1] > nodeCenter[1]) + int(p[2] > nodeCenter[2]);
			if (!n->existsChild(selectedChild)) break;
			nodeLevel++;
			float hs = getHalfSide(nodeLevel);
			nodeCenter[0] += (p[0] > nodeCenter[0]) ? hs : -1.f * hs;
			nodeCenter[1] += (p[1] > nodeCenter[1]) ? hs : -1.f * hs;
			nodeCenter[2] += (p[2] > nodeCenter[2]) ? hs : -1.f * hs;
			mX = n->getChildMirrrorBit(0, selectedChild);
			mY = n->getChildMirrrorBit(1, selectedChild);
			mZ = n->getChildMirrrorBit(2, selectedChild);

			if (nodeLevel < _levels) n = &_data[nodeLevel][n->children[selectedChild]];
		}
		return nodeLevel;
	}

	return -1;
}

