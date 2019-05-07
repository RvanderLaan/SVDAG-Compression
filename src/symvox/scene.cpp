
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

#include "scene.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ctime>
#include <stdio.h>
#include <sys/stat.h>

#if BUILD_LASPARSER
#include <liblas/liblas.hpp>
#endif

#ifndef _MSC_VER
#	include <string.h>
#	define strtok_s strtok_r
#endif

#define ETA(x) ((clock() - x) / float(CLOCKS_PER_SEC))

Scene::Scene() {
	_nTris = 0;
	_bbox.to_empty();

	TMaterial defaultMat;
    sprintf(defaultMat.name, "Voxelator Default Mat");
	defaultMat.ambientColor = sl::color3f(0.1, 0.1, 0.1);
	defaultMat.diffuseColor = sl::color3f(0.9, 0.9, 0.9);
	defaultMat.specColor = sl::color3f(0.2, 0.2, 0.2);
	_materials.push_back(defaultMat);
}
void Scene::loadObj(std::string fileName, bool tryLoadBinCache, bool loadMaterials, bool loadNormals, bool recomputeNormals, bool doBuildTriVector) {
	
	printf("* Loading '%s'...\n", fileName.c_str());

	if (tryLoadBinCache) {
		std::string binCacheFilename = fileName + ".bincache";
		std::ifstream ifsBinCache(binCacheFilename);
		if (ifsBinCache.good()) {
			printf("\t- Binary cache found! Loading '%s'... ", binCacheFilename.c_str());
			loadBinObj(binCacheFilename);
			printf("OK!\n");
			if (doBuildTriVector) {
				printf("\t- Building tri vector... ");
                buildTriVector(false);
				printf("OK!\n");
			}
			printf("\t- Loaded %s triangles\n", sl::human_readable_quantity(_triangles.size()).c_str());
			printf("\t- Bbox: [%.3f %.3f %.3f]  [%.3f %.3f %.3f]\n", _bbox[0][0], _bbox[0][1], _bbox[0][2], _bbox[1][0], _bbox[1][1], _bbox[1][2]);
			return;
		}
	}
	printf("\t- Reading ASCII OBJ (reading lines)...");
	FILE * file = fopen(fileName.c_str(), "r");
	if (file == NULL) {
		std::cout << "ERROR: Scene:loadOBJ: Can't open file " << fileName << std::endl;
		exit(1);
	}

	// Size estimation
#if 0
	fseek(file, 0L, SEEK_END);
	size_t totalSize = ftell(file);
	fseek(file, 0L, SEEK_SET);
	_indexedTris.reserve(totalSize / 10000);
	_vertices.reserve(totalSize / 10000);
	_normals.reserve(totalSize / 10000);
#endif

	clock_t tini = clock();
	_bbox.to_empty();

	_nTris = 0;
	std::size_t currentMat = 0; // default material
	size_t nLine = 0, recomputedNormals = 0;
	const int len = 10000;
	char str[len];
	while (fgets(str, len, file) != NULL) {
		nLine++;
        if (nLine % 1000000L == 0) {
			printf(".%zuM.", nLine / 1000000L); fflush(stdout);
			//printf("\b\b\b\b\b%3d %%", int(100 * (float(ftell(file)) / float(totalSize)))); fflush(stdout);
			//printf("\t[%.2fs] Loaded %iM lines (%.2fM tris, %.2fM verts)\n", ETA(tini), nRead / 1000000, _nTris/1000000.0f, _vertices.size()/1000000.0f);
			//printf("\tBbox: [%.4f, %.4f, %.4f]  [%.4f, %.4f, %.4f]\n\n", _bbMin.x, _bbMin.y, _bbMin.z, _bbMax.x, _bbMax.y, _bbMax.z);
        }

		switch (str[0]) {
		case 'v': { //vertex data
			switch (str[1]) {
			case 'n': // normal
			{
				if (loadNormals) {
					sl::vector3f v;
					sscanf(str + 2, "%f %f %f", &v[0], &v[1], &v[2]);
					_normals.push_back(v);
				}
				break;
			}
            case 't': // tex coord
			{
                // Yes dealing with this
                if (loadMaterials) {
                    sl::vector2f v;
                    sscanf(str + 2, "%f %f", &v[0], &v[1]);
                    _texCoords.push_back(v);
                }
				break;
			}
			default: //vertex
			{
				sl::point3f p;
				int k = sscanf(str + 1, "%f %f %f", &p[0], &p[1], &p[2]);
				if (k != 3) {
					printf("Read in line %zu a bad vertex format\n", nLine);
					break;
				}
				_bbox.merge(p);
				_vertices.push_back(p);
				break;
			}
			}
			break;
		}
		case 'f': { // face
		    // Format doc: https://en.wikipedia.org/wiki/Wavefront_.obj_file
			std::vector<std::string> triples;
			char *next_token;
			char *pch = strtok_s(str + 1, " ", &next_token);
			while (pch) {
				std::string s(pch);
				triples.push_back(s);
				pch = strtok_s(NULL, " ", &next_token);
			}
			// throw away last symbol (\n)
			triples.back().resize(triples.back().size() - 1);
			
			TIndexedTri itri;
			itri.material = currentMat;
			bool noNormal = false;
			for (int i = 0; i < 3; ++i) {
				int idxV=0, idxT=0, idxN=0;
                bool hasNormWithoutTex = triples[i].find("//") != std::string::npos;
				int n = hasNormWithoutTex
				        ? sscanf(triples[i].c_str(), "%d//%d", &idxV, &idxN)
				        : sscanf(triples[i].c_str(), "%d/%d/%d", &idxV, &idxT, &idxN);
				// negative indices
				if (idxV < 0) idxV = int(_vertices.size()) + idxV + 1;
				if (idxT < 0) idxT = int(_texCoords.size()) + idxT + 1;
				if (idxN < 0) idxN = int(_normals.size()) + idxN + 1;
				if (n == 0) {
					std::cout << "ERROR:Scene:loadObj: problem parsing face instruction" << std::endl;
				} else if (n == 1) { // just vertex info
					itri.verticesIdx[i] = idxV - 1;
					noNormal = true;
                    itri.texCoordIdx[i] = 0; // all texCoordIdx of 0 indicates that this triangle does not have a tex
				} else if (n == 2 && hasNormWithoutTex) { // vertex and normal info
					itri.verticesIdx[i] = idxV - 1;
					if (!recomputeNormals) itri.normalsIdx[i] = idxN - 1;
                    itri.texCoordIdx[i] = 0;
                } else if (n == 2) { // vertex and tex coord info
                    itri.verticesIdx[i] = idxV - 1;
                    if (loadMaterials) itri.texCoordIdx[i] = idxT - 1;
                    noNormal = true;
				} else { // vertex, normal and tex info
                    itri.verticesIdx[i] = idxV - 1;
                    if (!recomputeNormals) itri.normalsIdx[i] = idxN - 1;
                    if (loadMaterials) itri.texCoordIdx[i] = idxT - 1;
				}
			}

			if (loadNormals && (noNormal || recomputeNormals)) {
				sl::vector3f e1 = _vertices[itri.verticesIdx[1]] - _vertices[itri.verticesIdx[0]];
				sl::vector3f e2 = _vertices[itri.verticesIdx[2]] - _vertices[itri.verticesIdx[0]];
				sl::vector3f n = e1.cross(e2);
				itri.normalsIdx[0] = itri.normalsIdx[1] = itri.normalsIdx[2] = _normals.size();
				_normals.push_back(n);
				recomputedNormals++;
			}
			_indexedTris.push_back(itri);
			_nTris++;

			if (triples.size() > 4) { // if it has 4 indices (a quad)
				for (int i = 0; i < 3; ++i) {
					int idxV = 0, idxN = 0;
					int n = sscanf(triples[(i + 2) % 4].c_str(), "%d//%d", &idxV, &idxN);
					// negative indices
					if (idxV < 0) idxV = int(_vertices.size()) + idxV + 1;
					if (idxN < 0) idxN = int(_normals.size()) + idxN + 1;
					if (n == 0) {
						std::cout << "ERROR:Scene:loadObj: problem parsing face instruction" << std::endl;
					}
					else if (n == 1) { // just vertex info
						itri.verticesIdx[i] = idxV - 1;
						noNormal = true;
					}
					else if (n == 2) { // vertex and normal info
						itri.verticesIdx[i] = idxV - 1;
						if (!recomputeNormals) itri.normalsIdx[i] = idxN - 1;
					}
				}

				if (loadNormals && (noNormal || recomputeNormals)) {
					sl::vector3f e1 = _vertices[itri.verticesIdx[1]] - _vertices[itri.verticesIdx[0]];
					sl::vector3f e2 = _vertices[itri.verticesIdx[2]] - _vertices[itri.verticesIdx[0]];
					sl::vector3f n = e1.cross(e2);
					itri.normalsIdx[0] = itri.normalsIdx[1] = itri.normalsIdx[2] = _normals.size();
					_normals.push_back(n);
					recomputedNormals++;
				}
				_indexedTris.push_back(itri);
				_nTris++;

                // TODO: Fix texture coordinates for quads. Not dealing with this now, just assume the model is triangulated
            }


			break;
		}
		case 'm': { // mtllib
			if (loadMaterials) {
				std::string path = sl::pathname_directory(fileName);
				std::string mtlFileName(str + 7);
				// throw away linebreak character
				mtlFileName.resize(mtlFileName.size() - 1);
				mtlFileName.erase(std::remove(mtlFileName.begin(), mtlFileName.end(), '\r'), mtlFileName.end());
				loadObjMtlLib(path + "/" + mtlFileName);
			}
			break;
		}
		case 'u': { // usemtl
			if (loadMaterials) {
				std::string mtlName(str + 7);
				mtlName.resize(mtlName.size() - 1);

				bool found = false;
				for (std::size_t i = 0; i < _materials.size(); i++) {
                    if (strncmp(_materials[i].name, mtlName.c_str(), mtlName.size()) == 0) {
						currentMat = i;
						found = true;
						break;
					}
				}
				if (!found) std::cout << "WARNING:Scene:loadObj(): Material " << mtlName << " not found in library" << std::endl;
			}
		}
		}
	}

	printf(" OK! [%.2f s]    \n", ETA(tini));
	printf("\t- Loaded %zu triangles (%zu recomputed normals, %zu materials) in %.2fM obj lines\n", _nTris, recomputedNormals, _materials.size(), nLine / 1000000.0f);
	printf("\t- Bbox: [%.3f %.3f %.3f]  [%.3f %.3f %.3f]\n", _bbox[0][0], _bbox[0][1], _bbox[0][2], _bbox[1][0], _bbox[1][1], _bbox[1][2]);
	fflush(stdout);

	_vertices.shrink_to_fit();
	_normals.shrink_to_fit();
	_materials.shrink_to_fit();
	_indexedTris.shrink_to_fit();

	printf("\t- Saving '%s' binary cache... ", (fileName + ".bincache").c_str());
	saveBinObj(fileName + ".bincache");
	printf("OK!\n");

	fclose(file);

	if (doBuildTriVector) {
		printf("\t- Building tri vector... ");
		buildTriVector(true);
		printf("OK!\n");
	}
}


bool Scene::loadObjMtlLib(std::string fileName) {

	std::ifstream ifs;
	ifs.open(fileName);
	printf(" [Loading %s... ", fileName.c_str());
	if (!ifs.is_open()) {
		printf("WARNING: Couldn't open file.\n");
		return false;
	}

    std::string line, cmd, tmp;
    TMaterial * mat = NULL;
	std::size_t n = 0;
	while (!ifs.eof()) {
		std::getline(ifs, line);
		if (line.empty() || line.find_first_of("#") != std::string::npos)
			continue;
		std::stringstream ss(line);
		ss >> cmd;
        if (cmd == "newmtl") {
            if (mat) _materials.push_back(*mat);
            mat = new TMaterial();
            ss >> tmp;
            strcpy(mat->name, tmp.c_str());
		}
		else if (cmd == "Ka") {
			sl::color3f c;
			ss >> c[0] >> c[1] >> c[2];
            mat->ambientColor = c;
		}
		else if (cmd == "Kd") {
			sl::color3f c;
			ss >> c[0] >> c[1] >> c[2];
            mat->diffuseColor = c;
		}
		else if (cmd == "Ks") {
			sl::color3f c;
			ss >> c[0] >> c[1] >> c[2];
            mat->specColor = c;
		}
		else if (cmd == "Ns") {
            ss >> mat->specCoef;
		}
		else if (cmd == "d" || cmd == "Ts") {
            ss >> mat->transparency;
		}
        else if (cmd == "map_Kd") { // diffuse texture
            ss >> tmp;
            strcpy(mat->texture, tmp.c_str());
        }
	}
    if (mat)  _materials.push_back(*mat);
	printf("OK! (%zu materials loaded)] ", _materials.size());
	ifs.close();
	return true;
}

void Scene::saveBinObj(std::string filename) {

	BinObjHeader header;
	header.nVertices = _vertices.size();
	header.nNormals = _normals.size();
	header.nMaterials = _materials.size();
    header.nTexCoords = _texCoords.size();
	header.nIndexedTris = _indexedTris.size();
	header.bboxMin[0] = _bbox[0][0];
	header.bboxMin[1] = _bbox[0][1];
	header.bboxMin[2] = _bbox[0][2];
	header.bboxMax[0] = _bbox[1][0];
	header.bboxMax[1] = _bbox[1][1];
	header.bboxMax[2] = _bbox[1][2];
	
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	file.write((char *)(&header), sizeof(BinObjHeader));
	
	file.write((char*)&_vertices[0], _vertices.size() * sizeof(sl::point3f));
	if (_normals.size() > 0) file.write((char*)&_normals[0], _normals.size() * sizeof(sl::vector3f));

	file.write((char*)&_materials[0], _materials.size() * sizeof(TMaterial));

    if (_texCoords.size() > 0) file.write((char*)&_texCoords[0], _texCoords.size() * sizeof(sl::vector2f));
	file.write((char*)&_indexedTris[0], _indexedTris.size() * sizeof(TIndexedTri));

	file.close();	
}

void Scene::loadBinObj(std::string filename) {

	std::ifstream file(filename, std::ios::in | std::ios::binary);
	BinObjHeader header;
	file.read((char *)(&header), sizeof(BinObjHeader));
	
	memcpy(&_bbox[0], &header.bboxMin, 3 * sizeof(float));
	memcpy(&_bbox[1], &header.bboxMax, 3 * sizeof(float));
	_nTris = header.nIndexedTris;

	_vertices.resize(header.nVertices);
	file.read((char *)(&_vertices[0]), header.nVertices*sizeof(sl::point3f));

	_normals.resize(header.nNormals);
	if (header.nNormals > 0) file.read((char *)(&_normals[0]), header.nNormals*sizeof(sl::vector3f));

	_materials.resize(header.nMaterials);
	file.read((char *)(&_materials[0]), header.nMaterials*sizeof(TMaterial));

    _texCoords.resize(header.nTexCoords);
    if (header.nTexCoords > 0) file.read((char *)(&_texCoords[0]), header.nTexCoords*sizeof(sl::vector2f));

	_indexedTris.resize(header.nIndexedTris);
	file.read((char *)(&_indexedTris[0]), header.nIndexedTris*sizeof(TIndexedTri));

	file.close();
}

void Scene::testOutput(std::string filename)
{
	FILE *f;
	f = fopen(filename.c_str(), "w");

	for (unsigned int i = 0; i < _triangles.size()/3; i++) {
		sl::color3f mat = _materials[_indexedTris[i].material].diffuseColor;

		fprintf(f, "v %f %f %f %f %f %f\n", _triangles[i * 3 + 0][0], _triangles[i * 3 + 0][1], _triangles[i * 3 + 0][2], mat[0], mat[1], mat[2]);
		fprintf(f, "v %f %f %f %f %f %f\n", _triangles[i * 3 + 1][0], _triangles[i * 3 + 1][1], _triangles[i * 3 + 1][2], mat[0], mat[1], mat[2]);
		fprintf(f, "v %f %f %f %f %f %f\n", _triangles[i * 3 + 2][0], _triangles[i * 3 + 2][1], _triangles[i * 3 + 2][2], mat[0], mat[1], mat[2]);

		fprintf(f, "f %i %i %i\n", i * 3 + 1, i * 3 + 2, i * 3 + 3);
	}
	
	fclose(f);
}

void Scene::buildTriVector(bool clearOtherData) {
	_triangles.clear();
	for (TIndexedTri& itri : _indexedTris) {
		if (itri.verticesIdx[0] < _vertices.size() &&
			itri.verticesIdx[1] < _vertices.size() &&
			itri.verticesIdx[2] < _vertices.size()) {
			_triangles.push_back(_vertices[itri.verticesIdx[0]]);
			_triangles.push_back(_vertices[itri.verticesIdx[1]]);
			_triangles.push_back(_vertices[itri.verticesIdx[2]]);
		}
	}

	if (clearOtherData) {
		_vertices.clear(); _vertices.shrink_to_fit();
		_normals.clear(); _normals.shrink_to_fit();
        //_materials.clear(); _materials.shrink_to_fit(); // causes some problems with building the attr dag
        //_texCoords.clear(); _texCoords.shrink_to_fit();
		//_indexedTris.clear(); _indexedTris.shrink_to_fit();
	}
}

void Scene::loadLas(std::string fileName) {
	#if BUILD_LASPARSER

    printf("* Loading '%s'...\n", fileName.c_str());

    // Parsing LAS based on https://liblas.org/tutorial/cpp.html
    std::ifstream ifs;
    ifs.open(fileName, std::ios::in | std::ios::binary);

    liblas::ReaderFactory f;
    liblas::Reader reader = f.CreateWithStream(ifs);

    // After the reader has been created, you can access members of the Public Header Block
    liblas::Header const& header = reader.GetHeader();

    std::cout << "Compressed: " << ((header.Compressed() == true) ? "true" : "false");
    std::cout << "Signature: " << header.GetFileSignature() << '\n';
    std::cout << "Points count: " << header.GetPointRecordsCount() << '\n';

    _bbox[0] = sl::point3f(header.GetMinX(), header.GetMinY(), header.GetMinZ());
    _bbox[1] = sl::point3f(header.GetMaxX(), header.GetMaxY(), header.GetMaxZ());
    printf("\t- Bbox: [%.3f %.3f %.3f]  [%.3f %.3f %.3f]\n", _bbox[0][0], _bbox[0][1], _bbox[0][2], _bbox[1][0], _bbox[1][1], _bbox[1][2]);
	#endif
}
