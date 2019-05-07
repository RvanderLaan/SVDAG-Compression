
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

#include <vector>
#include <string>

#include <sl/fixed_size_point.hpp>
#include <sl/fixed_size_vector.hpp>
#include <sl/axis_aligned_box.hpp>

#include <CImg.h>

using namespace cimg_library;

class Scene {
public:
	typedef struct _Material {
        char name[256];
		sl::color3f diffuseColor;
		sl::color3f specColor;
		sl::color3f ambientColor;
		float specCoef;
		float transparency;
        char texture[256];
	} TMaterial;

	typedef struct _IndexedTri {
		std::size_t verticesIdx[3];
		std::size_t normalsIdx[3];
        std::size_t texCoordIdx[3];
		std::size_t material;
	} TIndexedTri;

	typedef struct _Triangle {
		sl::point3f v0, v1, v2;
		sl::vector3f n0, n1, n2;
		//sl::color3f c0, c1, c2;
	} TTriangle;

public:
	Scene();
	
	void loadObj(std::string fileName, bool tryLoadBinCache = true, bool loadMaterials = true, bool loadNormals = false, bool recomputeNormals = false, bool doBuildTriVector = false);
    void loadLas(std::string fileName);

	inline std::vector<sl::point3f> * getVertices() { return &_vertices; }
	inline std::vector<sl::vector3f> * getNormals() { return &_normals; }
	inline std::vector<TMaterial> * getMaterials() { return &_materials; }
    inline std::vector<sl::vector2f> * getTexCoords() { return &_texCoords; }
	inline std::vector<TIndexedTri> * getIndexedGeom() { return &_indexedTris; }

	inline std::size_t getNTriangles() { return _indexedTris.size()>0 ? _indexedTris.size() : _triangles.size(); }
	inline std::size_t getNVertices() { return _vertices.size(); }

	inline void getBounds(sl::point3f & min, sl::point3f & max) const { min = _bbox[0]; max = _bbox[1]; }
	inline sl::aabox3f getAABB() const { return _bbox; }
	
	inline void getTriangleVertices(std::size_t idTri, sl::point3f & v0, sl::point3f & v1, sl::point3f & v2) {
		v0 = _vertices[_indexedTris[idTri].verticesIdx[0]];
		v1 = _vertices[_indexedTris[idTri].verticesIdx[1]];
		v2 = _vertices[_indexedTris[idTri].verticesIdx[2]];
	}
	inline void getTriangleNormals(std::size_t idTri, sl::vector3f & n0, sl::vector3f & n1, sl::vector3f & n2) {
		n0 = _normals[_indexedTris[idTri].normalsIdx[0]];
		n1 = _normals[_indexedTris[idTri].normalsIdx[1]];
		n2 = _normals[_indexedTris[idTri].normalsIdx[2]];
	}
	inline void getTriangleColor(std::size_t idTri, sl::color3f & c) {
		c = _materials[getTriangleMaterialId(idTri)].diffuseColor;
	}
	inline size_t getTriangleMaterialId(std::size_t idTri) {
		return (idTri < _indexedTris.size()) ? _indexedTris[idTri].material : 0;
    }
    inline void getTriangleTexCoords(std::size_t idTri, sl::vector2f & t0, sl::vector2f & t1, sl::vector2f & t2) {
        t0 = _texCoords[_indexedTris[idTri].texCoordIdx[0]];
        t1 = _texCoords[_indexedTris[idTri].texCoordIdx[1]];
        t2 = _texCoords[_indexedTris[idTri].texCoordIdx[2]];
    }
    inline bool isTriangleTextured(std::size_t idTri) {
	    return _indexedTris[idTri].texCoordIdx[0] == -1
	        ||  _indexedTris[idTri].texCoordIdx[1] == -1
	        ||  _indexedTris[idTri].texCoordIdx[2] == -1;
	}
    inline void getTexColor(const std::string texName, const sl::vector2f & uv, sl::color3f & c) {
	    if (_textures.count(texName) == 0) {
	        printf("Texture not found: %s", texName.c_str());
	        return;
	    }
        auto texture = _textures[texName];
	    int w = texture.width();
	    int h = texture.height();
	    int x = uv[0] * w;
	    int y = uv[1] * h;
	    // Todo: this assumes the texture is RGB
	    unsigned char r = texture(x, y, 0);
	    unsigned char g = texture(x, y, 1);
	    unsigned char b = texture(x, y, 2);

	    c[0] = (float) r;
	    c[1] = (float) g;
	    c[2] = (float) b; // todo: / 255?
    }

	void saveBinObj(std::string filename);
	void loadBinObj(std::string filename);

	inline const sl::point3f * getTrianglePtr(const std::size_t id) const { return (id*3 < _triangles.size()) ? &_triangles[id*3] : NULL; }
	inline std::size_t getNRawTriangles() { return _triangles.size()/3; }

	void testOutput(std::string filename);

	inline void loadTextures(std::string path) {
	    _materials.clear();
	    for (auto material : _materials) {
	        std::string texName(material.texture);
	        if (texName.size() > 3 && this->_textures.count(texName) == 0) {
	            // TODO: Check if absolute or relative path - now assumes relative
	            printf("Loading texture '%s'...\n", texName.c_str());
                std::string texPath = path + "/" + texName;
                std::replace(texPath.begin(), texPath.end(), '\\', '/');
                CImg<unsigned char> texture(texPath.c_str());
                this->_textures[texName] = texture;
	        }
	    }
	}

private:
	std::vector <sl::point3f> _vertices;
	std::vector <sl::vector3f> _normals;
	std::vector < TMaterial > _materials;
    std::vector <sl::vector2f> _texCoords;
	std::vector < TIndexedTri > _indexedTris;

	// optional
	void buildTriVector(bool clearOtherData = false);
	std::vector <sl::point3f> _triangles;

	bool loadObjMtlLib(std::string fileName);

	std::size_t _nTris;
	sl::aabox3f _bbox;

    std::map<std::string, CImg<unsigned char>> _textures;

private:
	typedef struct {
		size_t nVertices;
		size_t nNormals;
		size_t nMaterials;
        size_t nTexCoords;
		size_t nIndexedTris;
		float bboxMin[3];
		float bboxMax[3];
	} BinObjHeader;
};
