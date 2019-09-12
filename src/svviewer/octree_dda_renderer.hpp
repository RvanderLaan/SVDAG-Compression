
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

#include <sl/linear_map_factory.hpp>
#include <sl/clock.hpp>
#include <symvox/encoded_octree.hpp>
#include "renderer.hpp"
#include "glslprogram.hpp"
#include "renderer_monitor.hpp"
#include "screenquadrenderer.hpp"
#include <GL/glew.h>

class OctreeDDARenderer : public Renderer {
public:
	typedef enum {
		SVO,
		SVDAG,
		USSVDAG,
		SSVDAG 
	} OctreeFormat;

	typedef enum {
		VIEWER,
		DEPTH,
		SHADOW,
		AO
	} RenderMode;

public:
	OctreeDDARenderer(EncodedOctree * eo);
	virtual void uploadData(bool init = false);
	virtual void init();
	virtual void draw(GLint renderBuffer = 0);
	virtual void clearState();
	virtual void resetState();
	virtual void end();
	void clearVoxel(int position, int level);
	virtual unsigned int getNLevels() { return _encodedOctree->getNLevels(); }
	inline void selectRenderMode(RenderMode mode) { _selectedRenderMode = mode; }
	inline RenderMode getSelectedRenderMode() { return _selectedRenderMode; }

	virtual void clearVoxel(int position);
	inline void setSelectedVoxelIndex(int i) { _selectedVoxelIndex = i; }
	inline int getSelectedVoxelIndex() { return _selectedVoxelIndex; }
	inline bool getRandomColors() { return _randomColors; }
	inline void toggleRandomColors() { _randomColors = !_randomColors; }
	inline bool getShadowsEnabled() { return _enableShadows; }
	inline void toggleShadowsEnabled() { _enableShadows = !_enableShadows; }

	inline sl::aabox3f getSceneBBox() { return _encodedOctree->getSceneBBox(); }
	inline void setEncodedOctree(EncodedOctree * eo) { _encodedOctree = eo; }

	inline void setDrawLevel(unsigned int dl) { _drawLevel = sl::max(1U, sl::min(dl, getNLevels())); }
	inline void incDrawLevel(int k = 1) { setDrawLevel(_drawLevel + k); }
	inline void decDrawLevel(int k = 1) { if (_drawLevel - k >= 1) setDrawLevel(_drawLevel - k); }
	inline unsigned int getDrawLevel() { return _drawLevel; }

	inline void incPixelTolerance() { _pixelTolerance *= 2.0f; _projectionFactor = getProjectionFactor(_pixelTolerance); }
	inline void decPixelTolerance() { _pixelTolerance /= 2.0f; _projectionFactor = getProjectionFactor(_pixelTolerance); }
	inline void setPixelTolerance(float f) { _pixelTolerance = f; _projectionFactor = getProjectionFactor(_pixelTolerance); }
	inline float getPixelTolerance() const { return _pixelTolerance; }
	inline void setFovH(float x) { _fovy = x; _projectionFactor = getProjectionFactor(_pixelTolerance); }
	inline float getFovH() { return _fovy; }
	inline float getProjectionFactor() { return _projectionFactor; }

	inline void setGPUTraversalMaxIters(unsigned int mi) { _gpuTraversalMaxIters = sl::max(1U, mi); }
	inline unsigned int getGPUTraversalMaxIters() { return _gpuTraversalMaxIters; }
	inline void incGPUTraversalMaxIters(unsigned int k = 1) { setGPUTraversalMaxIters(_gpuTraversalMaxIters + k); }
	inline void decGPUTraversalMaxIters(unsigned int k = 1) { setGPUTraversalMaxIters(_gpuTraversalMaxIters - k); }

	inline void nextViewerRenderMode() { _viewerRenderMode = (_viewerRenderMode < 3) ? _viewerRenderMode + 1 : 0; }
	inline void setViewerRenderMode(int mode) { _viewerRenderMode = mode; }
	inline int getViewerRenderMode() { return _viewerRenderMode; }
	inline void setLightPos(sl::point3f p) { _lightPos = p; }
	inline sl::point3f getLightPos() { return _lightPos; }

	inline void toggleUseMinDepthOptimization() { _useMinDepthOptimization = !_useMinDepthOptimization; }
	inline bool getUseMinDepthOptimization() { return _useMinDepthOptimization; }

	void generateHSSamples(int n);
	inline void setShadowsInputTexs(GLuint depthTex, GLuint posTex, GLuint normTex) {
		_shadowsInputDepthTex = depthTex;
		_shadowsInputPosTex = posTex;
		_shadowsInputNormTex = normTex;
	}

	std::string getOctreeFormatName();
	std::string getRenderModeName(RenderMode mode);
	float getProjectionFactor(float pixelTolerance, float screenDivisor = 1.0f);

protected:
	EncodedOctree* _encodedOctree;
	ScreenQuadRenderer _sqr;
	GLSLProgram _program[4];
	RenderMode _selectedRenderMode;
	OctreeFormat _octreeFormat;
	GLuint _glTex[4], _glBuf[2], _glFBO;
	
	GLuint _shadowsInputDepthTex, _shadowsInputPosTex, _shadowsInputNormTex;

	bool _sdag4UseTex3D;

    bool _useMinDepthOptimization{ true };
	
	// params mode VIEWER
	int _viewerRenderMode;
	int _selectedVoxelIndex;
	bool _randomColors;
	bool _enableShadows;

	// params mode SHADOW
	sl::point3f _lightPos;

	// parms mode AO
	std::vector<sl::vector2f> _hsSamples2D;
	std::vector<sl::vector3f> _hsSamples3D;
	int _numAORays;
	float _lengthAORays;


	unsigned int _gpuTraversalMaxIters;
	int _drawLevel;
	float _pixelTolerance;
	float _fovy;
	float _projectionFactor;

	GLuint generateChildrenIndirectionTex();
};
