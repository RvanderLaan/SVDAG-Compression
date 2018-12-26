
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

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <sl/fixed_size_point.hpp>
#include <sl/rigid_body_map.hpp>
#include <sl/projective_map.hpp>
#include <vector>

class Camera {
public:
	typedef enum {
		Z_UP,
		Y_UP
	} TCoordinateSystem;

	typedef enum {
		FREE,
		WT_RECORDING,
		WT_PLAYING
	} TState;
	
	struct Config {
		sl::point3f pos;
		sl::point3f target;
		TCoordinateSystem coordSyst;
		sl::vector3f up;
		float fov;
		float nearClip, farClip;
	};

	 struct WTKey {
		 WTKey(sl::point3f p, sl::point3f t) : pos(p), target(t) {}
		 sl::point3f pos;
		 sl::point3f target;
	 };

	 typedef struct {
		 std::vector<WTKey> keys;
		 size_t currentKey;
		 std::string filename;
		 void reset() { keys.clear(); currentKey = 0; filename.clear(); }
	 } WalkThrough;



public:
	Camera(GLFWwindow * win);
	void update(bool sticky = true);
	
	void setInitCamera(sl::point3f pos, sl::point3f target, TCoordinateSystem coordSyst, float fov, float nearClip = 0.1f, float farClip = 1000.0f);
	inline void setZoomFactor(float zf) { _zoomFactor = zf; }
	inline void setRotFactor(float rf) { _rotFactor = rf; }
	inline void setPanFactor(float pf) { _panFactor = pf; }
	inline void setWalkFactor(float wf) { _walkFactor = wf; }

	inline sl::projective_map3f getProjMatrix() { return _projMatrix; }
	inline sl::rigid_body_map3f getViewMatrix() { return _viewMatrix; }
	inline sl::projective_map3f getProjMatrixInv() { return _projMatrixInv; }
	inline sl::rigid_body_map3f getViewMatrixInv() { return _viewMatrixInv; }
	inline Config getCurrentConfig() { return _cam; }

	void recordWalkthrough(std::string filename = "");
	void loadWalkthrough(std::string filename);
	bool playWalkthrough(int interpol = 1, bool rewind = false);
	void stopWalkthrough();
	inline bool isPlayingWalkthrough() { return (_state == WT_PLAYING); }
	
	void printControls();

private:
	typedef enum {
		NONE,
		ZOOMING,
		ROTATING,
		PANNING
	} TGesture;

	void updateMatrices();
	void reset();

	sl::rigid_body_map3f _viewMatrix, _viewMatrixInv;
	sl::projective_map3f _projMatrix, _projMatrixInv;
	GLFWwindow * _glfwWindow;
	float _winWidth, _winHeight;
	Config _cam, _initCam;
	float _rotRadius;
	sl::point2f _coords, _clickedCoords, _lastCoords;
	float _zoomFactor, _rotFactor, _panFactor, _walkFactor;
	TState _state{ FREE };
	TGesture _gestureState{ NONE };
	WalkThrough _walkthrough;
	int _wtInterpol { 1 };
};
