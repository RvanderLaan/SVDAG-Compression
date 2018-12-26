
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
#include <sl/linear_map_factory.hpp>
#include <sl/clock.hpp>
#include <symvox/scene.hpp>
#include <symvox/util.hpp>
#include "renderer_monitor.hpp"
#include "camera.hpp"

class Renderer {
public:

	Renderer() {
	}

	virtual void init() = 0;
	virtual void draw(GLint renderBuffer = 0) = 0;
	virtual void clearState() = 0;
	virtual void resetState() = 0;
	virtual void end() = 0;
	virtual sl::aabox3f getSceneBBox() = 0;


	inline void setCamera(Camera * cam) { _camera = cam; }
	inline void setScreenResolution(int width, int height) { _screenRes[0] = float(width); _screenRes[1] = float(height); }
	inline RendererMonitor::FrameStats getStats() { return _rendererMonitor.getStats(); }
	inline void toggleRenderingStats() { if (_rendererMonitor.isEnabled()) _rendererMonitor.stop(); else _rendererMonitor.start(); }

	inline void saveLastFrameBMP(std::string filename) {
		printf("* Saving last FrameBuffer in '%s'... ", filename.c_str());
		int resX = (int)_screenRes[0], resY = (int)_screenRes[1];
		sl::uint8_t * pixels = new sl::uint8_t[3 * resX * resY];
		glReadPixels(0, 0, resX, resY, GL_BGR, GL_UNSIGNED_BYTE, pixels);
		writeBMP(filename, resX, resY, pixels, true);
		printf("OK!\n");
	}

protected:
	Camera * _camera;
	sl::point2f _screenRes;
	RendererMonitor _rendererMonitor;

#define printfGLError() checkGlError(__FILE__, __LINE__)
public:
	inline bool checkGlError(std::string file, unsigned int line) {
#if 1
		GLenum err = glGetError();
		if (err != GL_NO_ERROR) {
			printf("-> GL_ERROR %i: %s [%s, %u]\n", err, gluErrorString(err), file.c_str(), line);
			return false;
		}
#endif
		return true;
	}
	
	inline bool checkCurrentFrameBuffer() {
		GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		if (status != GL_FRAMEBUFFER_COMPLETE) {
			printf("FRAMEBUFFER ERROR! -> "); printfGLError();
			return false;
		}
		return true;
	}
};
