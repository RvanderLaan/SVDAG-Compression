
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

class RendererMonitor {

public:
	struct FrameStats {
		FrameStats() : time(0), samples(0) {}
		FrameStats(float time, unsigned int samples) : time(time), samples(samples) {}
		float time;
		unsigned int samples;
	};

public:
	RendererMonitor();
	~RendererMonitor();
	void start(unsigned int bufferSize = 5);
	void stop();
	void beginFrame();
	void endFrame();
	inline FrameStats getStats() { return _lastStats; }
	inline bool isEnabled() { return _enabled; }

private:
	unsigned int _bufferSize;
	FrameStats _lastStats;
	bool _enabled;
	unsigned int _frameCount;
	unsigned int _currentQuery;
	GLuint * _timeQueries;
	GLuint * _samplesQueries;

};
