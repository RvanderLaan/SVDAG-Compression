
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

#include "renderer_monitor.hpp"
#include <iostream>

RendererMonitor::RendererMonitor() {
	_bufferSize = 0;
	_enabled = false;
	_frameCount = 0;
}

RendererMonitor::~RendererMonitor() {
	if(_enabled) stop();
}

void RendererMonitor::start(unsigned int bufferSize) {
	_bufferSize = bufferSize;
	_enabled = true;
	_frameCount = 0;
	_lastStats = FrameStats();
	_timeQueries = new GLuint[_bufferSize];
	_samplesQueries = new GLuint[_bufferSize];
	glGenQueries(_bufferSize, _timeQueries);
	glGenQueries(_bufferSize, _samplesQueries);
}

void RendererMonitor::stop() {
	glDeleteQueries(_bufferSize, _timeQueries);
	glDeleteQueries(_bufferSize, _samplesQueries);
	delete[] _timeQueries;
	delete[] _samplesQueries;
	_frameCount = 0;
	_enabled = false;
}

void RendererMonitor::beginFrame() {
	if (!_enabled) return;
	const int i = _frameCount % _bufferSize;
	glBeginQuery(GL_TIME_ELAPSED, _timeQueries[i]);
	glBeginQuery(GL_SAMPLES_PASSED, _samplesQueries[i]);
}

void RendererMonitor::endFrame() {
	if (!_enabled) return;
	glEndQuery(GL_TIME_ELAPSED);
	glEndQuery(GL_SAMPLES_PASSED);
	_frameCount++;
	if (_frameCount < _bufferSize) return;

	unsigned int count = 0;
	bool found = false;
	do {
		const unsigned int idx = (_frameCount + count) % _bufferSize;
		GLint foundTQ = 0, foundSQ = 0;
		glGetQueryObjectiv(_timeQueries[idx], GL_QUERY_RESULT_AVAILABLE, &foundTQ);
		glGetQueryObjectiv(_samplesQueries[idx], GL_QUERY_RESULT_AVAILABLE, &foundSQ);
		if (!foundTQ || !foundSQ) break;
		count++;
	} while (count < _bufferSize);

	if (count >= _bufferSize) {
		printf("WARNING RendererMonitor:endFrame: Stat GL Queries sooo slow\n");
		return;
	}

	unsigned int idx = (_frameCount + count - 1) % _bufferSize;

	//printf("%i\n", _bufferSize - count);
	GLuint64 tmp;
	glGetQueryObjectui64v(_timeQueries[idx], GL_QUERY_RESULT, &tmp);
	_lastStats.time = tmp / 1000000.0f;
	glGetQueryObjectui64v(_samplesQueries[idx], GL_QUERY_RESULT, &tmp);
	_lastStats.samples = (unsigned int)tmp;
}

