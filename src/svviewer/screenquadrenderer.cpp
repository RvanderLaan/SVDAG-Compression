
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


#include "screenquadrenderer.hpp"


const GLfloat ScreenQuadRenderer::_vboData[] = { 
	-1.0f, -1.0f, 0.0f,
	 1.0f, -1.0f, 0.0f,
	-1.0f,  1.0f, 0.0f,
	 1.0f,  1.0f, 0.0f,
};

void ScreenQuadRenderer::initGL() {

	glGenVertexArrays(1, &_glVertexArray);
	glBindVertexArray(_glVertexArray);

	glGenBuffers(1, &_glvbo);
	glBindBuffer(GL_ARRAY_BUFFER, _glvbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(_vboData), _vboData, GL_STATIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glBindVertexArray(0);
}

ScreenQuadRenderer::~ScreenQuadRenderer() {
	glDeleteBuffers(1, &_glvbo);
	glDeleteVertexArrays(1, &_glVertexArray);
}

void ScreenQuadRenderer::draw() {

	glBindVertexArray(_glVertexArray);

	glBindBuffer(GL_ARRAY_BUFFER, _glvbo);
	
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);                 

	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4); 

	glDisableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}