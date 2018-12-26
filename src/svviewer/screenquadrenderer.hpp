
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
#include <string>

class ScreenQuadRenderer {

public:
	ScreenQuadRenderer() {}
	void initGL();
	~ScreenQuadRenderer();
	void draw();
	inline std::string getDefaultVertexProgram() const {
		return
			"#version 420 core											\n"
			"layout(location = 0) in vec3 vertexPosition_modelspace;	\n"
			"void main() {												\n"
			"	gl_Position.xyz = vertexPosition_modelspace;			\n"
			"	gl_Position.w = 1.0;									\n"
			"}															\n";
	}

protected:
	static const GLfloat _vboData[];
	GLuint _glVertexArray;
	GLuint _glvbo;
};
