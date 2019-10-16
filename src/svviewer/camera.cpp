
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


#include "camera.hpp"
#define _USE_MATH_DEFINES
#include <sl/linear_map_factory.hpp>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>

#include "imgui.h"


Camera::Camera(GLFWwindow * win) {
	_glfwWindow = win;
	int tmpWidth, tmpHeight;
	glfwGetWindowSize(_glfwWindow, &tmpWidth, &tmpHeight);
	_winWidth = float(tmpWidth);
	_winHeight = float(tmpHeight);
	_zoomFactor = 1.0f;
	_rotFactor = 1.0f;
	_panFactor = 0.6f;
	_walkFactor = 0.1f;
	setInitCamera(sl::point3f(0, 0, 10.0f), sl::point3f(0, 0, 0), Y_UP, 45.0f);
	updateMatrices();
#if 0
	_currentWTKey = 0;
#endif
}

void Camera::scroll_callback(double xoffset, double yoffset) {
	_walkFactor *= yoffset > 0 ? 1.1 : 0.9;
	printf("Walkfactor: %f\n", _walkFactor);
}

void Camera::update(bool sticky) {
    ImGuiIO& io = ImGui::GetIO();

    if (_state == WT_PLAYING) {
		_cam.pos = sl::point3f::zero();
		_cam.target = sl::point3f::zero();
		for (int i = 0; i < _wtInterpol; ++i) {
			_cam.pos += _walkthrough.keys[_walkthrough.currentKey+i].pos.as_vector();
			_cam.target += _walkthrough.keys[_walkthrough.currentKey+i].target.as_vector();
		}
		_cam.pos.as_vector() /= (float)_wtInterpol;
		_cam.target.as_vector() /= (float)_wtInterpol;

		updateMatrices();
		_walkthrough.currentKey++;
		if (_walkthrough.currentKey > _walkthrough.keys.size()-_wtInterpol) {
			_walkthrough.currentKey = 0;
			_state = FREE;
		}
		return;
	}

	double tmpX, tmpY;
	glfwGetCursorPos(_glfwWindow, &tmpX, &tmpY);
	_coords[0] = (float)tmpX; _coords[1] = (float)tmpY;

	if (!io.WantCaptureKeyboard) {
        if(glfwGetKey(_glfwWindow, GLFW_KEY_SPACE) == GLFW_PRESS) {
            reset();
            return;
        }

        if (glfwGetKey(_glfwWindow, GLFW_KEY_U) == GLFW_PRESS) {
            if (_cam.coordSyst == Z_UP) {
                _initCam.coordSyst = Y_UP;
                _initCam.up = sl::vector3f(0, 1, 0);
            }
            else {
                _initCam.coordSyst = Z_UP;
                _initCam.up = sl::vector3f(0, 0, 1);
            }
			if (_initCam.pos.distance_to(_cam.pos) > _walkFactor) {
            	reset();
			}
            return;
        }

		float speedMod = 1;
		if (glfwGetKey(_glfwWindow, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
			speedMod = 5;
		}

        if (glfwGetKey(_glfwWindow, GLFW_KEY_C) == GLFW_PRESS) {
            printf("Camera POS & TARGET: %.3f %.3f %.3f %.3f %.3f %.3f\n",
                _cam.pos[0], _cam.pos[1], _cam.pos[2],
                _cam.target[0], _cam.target[1], _cam.target[2]);
        }

        if (glfwGetKey(_glfwWindow, GLFW_KEY_W) == GLFW_PRESS) {
            sl::vector3f dir = (_cam.target - _cam.pos).ok_normalized();
            _cam.pos += dir * _walkFactor * speedMod;
            _cam.target += dir * _walkFactor * speedMod;
            updateMatrices();
        }
        if (glfwGetKey(_glfwWindow, GLFW_KEY_S) == GLFW_PRESS) {
            sl::vector3f dir = (_cam.target - _cam.pos).ok_normalized();
            _cam.pos -= dir * _walkFactor * speedMod;
            _cam.target -= dir * _walkFactor * speedMod;
            updateMatrices();
        }
        if (glfwGetKey(_glfwWindow, GLFW_KEY_A) == GLFW_PRESS) {
            sl::vector3f dir = (_cam.target - _cam.pos).ok_normalized();
            dir = dir.cross(_cam.up);
            _cam.pos -= dir * _walkFactor * speedMod;
            _cam.target -= dir * _walkFactor * speedMod;
            updateMatrices();
        }
        if (glfwGetKey(_glfwWindow, GLFW_KEY_D) == GLFW_PRESS) {
            sl::vector3f dir = (_cam.target - _cam.pos).ok_normalized();
            dir = -dir.cross(_cam.up);
            _cam.pos -= dir * _walkFactor * speedMod;
            _cam.target -= dir * _walkFactor * speedMod;
            updateMatrices();
        }
        if (glfwGetKey(_glfwWindow, GLFW_KEY_Q) == GLFW_PRESS) {
          sl::vector3f dir = (_cam.target - _cam.pos);
          float delta = 2.0f * 3.14f / 180.0f * 0.5;
          sl::rigid_body_map3f r = sl::linear_map_factory3f::rotation(_cam.up, delta);
          _cam.target = _cam.pos + r * dir;
          updateMatrices();
        }
        if (glfwGetKey(_glfwWindow, GLFW_KEY_E) == GLFW_PRESS) {
          sl::vector3f dir = (_cam.target - _cam.pos);
          float delta = -2.0f * 3.14f / 180.0f * 0.5;
          sl::rigid_body_map3f r = sl::linear_map_factory3f::rotation(_cam.up, delta);
          _cam.target = _cam.pos + r * dir;
          updateMatrices();
        }
    }

    if (io.WantCaptureMouse) {
        return;
    }
	if (!ImGui::IsAnyWindowHovered() && !ImGui::IsAnyItemHovered()) {
		if(glfwGetMouseButton(_glfwWindow, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS) {
			if(_gestureState != ROTATING) {
				_clickedCoords = _coords;
				_gestureState = ROTATING;
			}
		} else if(glfwGetMouseButton(_glfwWindow, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS) {
			if(_gestureState != ZOOMING) {
				_clickedCoords = _coords;
				_gestureState = ZOOMING;
			}
		} else if(glfwGetMouseButton(_glfwWindow, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS) {
			if(_gestureState != PANNING) {
				_clickedCoords = _coords;
				_gestureState = PANNING;
			}
		} else {
			_gestureState = NONE;
		}
	}

	switch(_gestureState) {
	case ROTATING:
		{
			sl::vector3f p,n;
			p = _cam.pos-_cam.target;
			float phi = (_cam.coordSyst == Y_UP) ? atan2(p[2], p[0]) : atan2(p[1], p[0]);
			float theta = (_cam.coordSyst == Y_UP) ? acos(p[1] / _rotRadius) : acos(p[2] / _rotRadius);
			float dPhi = (_coords[0]-_clickedCoords[0])/_winWidth;
			float dTheta = (_clickedCoords[1]-_coords[1])/_winHeight;
			if (_cam.coordSyst != Y_UP) dPhi = -dPhi;
			phi += dPhi * _rotFactor;
			if((theta>0.1 && dTheta<0) || (theta<M_PI-0.1f && dTheta>0))
				theta += dTheta * _rotFactor;

			if (_cam.coordSyst == Y_UP) {
				n[0] = _rotRadius * sin(theta) * cos(phi);
				n[1] = _rotRadius * cos(theta);
				n[2] = _rotRadius * sin(theta) * sin(phi);
			} else {
				n[0] = _rotRadius * sin(theta) * cos(phi);
				n[1] = _rotRadius * sin(theta) * sin(phi);
				n[2] = _rotRadius * cos(theta);
			}
			if (_camController == ORBIT) {
				_cam.pos += n-p;
			} else {
				sl::vector3f dir = (_cam.target - _cam.pos);

				float deltaX = 2.0f * 3.14f / 180.0f * dPhi * _rotFactor * -100.0f;
				float deltaY = 2.0f * 3.14f / 180.0f * dTheta * _rotFactor * -100.0f;

				if (_cam.coordSyst == Z_UP) {
					deltaX *= -1;
				}

				sl::rigid_body_map3f rX = sl::linear_map_factory3f::rotation(_cam.up, deltaX);

				sl::vector3f camForward = dir.ok_normalized();
				sl::vector3f camRight = _cam.up.cross(camForward);
				sl::rigid_body_map3f rY = sl::linear_map_factory3f::rotation(camRight, deltaY);
				
				_cam.target = _cam.pos + rX * rY * dir;
				updateMatrices();
			}
		}
		break;
	case ZOOMING:
		_cam.pos += float((_coords[1]-_clickedCoords[1])/_winHeight) * (_cam.target-_cam.pos) * _zoomFactor;
		_rotRadius = (_cam.target-_cam.pos).two_norm();
		break;
		
	case PANNING:
		sl::vector3f d(0,0,0);
		sl::vector3f camDir = (_cam.target-_cam.pos).ok_normalized();
		sl::vector3f right = camDir.cross(_cam.up);
		sl::vector3f actualUp = right.cross(camDir);
		d += right * (_clickedCoords[0]-_coords[0])/_winWidth;
		d += actualUp * (_coords[1] - _clickedCoords[1])/_winHeight;
		d *= _panFactor * (_cam.target - _cam.pos).two_norm();
		_cam.pos += d;
		_cam.target += d;
		break;
	}

	if(_gestureState!=NONE) updateMatrices();
	if(sticky) _clickedCoords = _coords;

	if (_state == WT_RECORDING) _walkthrough.keys.push_back(WTKey(_cam.pos, _cam.target));
}

void Camera::setInitCamera(sl::point3f pos, sl::point3f target, TCoordinateSystem coordSyst, float fov, float nearClip, float farClip) {
	_initCam.pos = pos;
	_initCam.target = target;
	_initCam.coordSyst = coordSyst;
	_initCam.fov = fov;
	_initCam.nearClip = nearClip;
	_initCam.farClip = farClip;
	if (coordSyst == Z_UP)
		_initCam.up = sl::vector3f(0, 0, 1);
	else
		_initCam.up = sl::vector3f(0, 1, 0);
	reset();
}


void Camera::updateMatrices() {
	_viewMatrix = sl::linear_map_factory3f::lookat(_cam.pos, _cam.target, _cam.up);
	_projMatrix = sl::linear_map_factory3f::perspective(float(M_PI/180.0f)*_cam.fov, _winWidth/_winHeight, _cam.nearClip, _cam.farClip);
	_viewMatrixInv = _viewMatrix.inverse();
	_projMatrixInv = _projMatrix.inverse();
}

void Camera::reset() {
	_cam = _initCam;
	_rotRadius = (_cam.target - _cam.pos).two_norm();
	updateMatrices();
}

void Camera::recordWalkthrough(std::string filename) {
	if (_state == WT_RECORDING) 
		printf("* Camera: RECORDING (resetting current record)\n");
	else
		printf("* Camera: RECORDING\n");
	_state = WT_RECORDING;
	_walkthrough.reset();

	if (filename.empty()) {
		time_t t = time(0);
		struct tm * now = localtime(&t);
		_walkthrough.filename = "walkthrough_"
			+ std::to_string(now->tm_year + 1900) + "_"
			+ std::to_string(now->tm_mon + 1) + "_"
			+ std::to_string(now->tm_mday) + "_"
			+ std::to_string(now->tm_hour) + "_"
			+ std::to_string(now->tm_min) + "_"
			+ std::to_string(now->tm_sec) + ".txt";
	}
	else {
		_walkthrough.filename = filename;
	}
}


bool Camera::playWalkthrough(int interpol, bool rewind) {
	if (_walkthrough.keys.empty()) return false;
	printf("* Camera: PLAYING\n");
	_state = WT_PLAYING;
	_wtInterpol = interpol;
	if(rewind) _walkthrough.currentKey = 0;
	return true;
}

void Camera::stopWalkthrough() {
	printf("* Camera: STOP");
	if (_state == WT_RECORDING) {
		printf(" (Saving '%s'... ", _walkthrough.filename.c_str());
		FILE * f = fopen(_walkthrough.filename.c_str(), "w");
		if(f == NULL) {
			std::cout << "WARNING: Camera::stopWalkthrough() Can't open file " << _walkthrough.filename << std::endl;
			_state = FREE;
			return;
		}
		fprintf(f, "%i %f %f %f\n", _cam.coordSyst, _cam.fov, _cam.nearClip, _cam.farClip);
		for (auto i : _walkthrough.keys)
			fprintf(f, "%f %f %f %f %f %f\n", i.pos[0], i.pos[1], i.pos[2], i.target[0], i.target[1], i.target[2]);
		fclose(f);
		printf("OK!\n");
	}
	_walkthrough.currentKey = 0;
	_state = FREE;
	printf("\n");
}

void Camera::printControls() {
	printf(
		"* Camera Controls: \n"
		"\t- Mouse Left Button   : orbit\n"
		"\t- Mouse Right Button  : zoom\n"
		"\t- Mouse Middle Button : pam\n"
		"\t- W,S,A,D : forward, backward, lateral\n"
		"\t- Q,E     : rotate over UP vector\n"
		"\t- U       : switch UP vector betweem Y and Z axes\n"
		"\t- SPACE   : reset camera\n"
		"\t- C       : print camera info\n"
	);
}


void Camera::loadWalkthrough(std::string filename) {
	std::ifstream ifs;
	ifs.open(filename);
	printf("* Camera: Loading '%s'... ", filename.c_str());
	if (!ifs.good()) {
		std::cout << "WARNING: Camera::loadWalkthrough() Can't open file " << _walkthrough.filename << std::endl;
		return;
	}

	_walkthrough.reset();
	std::string line, tmp;

	do { std::getline(ifs, line); } while (!ifs.eof() && line.empty());
	sscanf(line.c_str(), "%i %f %f %f", &_cam.coordSyst, &_cam.fov, &_cam.nearClip, &_cam.farClip);
	

	if (_cam.coordSyst == Z_UP) {
		printf("Z UP\n");
		_cam.up = sl::vector3f(0, 0, 1);
	}
	else {
		printf("Y UP\n");
		_cam.up = sl::vector3f(0, 1, 0);
	}

	while (!ifs.eof()) {
		
		std::getline(ifs, line);
		if (line.empty()) continue;
		sl::point3f pos, target;
		if (sscanf(line.c_str(), "%f %f %f %f %f %f", &pos[0], &pos[1], &pos[2], &target[0], &target[1], &target[2]) == 6) {
			_walkthrough.keys.push_back(WTKey(pos, target));
		}
		else {
			std::cout << "Camera:loadWalkthrough: WARNING: In file " << filename << " couldn't parse line \"" << line << "\"" << std::endl;
		}
	}
	ifs.close();
	printf("OK!\n");
}
