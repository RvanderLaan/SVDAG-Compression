
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

/// TRIANGLE - BOX Intersection Code ------------------------------------
/// by Tomas Akenine-M�ller
/// see http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/
/////////////////////////////////////////////////////////////////////////

#include <symvox/test_triangle_box.hpp>

#define X 0
#define Y 1
#define Z 2

#define FINDMINMAX(x0,x1,x2,min,max) \
min = max = x0;			\
if (x1<min) min = x1;	\
if (x1>max) max = x1;	\
if (x2<min) min = x2;	\
if (x2>max) max = x2;

bool planeBoxOverlap(const sl::vector3d normal, const sl::point3d vert, const double maxbox) {
	int q;
	sl::vector3d vmin, vmax;
	double v;
	for (q = X; q <= Z; q++) {
		v = vert[q];
		if (normal[q] > 0.0f) {
			vmin[q] = -maxbox - v;
			vmax[q] = maxbox - v;
		}
		else {
			vmin[q] = maxbox - v;
			vmax[q] = -maxbox - v;
		}
	}
	if (normal.dot(vmin) > 0.0) return false;
	if (normal.dot(vmax) >= 0.0) return true;
	return false;
}


/*========================== X-tests ==========================*/
#define AXISTEST_X01(a, b, fa, fb)								\
p0 = a*v0[Y] - b*v0[Z];											\
p2 = a*v2[Y] - b*v2[Z];											\
if (p0<p2) { min = p0; max = p2; } else { min = p2; max = p0; }	\
rad = (fa + fb) * boxHalfSide;									\
if (min>rad || max<-rad) return false;

#define AXISTEST_X2(a, b, fa, fb)								\
p0 = a*v0[Y] - b*v0[Z];											\
p1 = a*v1[Y] - b*v1[Z];											\
if (p0<p1) { min = p0; max = p1; } else { min = p1; max = p0; }	\
rad = (fa + fb) * boxHalfSide;									\
if (min>rad || max<-rad) return false;

/*========================== Y-tests ==========================*/
#define AXISTEST_Y02(a, b, fa, fb)								\
p0 = -a*v0[X] + b*v0[Z];										\
p2 = -a*v2[X] + b*v2[Z];										\
if (p0<p2) { min = p0; max = p2; } else { min = p2; max = p0; }	\
rad = (fa + fb) * boxHalfSide;									\
if (min>rad || max<-rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)								\
p0 = -a*v0[X] + b*v0[Z];										\
p1 = -a*v1[X] + b*v1[Z];										\
if (p0<p1) { min = p0; max = p1; } else { min = p1; max = p0; }	\
rad = (fa + fb) * boxHalfSide;									\
if (min>rad || max<-rad) return false;

/*========================== Z-tests ==========================*/
#define AXISTEST_Z12(a, b, fa, fb)								\
p1 = a*v1[X] - b*v1[Y];											\
p2 = a*v2[X] - b*v2[Y];											\
if (p2<p1) { min = p2; max = p1; } else { min = p1; max = p2; }	\
rad = (fa + fb) * boxHalfSide;									\
if (min>rad || max<-rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)								\
p0 = a*v0[X] - b*v0[Y];											\
p1 = a*v1[X] - b*v1[Y];											\
if (p0<p1) { min = p0; max = p1; } else { min = p1; max = p0; }	\
rad = (fa + fb) * boxHalfSide;									\
if (min>rad || max<-rad) return false;


bool testTriBox(const sl::point3d boxCenter, const double boxHalfSide, const sl::point3f *tri) {

	/*    Use separating axis theorem to test overlap between triangle and box */
	/*    need to test for overlap in these directions: */
	/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
	/*       we do not even need to test these) */
	/*    2) normal of the triangle */
	/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
	/*       this gives 3x3=9 more tests */

	sl::point3d v0, v1, v2;
	double min, max, p0, p1, p2, rad, fex, fey, fez;
	sl::vector3d normal, e0, e1, e2;

	sl::point3d triD[3];
	triD[0] = sl::point3d(tri[0][0], tri[0][1], tri[0][2]);
	triD[1] = sl::point3d(tri[1][0], tri[1][1], tri[1][2]);
	triD[2] = sl::point3d(tri[2][0], tri[2][1], tri[2][2]);

	/* This is the fastest branch on Sun */
	/* move everything so that the boxcenter is in (0,0,0) */
	v0 = triD[0] - boxCenter.as_vector();
	v1 = triD[1] - boxCenter.as_vector();
	v2 = triD[2] - boxCenter.as_vector();

	/* compute triangle edges */
	e0 = v1 - v0;      /* tri edge 0 */
	e1 = v2 - v1;      /* tri edge 1 */
	e2 = v0 - v2;      /* tri edge 2 */

	/* Bullet 3:  */
	/*  test the 9 tests first (this was faster) */
	fex = fabs(e0[X]);
	fey = fabs(e0[Y]);
	fez = fabs(e0[Z]);
	AXISTEST_X01(e0[Z], e0[Y], fez, fey);
	AXISTEST_Y02(e0[Z], e0[X], fez, fex);
	AXISTEST_Z12(e0[Y], e0[X], fey, fex);

	fex = fabs(e1[X]);
	fey = fabs(e1[Y]);
	fez = fabs(e1[Z]);
	AXISTEST_X01(e1[Z], e1[Y], fez, fey);
	AXISTEST_Y02(e1[Z], e1[X], fez, fex);
	AXISTEST_Z0(e1[Y], e1[X], fey, fex);

	fex = fabs(e2[X]);
	fey = fabs(e2[Y]);
	fez = fabs(e2[Z]);
	AXISTEST_X2(e2[Z], e2[Y], fez, fey);
	AXISTEST_Y1(e2[Z], e2[X], fez, fex);
	AXISTEST_Z12(e2[Y], e2[X], fey, fex);

	/* Bullet 1: */
	/*  first test overlap in the {x,y,z}-directions */
	/*  find min, max of the triangle each direction, and test for overlap in */
	/*  that direction -- this is equivalent to testing a minimal AABB around */
	/*  the triangle against the AABB */

	/* test in X-direction */
	FINDMINMAX(v0[X], v1[X], v2[X], min, max);
	if (min > boxHalfSide || max < -boxHalfSide) return false;

	/* test in Y-direction */
	FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);
	if (min > boxHalfSide || max < -boxHalfSide) return false;

	/* test in Z-direction */
	FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);
	if (min > boxHalfSide || max < -boxHalfSide) return false;

	/* Bullet 2: */
	/*  test if the box intersects the plane of the triangle */
	/*  compute plane equation of triangle: normal*x+d=0 */
	normal = e0.cross(e1);
	if (!planeBoxOverlap(normal, v0, boxHalfSide)) return false;

	return true;   /* box and triangle overlaps */

}

float clamp(float n, float lower, float upper) {
    return std::max(lower, std::min(n, upper));
}

/** From https://www.gamedev.net/forums/topic/552906-closest-point-on-triangle/ **/
sl::point3f closestPointOnTri(const sl::point3d p, const sl::point3f *triangle) {
	//sl::vector3f pos = sl::conv_to<sl::vector3f>::from(p.as_vector());
	sl::vector3f pos = sl::vector3f(p[0], p[1], p[2]);

    sl::vector3f edge0 = triangle[1] - triangle[0];
    sl::vector3f edge1 = triangle[2] - triangle[0];
    sl::vector3f v0 = triangle[0].as_vector() - pos;

    float a = edge0.dot( edge0 );
    float b = edge0.dot( edge1 );
    float c = edge1.dot( edge1 );
    float d = edge0.dot( v0 );
    float e = edge1.dot( v0 );

    float det = a*c - b*b;
    float s = b*e - c*d;
    float t = b*d - a*e;

    if ( s + t < det ) {
        if ( s < 0.f ) {
            if ( t < 0.f ) {
                if ( d < 0.f ) {
                    s = clamp( -d/a, 0.f, 1.f );
                    t = 0.f;
                } else {
                    s = 0.f;
                    t = clamp( -e/c, 0.f, 1.f );
                }
            } else {
                s = 0.f;
                t = clamp( -e/c, 0.f, 1.f );
            }
        } else if ( t < 0.f ) {
            s = clamp( -d/a, 0.f, 1.f );
            t = 0.f;
        } else {
            float invDet = 1.f / det;
            s *= invDet;
            t *= invDet;
        }
    } else {
        if ( s < 0.f ) {
            float tmp0 = b+d;
            float tmp1 = c+e;
            if ( tmp1 > tmp0 ) {
                float numer = tmp1 - tmp0;
                float denom = a-2*b+c;
                s = clamp( numer/denom, 0.f, 1.f );
                t = 1-s;
            } else {
                t = clamp( -e/c, 0.f, 1.f );
                s = 0.f;
            }
        }
        else if ( t < 0.f ) {
            if ( a+d > b+e ) {
                float numer = c+e-b-d;
                float denom = a-2*b+c;
                s = clamp( numer/denom, 0.f, 1.f );
                t = 1-s;
            } else {
                s = clamp( -e/c, 0.f, 1.f );
                t = 0.f;
            }
        } else {
            float numer = c+e-b-d;
            float denom = a-2*b+c;
            s = clamp( numer/denom, 0.f, 1.f );
            t = 1.f - s;
        }
    }

    return triangle[0] + s * edge0 + t * edge1;
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void barycentric(sl::point3d p, const sl::point3f *tri, float &u, float &v, float &w)
{
    sl::point3d a, b, c;
    //a = sl::conv_to< sl::point3d >::from(tri[0]);
    //b = sl::conv_to< sl::point3d >::from(tri[1]);
    //c = sl::conv_to< sl::point3d >::from(tri[2]);

	a = sl::point3d(tri[0][0], tri[0][1], tri[0][2]);
	b = sl::point3d(tri[1][0], tri[1][1], tri[1][2]);
	c = sl::point3d(tri[2][0], tri[2][1], tri[2][2]);

    sl::vector3d v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}
