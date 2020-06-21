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


#version 440 core

// Defined from outside:
// SVDAG || USSVDAG || SSVDAG
// SSVDAG_TEX3D
// VIEWER_MODE || DEPTH_MODE || SHADOW_MODE || AO_MODE
// N_HS_SAMPLES

uniform mat4 viewMatInv;
uniform mat4 projMatInv;
uniform vec2 screenRes;

uniform vec3 sceneBBoxMin;
uniform vec3 sceneBBoxMax;
uniform vec3 sceneCenter;
uniform float rootHalfSide;

uniform uint maxIters;
uniform uint drawLevel;
uniform float projectionFactor;

#define TEST_CHILD_INDIR_TEX 1

#if (SVDAG || USSVDAG)
uniform usamplerBuffer nodes;
#elif (SSVDAG)
#if SSVDAG_TEX3D
uniform usampler3D innerNodes;
#else
uniform usamplerBuffer innerNodes;
#endif
uniform usamplerBuffer leafNodes;
uniform uint levelOffsets[INNER_LEVELS]; // INNER_LEVELS defined from app

#if TEST_CHILD_INDIR_TEX
uniform isamplerBuffer childIndir;
#endif 

#endif

#if VIEWER_MODE
uniform int viewerRenderMode;
uniform uint selectedVoxelIndex;
uniform bool randomColors;
uniform vec3 lightPos;
uniform bool enableShadows;
uniform float normalEpsilon;
uniform int normSamples;
#elif SHADOW_MODE
uniform sampler2D hitPosTex;
uniform sampler2D hitNormTex;
uniform sampler2D depthTex;
uniform vec3 lightPos;
#elif AO_MODE
uniform sampler2D hitPosTex;
uniform sampler2D hitNormTex;
uniform sampler2D depthTex;
uniform vec3 hsSamples[N_HS_SAMPLES];
uniform int numAORays;
uniform float lengthAORays;
#endif

#if (VIEWER_MODE && (SSVDAG || USSVDAG))
uniform bool freqColors;
#endif

#if (VIEWER_MODE || DEPTH_MODE)
uniform sampler2D minDepthTex;
uniform bool useMinDepthTex;
#endif

#if VIEWER_MODE
layout(location = 0) out vec3 color;
#elif DEPTH_MODE
layout(location = 0) out vec3 output_t;
#elif (SHADOW_MODE || AO_MODE)
layout(location = 0) out float output_t;
#endif

// TODO: Shaders probably work better with floats isntead of ints
//https://stackoverflow.com/questions/1737726/how-to-perform-rgb-yuv-conversion-in-c-c
#define CLIP(X) ( (X) > 255 ? 255 : (X) < 0 ? 0 : X)
#define CYCbCr2R(Y, Cb, Cr) CLIP( Y + ( 91881 * Cr >> 16 ) - 179 )
#define CYCbCr2G(Y, Cb, Cr) CLIP( Y - (( 22544 * Cb + 46793 * Cr ) >> 16) + 135)
#define CYCbCr2B(Y, Cb, Cr) CLIP( Y + (116129 * Cb >> 16 ) - 226 )

struct Ray {
	vec3 o;
	vec3 d;
};

struct traversal_status {
	float t_current;
	int node_index;
	uint hdr;
	ivec3 mirror_mask;
	uvec2 leaf_data;
	
	ivec3 idx; 
	ivec3 local_idx;  uint child_linear_index; // Could be merged
	vec3 t_next_crossing; // exit of current voxel on x,y,z
	vec3 inv_ray_d;
	ivec3 delta_idx;
	int current_node_size;
	
	float cell_size;
	uint level;

	vec3 attr_sum;
};

//////////////////////////// STACK STUFF ////////////////////

#define MAX_STACK_DEPTH (INNER_LEVELS+1)

ivec3 stack[MAX_STACK_DEPTH]; // stack of ivec3(nodeIndex, hdr, mask) for each traversal level 
uint stack_size = 0;

vec3 get_header_attrs(uint hdr) {
	vec3 res = vec3(0);
	res.x = ((hdr & 0XFF000000) >> 24) / 255.0;
	res.y = ((hdr & 0X00FF0000) >> 16) / 255.0;
	res.z = ((hdr & 0X0000FF00) >>  8) / 255.0;
	return res;
}

void stack_push(in const int node, in const uint hdr, in const ivec3 mirror_mask, in const uint level) {
	const int mask = mirror_mask.x | (mirror_mask.y << 1) | (mirror_mask.z << 2) | int(level << 3);
	stack[stack_size] = ivec3(node, hdr, mask);
	++stack_size;
}

void stack_pop_in(out int node, out uint hdr, out ivec3 mirror_mask, out uint level) {
	--stack_size;
	const ivec3 node_mask = stack[stack_size];
	node = node_mask.x;
	hdr = node_mask.y;
	const int mask = node_mask.z;
	mirror_mask = ivec3(mask & 1, (mask >> 1) & 1, (mask >> 2) & 1);
	level = uint(mask >> 3) & 255;
}

bool stack_is_empty() {
	return stack_size == 0;
}

uint voxel_to_linear_idx(in const ivec3 mirror_mask, in ivec3 idx, in const int sz) {
	idx = (ivec3(1) - 2 * mirror_mask) * idx + mirror_mask * (sz-1);
	return uint(idx.z + sz * (idx.y + sz * idx.x));
}

#if SVDAG // ----------------------------------------------------------------------

#define LEAF_SIZE 2

uint myFetch(in const int idx) {
	// On AMD GPUs this causes some problems for files > 32 MM. Took some time to find this out
	// The idx is signed, so cast it to unsigned
	// const uvec4 tmp = texelFetch(nodes, idx/4);
	const uvec4 tmp = texelFetch(nodes, int(uint(idx)/4));

	const int selected = idx%4;
	uint result;
	if      (selected == 0) result = tmp.x;
	else if (selected == 1) result = tmp.y;
	else if (selected == 2) result = tmp.z;
	else if (selected == 3) result = tmp.w;
	return result;
}

bool fetch_voxel_bit(in const traversal_status ts) {
	return (ts.hdr & (1 << ts.child_linear_index)) != 0;
}

void fetch_data(inout traversal_status ts) {
	ts.hdr = myFetch(ts.node_index);
}

void fetch_child_index_in(inout traversal_status ts) {
	
	const int childPtrPos = bitCount((ts.hdr & 0xFF) >> ts.child_linear_index);
	ts.node_index = int(myFetch(ts.node_index + childPtrPos));
}

#elif USSVDAG // ------------------------------------------------------------------

#define LEAF_SIZE 2

uint myFetch(in const int idx) {
	const uvec4 tmp = texelFetch(nodes, int(uint(idx)/4));
	const int selected = idx%4;
	uint result;
	if      (selected == 0) result = tmp.x;
	else if (selected == 1) result = tmp.y;
	else if (selected == 2) result = tmp.z;
	else if (selected == 3) result = tmp.w;
	return result;
}

bool fetch_voxel_bit(in const traversal_status ts) {
	return (ts.hdr & (1 << ts.child_linear_index)) != 0;
}

void fetch_data(inout traversal_status ts) {
	ts.hdr = myFetch(ts.node_index);
}

void fetch_child_index_in(inout traversal_status ts) {
	
	const int childPtrPos = bitCount((ts.hdr & 0xFF) >> ts.child_linear_index);
	ts.node_index = int(myFetch(ts.node_index + childPtrPos));
	ivec3 m = ivec3(0,0,0); // FIXME optimize
	if((ts.hdr & (1 << (ts.child_linear_index +  8))) != 0) m.x = 1;
	if((ts.hdr & (1 << (ts.child_linear_index + 16))) != 0) m.y = 1;
	if((ts.hdr & (1 << (ts.child_linear_index + 24))) != 0) m.z = 1;
	ts.mirror_mask ^= m;
}

#elif SSVDAG // -----------------------------------------------------------------------

#define LEAF_SIZE 4

#if SSVDAG_TEX3D

ivec2 myInnerFetch(in const int idx) {
	ivec2 result;
	result.x = int(texelFetch( innerNodes, 
		ivec3(
			idx % TEX3D_SIZE,
			(idx/TEX3D_SIZE) % TEX3D_SIZE,
			idx/(TEX3D_SIZE_POW2)
			),
		0).x);
	const int idx2 = idx + 1;
	result.y = int(texelFetch( innerNodes, 
		ivec3(
			idx2 % TEX3D_SIZE,
			(idx2/TEX3D_SIZE) % TEX3D_SIZE,
			idx2/(TEX3D_SIZE_POW2)
			),
		0).x);
	
	return result;
}

uint myInnerHeaderFetch(in const int idx) {
	return texelFetch( innerNodes, 
		ivec3(
			idx % TEX3D_SIZE,
			(idx/TEX3D_SIZE) % TEX3D_SIZE,
			idx/(TEX3D_SIZE_POW2)
			),
		0).x;
}
#else
ivec2 myInnerFetch(in const int idx) {
	const int realIdx = int(uint(idx)/4);
	const ivec4 tmp = ivec4(texelFetch(innerNodes, realIdx));
	const int selected = idx - realIdx * 4;
	ivec2 result;
	if      (selected == 0) result = tmp.xy;
	else if (selected == 1) result = tmp.yz;
	else if (selected == 2) result = tmp.zw;
	else if (selected == 3) {
		result.x = tmp.w;
		result.y = int(texelFetch(innerNodes, realIdx+1).x);
	}
	return result;
}

uint myInnerHeaderFetch(in const int idx) {
	const int realIdx = int(uint(idx)/4);
	const uvec4 tmp = texelFetch(innerNodes, realIdx);
	const int selected = idx - realIdx * 4;
	uint result;
	if      (selected == 0) result = tmp.x;
	else if (selected == 1) result = tmp.y;
	else if (selected == 2) result = tmp.z;
	else if (selected == 3) result = tmp.w;
	return result;
}

#endif
int get_child_mask(in const uint hdr, in const uint child_id) {
	return int((hdr >> (child_id<<1)) & 3U);
}

bool leaf_fetch_voxel_bit(in const traversal_status ts) {
	const uint word  = (ts.child_linear_index & 32);
	const uint bit_mask = 1 << (ts.child_linear_index & 31);
	const uint value = (word == 0 ? ts.leaf_data.x : ts.leaf_data.y) & bit_mask;
	return value != 0;
}

bool fetch_voxel_bit(in const traversal_status ts) {
	// for inner nodes, true if at least 1 of the 2 corresponding mask bits != 0
	// For leafs if corresponding voxel bit is not empty
	return (ts.level < INNER_LEVELS ?
		(get_child_mask(ts.hdr, ts.child_linear_index) != 0) :
		leaf_fetch_voxel_bit(ts));
}

void fetch_data(inout traversal_status ts) {
	if (ts.level < INNER_LEVELS)
		ts.hdr = myInnerHeaderFetch(ts.node_index + int(levelOffsets[ts.level]));
	else
		ts.leaf_data = texelFetch(leafNodes, ts.node_index).xy;
}

void fetch_child_index_in(inout traversal_status ts) {
	// Compute offset to get child pointer information
	int node_offset = 1 + int(levelOffsets[ts.level]) + ts.node_index;
#if TEST_CHILD_INDIR_TEX
	node_offset += texelFetch(childIndir, int(ts.hdr + 65536 * ts.child_linear_index)).x;
#else
	for (uint i = 7; i > ts.child_linear_index; --i) {
		node_offset += min(get_child_mask(ts.hdr, i), 2);
	}
#endif
	// Read child pointer and child mask
	//int child_index = int(texelFetch(innerNodes, node_offset).x);
	ivec2 child_index = myInnerFetch(node_offset);
	const int child_mask = get_child_mask(ts.hdr, ts.child_linear_index);
	
	// read mirror info and remove this info from childidx
	ivec3 mirror;
	mirror.x = ((child_index.x >> 13) & 1);
	mirror.y = ((child_index.x >> 14) & 1);
	mirror.z = ((child_index.x >> 15) & 1);
	child_index.x &= -57345; // ~(7 << 13);
	
	if (child_mask > 1) {
		//int tmp = int(texelFetch(innerNodes, node_offset+1).x);
		//const int tmp = int(myInnerFetch(node_offset+1));
		child_index.x = (child_index.x << 16) | child_index.y;
		const int x_bit = child_mask & 1;
		child_index.x |= (x_bit << 29);
	}
	
	// Update accumulated mirror mask, node and its header from next level
	ts.mirror_mask ^= mirror;
	ts.node_index = child_index.x;
}

#endif

///////////////////////////// DDA PRIMITIVES

bool in_bounds(in const ivec3 local_idx, in const int sz) { 
	const bvec3 cond0 = lessThan(local_idx, ivec3(sz));
	const bvec3 cond1 = lessThanEqual(ivec3(0,0,0), local_idx);
	return cond0.x && cond0.y && cond0.z && cond1.x && cond1.y && cond1.z;
}

vec2 intersectAABB(in const Ray r, in const vec3 aabbMin, in const vec3 aabbMax) {
  const vec3 t1 = (aabbMin - r.o)/r.d;
  const vec3 t2 = (aabbMax - r.o)/r.d;
  const vec3 tMin = min(t1, t2);
  const vec3 tMax = max(t1, t2);

  vec2 t = vec2(max(max(tMin.x, 0.0), max(tMin.y, tMin.z)), min(tMax.x, min(tMax.y, tMax.z)));
  return t;
}

bool resolution_ok(float t, float cell_size, float projection_factor) {
  return (cell_size * projection_factor) < t;
}

// ==========================================================================

void dda_init(in const Ray r, inout traversal_status ts) {
	// Init dda FIXME USE OCTREE POINT LOCATION
	const float voxel_eps = 1.0f/(256.*1024.);
	const vec3 p_a = r.o + (ts.t_current + voxel_eps) * r.d; // find current pos
	ts.idx = ivec3(p_a / ts.cell_size); // current global grid voxel
	
	// During initialization do not step back for dir < 0, because it would move of more than once cell
	const ivec3 delta_idx_conservative = max(ivec3(0), ts.delta_idx);
	const ivec3 idx_next = ts.idx + delta_idx_conservative;
	const vec3 p_next_a = idx_next * ts.cell_size;	// this should be the plane
	
	ts.t_next_crossing = (p_next_a - r.o) * ts.inv_ray_d;
	ts.local_idx = ts.idx % 2;
}

// https://stackoverflow.com/questions/24599502/is-there-a-built-in-function-in-glsl-for-and-or-is-there-some-optimized-method-f
bvec3 bvec3_and(const bvec3 one, const bvec3 two) {
	return bvec3(uvec3(one) & uvec3(two));
}

// Returns the direction of the step
ivec3 dda_next(inout traversal_status ts) {
	const bvec3 b1 = lessThan(ts.t_next_crossing.xyz, ts.t_next_crossing.yzx);
	const bvec3 b2 = lessThanEqual(ts.t_next_crossing.xyz, ts.t_next_crossing.zxy);
	const bvec3 mask = bvec3_and(b1, b2);
	const vec3 mask_v3 = vec3(mask); 			
	
	//All components of mask are false except the one components to the shortest t_next_crossing
	// which is the direction in which the step have to be done
	const ivec3 delta = ivec3(mask) * ts.delta_idx;
	ts.idx += delta;
	ts.local_idx += delta;
	
	ts.t_current = dot(mask_v3, ts.t_next_crossing);
	ts.t_next_crossing += mask_v3 * ts.cell_size * abs(ts.inv_ray_d);
	return delta;
}
 	     
ivec3 dda_next_delta_index(in const traversal_status ts) {
	const bvec3 b1 = lessThan(ts.t_next_crossing.xyz, ts.t_next_crossing.yzx);
	const bvec3 b2 = lessThanEqual(ts.t_next_crossing.xyz, ts.t_next_crossing.zxy);
	const bvec3 mask = bvec3_and(b1, b2);
	return ivec3(mask) * ts.delta_idx;
}

void up_in(in const Ray r, inout traversal_status ts) {
	uint delta_level = ts.level;
	stack_pop_in(ts.node_index, ts.hdr, ts.mirror_mask, ts.level);
	delta_level -= ts.level;
	
	ts.idx >>= delta_level; // always delta_level >= 1
	ts.cell_size *= (1 << delta_level);  
	ts.current_node_size = ts.level < INNER_LEVELS ? 2 : LEAF_SIZE;
	ts.local_idx = ts.idx & 1; 
	
	const ivec3 delta_idx_conservative = max(ivec3(0), ts.delta_idx);
	const ivec3 idx_next = ts.idx + delta_idx_conservative;
	const vec3 p_next_a = idx_next * ts.cell_size;	// this should be the plane
	ts.t_next_crossing = (p_next_a - r.o) * ts.inv_ray_d;
}

void go_down_one_level(in const Ray r, inout traversal_status ts) {
	++ts.level;
	ts.cell_size*=0.5;
	
	// Init ts idx, t_next_crossing, local_idx using octree point location
	const vec3 p_a = r.o + ts.t_current * r.d;		
	const vec3 p_center = (ts.idx * 2 + 1) * ts.cell_size;
	const bvec3 child_pos = lessThan(p_center, p_a);
	const ivec3 delta = ivec3(child_pos);
	ts.idx = ts.idx*2 + delta;
	
	const ivec3 delta_idx_conservative = max(ivec3(0), ts.delta_idx);
	const ivec3 idx_next = ts.idx + delta_idx_conservative;
	const vec3 p_next_a = idx_next * ts.cell_size;	// this should be the plane
	
	ts.t_next_crossing = (p_next_a - r.o) * ts.inv_ray_d;
	ts.local_idx = ts.idx & (ts.current_node_size-1);
}

void down_in(in const Ray r, inout traversal_status ts) {
	// Check/push next
	const ivec3 local_idx = ts.local_idx;
	const ivec3 delta = dda_next_delta_index(ts);    
	
	if (in_bounds(local_idx + delta, 2)) { 
		stack_push(ts.node_index, ts.hdr, ts.mirror_mask, ts.level);
		
	}
	  
	// Go down to next level: Fetch child index (store in node_idx)
	// and update accumulated_mirror_mask
	fetch_child_index_in(ts);
	
	go_down_one_level(r, ts);
	
	// Add to color sum
//	ts.attr_sum += get_header_attrs(ts.hdr);

	
	if (ts.level == INNER_LEVELS) {

		// GO TO LEAVES
		ts.current_node_size = LEAF_SIZE;
		int voxel_count = LEAF_SIZE / 2;
		while (voxel_count > 1) {
			go_down_one_level(r, ts);
			voxel_count >>= 1;
		}
	}
}


bool transform_ray(inout Ray r, inout vec2 t_min_max)  {
	const float epsilon = 1E-4f;
	const vec3 sign_rd = sign(r.d);
	
	// Move ray to LOCAL box
	const float scale = 1.0/(2.0 * rootHalfSide);
	const vec3 octree_min = sceneCenter - vec3(rootHalfSide);
	const vec3 octree_max = sceneCenter + vec3(rootHalfSide);
	r.o = r.o - octree_min;
	r.o *= scale;
	t_min_max *= scale;
	
	// avoid div by zero
	if (r.d.x * sign_rd.x < epsilon) r.d.x = sign_rd.x * epsilon;
	if (r.d.y * sign_rd.y < epsilon) r.d.y = sign_rd.y * epsilon;
	if (r.d.z * sign_rd.z < epsilon) r.d.z = sign_rd.z * epsilon;
	
	const vec3 clip_box_min = (sceneBBoxMin - octree_min) * scale; 
	const vec3 clip_box_max = (sceneBBoxMax - octree_min) * scale; 
	
	const vec2 t_intersection = intersectAABB(r, clip_box_min, clip_box_max);
	
	t_min_max.x = max(t_intersection.x, t_min_max.x + 1e-10);
	t_min_max.y = min(t_intersection.y, t_min_max.y);
	
	return t_intersection.x < t_intersection.y;
}


void init(inout Ray r, inout traversal_status ts) {
	ts.inv_ray_d = vec3(1.0/r.d);
	ts.delta_idx = ivec3(sign(r.d));
	
	// Level status
	ts.mirror_mask = ivec3(0,0,0);
	ts.level = 0;
	ts.cell_size = 0.5;
	
	// Step status
	dda_init(r, ts);
	ts.current_node_size = 2;
	
	ts.node_index = 0;
	fetch_data(ts);
	ts.child_linear_index = voxel_to_linear_idx(ts.mirror_mask, ts.local_idx, ts.current_node_size);
	
	ts.attr_sum = vec3(0.5); // get_header_attrs(ts.hdr); // root node color
}

// TRACE RAY
// returns vec3
//	X: intersection t
//		 >= 0 := intersection!
//		-1   := inside scene bbox, but no intersection
//		-2   := -2 out of t bounds
//		-3   := too many iterations used (> maxIters)
//		-4   := out of scene bbox
//	Y: level of the intersection (-1 => no intersection)
//	Z: num Iterations used.
//  W: node index (-1 => no intersection)

vec4 trace_ray(in Ray r, in vec2 t_min_max, const in float projection_factor, out vec3 norm, out traversal_status ts_out) {
	
	if (!transform_ray(r, t_min_max))
		return vec4(-4.0,0,0,-1); // out of scene Bbox
	
	const float scale = 2.0 * rootHalfSide;
	traversal_status ts;
	ts.t_current = t_min_max.x;
	init(r, ts);

	ivec3 stepDir = ivec3(0);
	vec3 avgStepDir = vec3(0);

	// TODO: Each traversal iteration, re-use the attributes of the parent
	// sum up the attributes, divide by the nr of iterations at the end?
	// (later, the luma can have higher influence than chroma)
	
	int iteration_count = 0;
	const uint max_level = min(INNER_LEVELS, drawLevel-1);
	do {
		bool full_voxel = fetch_voxel_bit(ts);
	  
		if (!full_voxel) {
			stepDir = dda_next(ts);
			avgStepDir = avgStepDir * 0.5 + 0.5 * stepDir;
			if (!in_bounds(ts.local_idx, ts.current_node_size)) {
				if (stack_is_empty()) {
					return vec4(-1.,0,float(iteration_count),-1); // inside scene BBox, but no intersection
				}
				ts.attr_sum -= get_header_attrs(ts.hdr) - vec3(0.5);
				up_in(r, ts);
			}
		} else {
			const bool hit = (ts.level >= max_level || resolution_ok(ts.t_current, ts.cell_size, projection_factor)); 
			if (hit) {
				// Todo: Sample last few steps as norm: nope doesn't work
				norm = -vec3(stepDir);
				//norm = -avgStepDir;
				
				// Add to color sum
//				ts.attr_sum += get_header_attrs(ts.hdr);
				
				ts_out = ts;

				return vec4(ts.t_current * scale, ts.level, float(iteration_count), ts.node_index);  // intersection
			} else {
				down_in(r, ts);
				fetch_data(ts);
				if (fetch_voxel_bit(ts)) {
					ts.attr_sum += get_header_attrs(ts.hdr) - vec3(0.5);
				}
			}
		}
	    
		ts.child_linear_index = voxel_to_linear_idx(ts.mirror_mask, ts.local_idx, ts.current_node_size);
		++iteration_count;
	} while ((ts.t_current < t_min_max.y) && (iteration_count < maxIters));
	
	if (iteration_count >= maxIters) return vec4(-3.,0,float(iteration_count), -1); // too much itarations
	return vec4(-2.,0,float(iteration_count), -1); // intersection out of t bounds
}

vec3 fromHomog(in const vec4 v) {
	return v.xyz/v.w;
}

Ray computeCameraRay(in const vec2 pixelScreenCoords) {
	const vec4 pixel_s0 = vec4(pixelScreenCoords.x, pixelScreenCoords.y, 0, 1);
	const vec4 pixel_s1 = vec4(pixelScreenCoords.x, pixelScreenCoords.y, 1, 1);
	
	const vec3 pixel_w0 = fromHomog(projMatInv * pixel_s0);
	const vec3 pixel_w1 = fromHomog(projMatInv * pixel_s1);
	
	Ray r;
	r.o = vec3(0,0,0);
	r.d = normalize(pixel_w1 - pixel_w0);
	
	const vec3 o_prime = vec3(viewMatInv * vec4(r.o, 1));
	const vec3 e_prime = vec3(viewMatInv * vec4(r.d, 1));
	r.o = o_prime;
	r.d = normalize(e_prime - o_prime);
	return r;
}

#if (VIEWER_MODE || DEPTH_MODE)

float getMinT(in const int delta) {
	
	const ivec2 p = ivec2((gl_FragCoord.xy - delta / 2) / delta);

	float tl = texelFetch(minDepthTex, ivec2(p.x, p.y), 0).x;
	float tr = texelFetch(minDepthTex, ivec2(p.x+1, p.y), 0).x;
	float bl = texelFetch(minDepthTex, ivec2(p.x, p.y+1), 0).x;
	float br = texelFetch(minDepthTex, ivec2(p.x+1, p.y+1), 0).x;

	return min(min(tl, tr), min(bl, br));
}

#endif

///////////////////////// MAIN STUFF ////////////////////////////////////


#if VIEWER_MODE
// https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
uint hash( uint x ) {
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}
float floatConstruct( uint m ) {
    const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
    const uint ieeeOne      = 0x3F800000u; // 1.0 in IEEE binary32

    m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
    m |= ieeeOne;                          // Add fractional part to 1.0

    float  f = uintBitsToFloat( m );       // Range [1:2]
    return f - 1.0;                        // Range [0:1]
}
// Pseudo-random value in half-open range [0:1].
float randomFloat( uint x ) { return floatConstruct(hash(x)); }
//vec3 randomColor(uint x) {
//	return vec3(randomFloat(x), randomFloat(x / 2), randomFloat(x / 3));
//}

// From https://gist.github.com/mikhailov-work/0d177465a8151eb6ede1768d51d476c7
vec3 TurboColormap(in float x) {
	const vec4 kRedVec4 = vec4(-0.05195877, 5.18000081, -30.94853351, 81.96403246);
	const vec4 kGreenVec4 = vec4(0.16207513, 0.17712472, 15.24091500, -36.50657960);
	const vec4 kBlueVec4 = vec4(0.55305649, 3.00913185, -5.46192616, -11.11819092);
	const vec2 kRedVec2 = vec2(-86.53476570, 30.23299484);
	const vec2 kGreenVec2 = vec2(25.95549545, -5.02738237);
	const vec2 kBlueVec2 = vec2(27.81927491, -14.87899417);

//	x = saturate(x); // not available on nvidia
	x = clamp(x, 0.0, 1.0);
	vec4 v4 = vec4( 1.0, x, x * x, x * x * x);
	vec2 v2 = v4.zw * v4.z;
	return vec3(
		dot(v4, kRedVec4)   + dot(v2, kRedVec2),
		dot(v4, kGreenVec4) + dot(v2, kGreenVec2),
		dot(v4, kBlueVec4)  + dot(v2, kBlueVec2)
	);
}

void main() {
	const vec2 screenCoords = (gl_FragCoord.xy/screenRes) * 2.0 - 1.0;
	Ray r = computeCameraRay(screenCoords);

#if 0 // DEBUG for seeing the low-res first-depth pre-pass
	if(useMinDepthTex) {
		//color = texture(minDepthTex, (gl_FragCoord.xy/screenRes) ).xxx/10000.0f;
		color = vec3(getMinT(8) / length(sceneBBoxMax-sceneBBoxMin));
		return;
	}
#endif
	const float epsilon = 1E-3f;
	vec2 t_min_max = vec2(
		useMinDepthTex ? clamp(getMinT(8) - 1E-4f, 0., 1e30f) : 0.,
		1e30); // subtract epsilon for getting proper normal
	vec3 hitNorm;
	traversal_status ts;
	vec4 result = trace_ray(r, t_min_max, projectionFactor, hitNorm, ts);

	// Hit position = camera origin + depth * ray direction
	vec3 hitPos = r.o + result.x * r.d;
	const float cellSize = 2. * rootHalfSide / pow(2, result.y);
	uint nodeIndex = uint(result.w);

//	const float localEpsilon = 1E-3f * result.y;
	
	if (result.x >= 0) // Intersection!!!
	{
		float t = 0;
		if(viewerRenderMode == 0) // ITERATIONS
			t =  1. - (result.z / float(maxIters));
		else if(viewerRenderMode == 1) { // DEPTH
			t = result.x / length(sceneBBoxMax-sceneBBoxMin);
			// gamma correction
			t = 1. - pow(t, 1. / 2.2);
		} else if (viewerRenderMode == 2) // VOXEL LEVELS
			t = log2(2. * rootHalfSide / result.y) / float(LEVELS);
		else if (viewerRenderMode == 3) { // PRETTY

			// ====Base color, based on depth====
			t = 1.0; // - result.x / length(sceneBBoxMax-sceneBBoxMin);


			// ====Voxel normal direction====
//			vec3 localHitPos = hitPos - sceneCenter; // local position, align to grid (through bbox center)
			// Problem: adding r.d * eps causes the ray to cross the boundary to another voxel on around the edges
//			vec3 voxCenter = localHitPos - mod(localHitPos + r.d * localEpsilon, cellSize) + cellSize / 2.0;
//			voxCenter += sceneCenter; // Local -> global position
//			vec3 hitNorm = normalize(voxCenter - hitPos + r.d * localEpsilon);
			
			// The normal vector now points from the voxel center to the hit position like it's a sphere
//			hitNorm = round(hitNorm); // rounding the normal gives a nice "tile" look

			// This converts the spherical normal to a cube normal, setting only its maximum value to 1 (or min to -1)
//			vec3 absHN = abs(hitNorm);
//			float maxAbsHN = max(max(absHN.x, absHN.y), absHN.z);
//			hitNorm = -sign(hitNorm) * vec3(greaterThanEqual(absHN, vec3(maxAbsHN)));

			// This normal computation has some artifacts due to floating point precision
			// Todo: proper normal derivation can be done by looking at the direction of the previous voxel during traversal
			// Yes this is much easier, the data is already there, got it from dda_next

			vec2 rnd3 = normalize(vec2(
				(nodeIndex % 100) / 100.f,
				((3 * nodeIndex) % 200) / 200.f
			));

			const float normalRayLength = 64 * cellSize;
			vec3 smoothNorm = hitNorm;
			float numNormSamplesHit = 1;
			float normSampleRot = 3.14159f * 2.f / float(normSamples);
			float initAngleOffset = mod(result.x * screenCoords.x * nodeIndex * 0.01 + nodeIndex * 0.001 * screenCoords.y, 1);
			traversal_status ts_ignore;
			// Sample multiple times per pixel, in the cone that surrounds the pixel projected into the scene
			for (int i = 0; i < normSamples; i++) {
				vec3 normSample;
				float rad = 0.1 + mod(i / float(normSamples) * nodeIndex * 0.001, 0.9);
				// Sample in a circle around the main ray, at a random radius
				float angle = i * normSampleRot + initAngleOffset;
				vec2 rayPixOffset = rad * vec2(cos(angle), sin(angle)) / screenRes;
				r = computeCameraRay(screenCoords + rayPixOffset);
				vec2 t_min_max = vec2(result.x - normalRayLength * 0.5f, result.x + normalRayLength * 0.5f); // vec2(0, 1e30);
				stack_size = 0;
				vec4 normSampleRes = trace_ray(r, t_min_max, projectionFactor, normSample, ts_ignore);
				if (normSampleRes.x >= 0 && distance(normSampleRes.x, result.x) <= normalRayLength) {
					numNormSamplesHit++;
					smoothNorm += normSample;
				}
			}
			// smoothNorm /= numNormSamplesHit;
			hitNorm = normalize(smoothNorm);

			// ====Diffuse shading====
			vec3 lightDir = normalize(lightPos - hitPos);
			t *= 0.5 + 0.5 * max(dot(hitNorm, lightDir), 0);

			// Light fall-off
			// t -= 0.5 * distance(hitPos, lightPos) / length(sceneBBoxMax-sceneBBoxMin);
			// Traversal depth (looks nice without beam opt)
			// t *=  1. - (result.z / float(maxIters));
		}
	
		color = vec3(t);
//		color = TurboColormap(t);
//		color = hitNorm * 0.5 + vec3(0.5);

		if (nodeIndex == selectedVoxelIndex) {
			// Highlight selected voxel index with blue
			color.r *= 0.5f;
			color.g *= 0.5f;
			color.b = 1.f;
		}
#if (SSVDAG || USSVDAG)
		else if (freqColors) {
			// Visualize ref count by dividing (sorted) index by num of nodes in level
			// Todo: Should make this a sperate render mode
			int hitLev = int(result.y);
			float levSize = 0;
			if (hitLev >  LEVELS - 3) levSize = textureSize(leafNodes); // leaf nodes
			else if (hitLev == LEVELS - 3) levSize = textureSize(innerNodes) * 4 - levelOffsets[hitLev]; // lev above leaves
			else 					  levSize = levelOffsets[hitLev + 1] - levelOffsets[hitLev];
			float indexInLevel = nodeIndex / levSize;
			color *= TurboColormap(1.0f - indexInLevel - 0.05f);
		}
#endif
		// Assign random colors based on the index of a node
		else if (randomColors) {
			vec3 c = ts.attr_sum * 255.; //  get_header_attrs(ts.hdr) * 255.; //
			int y = int(c.x), cr = int(c.y), cb = int(c.z);
			// yuv to rgb
			vec3 randomColor = vec3(
				CYCbCr2R(y, cr, cb),
				CYCbCr2G(y, cr, cb),
				CYCbCr2B(y, cr, cb)
			) / vec3(255.);

			//int hitLev = int(result.y);
			//vec3 randomColor = vec3(stack_size / 15.);
			//vec3 randomColor = get_header_attrs(stack[stack_size - 1].y);
			// Todo: we need to sample neighboring nodes as well to avoid blocky colors... hmm
			// vec3 randomColor = ts.attr_sum / float(stack_size); //result.y; // sum of attrs or all nodes / number of levels

			// old idea: store attr in pointers
			// new: store in header padding

//			randomColor = vec3(clamp(ts.hdr / 255.0, 0, 1));

//			if (drawLevel == LEVELS - 1) { // for attribute DAG: attributes are at leaf level
//				// This node represents 8 voxels, each one has different attributes
//				// How do I find out which voxel of this node is intersected? Lets just pick the first one for now
//					
//				ts.child_linear_index = voxel_to_linear_idx(ts.mirror_mask, ts.local_idx, ts.current_node_size);
//				fetch_child_index_in(ts);
//				fetch_data(ts);
//
//				randomColor = vec3(ts.hdr / 255.);
//			} else {
//				randomColor = normalize(vec3(
//					(nodeIndex % 100) / 100.f,
//					((3 * nodeIndex) % 200) / 200.f,
//					((2 * nodeIndex) % 300) / 300.f
//				));
//			}

			color *= randomColor;
		} 

		// ====Shadow====
//		if (result.x >= 0 && enableShadows) {
//			// Trace ray to light pos
//			const vec3 p = hitPos + hitNorm * cellSize * 0.5f; // center of voxel + half voxel size in dir of normal
//			const vec3 lightToP = p - lightPos;
//			const float lightToPLength = length(lightToP);
//
//			const float sProjFactor = lightToPLength / cellSize;
//
//			r.o = p;
//			r.d = -normalize(lightToP);
//			vec2 t_min_max = vec2(0, lightToPLength);
//			stack_size = 0; // not sure why it needs to be reset, but else it causes artifacts
//			traversal_status ts_ignore;
//			vec4 shd_result = trace_ray(r, t_min_max * 0.5, sProjFactor, hitNorm, ts_ignore);
//
//			// Add shadows
//			color *= (shd_result.x > 0) ? 0.8 : 1;
////				t = shd_result.x > 0 ? t : (1.0 - shd_result.x / length(sceneBBoxMax-sceneBBoxMin));
//		}
	}
	else {
		if (result.x == -1. || result.x == -2.
		// ) // inside BBox, but no intersection
		// { 
		// 	float t = (viewerRenderMode == 0) ? 1. - (result.z / float(maxIters)) : 1.;
		// 	color = vec3(t,t*0.5,1);
		// }
		// else if (
			|| result.x == -4.) // no intersection, out of BBox => background
		{ 
			color = vec3(1); // 0.5 * (screenCoords.y + 1) * vec3(0.3,0.2,0.9);
		}
		
		else if (result.x == -3.) // too many iterations
		{
			color = vec3(1.,0,0);
		}
		else // never should happen
		{
			discard; 
		}
	}
}

#elif DEPTH_MODE

void main() {
	const vec2 screenCoords = (gl_FragCoord.xy/screenRes) * 2.0 - 1.0;
	Ray r = computeCameraRay(screenCoords);	
	vec2 t_min_max = vec2(useMinDepthTex ? getMinT(8) : 0, 1e30);
	vec3 norm;
	traversal_status ts_ignore;
	vec4 result = trace_ray(r, t_min_max, projectionFactor, norm, ts_ignore);
	if (result.x > 0) // Intersection!!!
		output_t = result.xyz;
	else
		discard;
}

#elif SHADOW_MODE 

void main() {
	const ivec2 coord = ivec2(gl_FragCoord.xy);
	const vec3 hitPos = texelFetch(hitPosTex, coord, 0).xyz;
	if(hitPos == vec3(0,0,0)) discard;
	const float cellSize = texelFetch(depthTex, coord, 0).y;
	const vec3 hitNorm = texelFetch(hitNormTex, coord, 0).xyz;

	const vec3 p = hitPos + hitNorm * cellSize * 1.0;

	Ray r;
	r.o = lightPos;
	const vec3 lightToP = p - lightPos;
	const float lightToPLength = length(lightToP);
	r.d = normalize(lightToP);
	vec2 t_min_max = vec2(0, lightToPLength);
	const float projFactor = lightToPLength / cellSize;
	vec3 norm;
	traversal_status ts_ignore;
	vec4 result = trace_ray(r, t_min_max, projFactor, norm, ts_ignore);

	output_t = (result.x > 0) ? 0.0 : 1.0;
}

#elif AO_MODE 

void main() {
	const ivec2 coord = ivec2(gl_FragCoord.xy);
	vec3 hitPos = texelFetch(hitPosTex, coord, 0).xyz;
	if(hitPos == vec3(0,0,0)) discard;
	const vec3 hitNorm = texelFetch(hitNormTex, coord, 0).xyz;
	const float cellSize = texelFetch(depthTex, coord, 0).y;
	hitPos += hitNorm * 1e-3;

	Ray aoRay;
	aoRay.o = hitPos;

	//const vec3 rvec = gl_FragCoord.xyz; // TODO review this !!!
	const vec3 rvec = vec3(0,1,0); // TODO review this !!!
	//const vec3 rvec = normalize(vec3(gl_FragCoord.xy, gl_FragCoord.x * gl_FragCoord.y));
	//const vec3 rvec = normalize(vec3(gl_FragCoord.x, 0, gl_FragCoord.y));
	const vec3 tangent = normalize(rvec - hitNorm * dot(rvec, hitNorm));
	const vec3 bitangent = cross(hitNorm, tangent);
	//const mat3 tnb = mat3(tangent, bitangent, hitNorm);
	const mat3 tnb = mat3(tangent, hitNorm, bitangent);
	float largo = 0.5;
	const vec2 t_min_max = vec2(0,largo);
	const float projFactor = largo / cellSize;
	float k = 0;
	
	traversal_status ts_ignore;

	for(int i=0; i<numAORays; i++) {
		aoRay.d = normalize(tnb *  hsSamples[(i+int(gl_FragCoord.x*gl_FragCoord.y))%N_HS_SAMPLES]);
		//aoRay.d = normalize(tnb *  hsSamples[i]);
		vec3 norm;
		vec4 result = trace_ray(aoRay, t_min_max, projFactor, norm, ts_ignore);
		if(result.x > 0) k += 1.0;// - (result.x/0.3);
	}
	float visibility = (numAORays>0) ? 1.0 - (k/float(numAORays)) : 1.0;
	output_t = visibility * visibility;
}

#endif
