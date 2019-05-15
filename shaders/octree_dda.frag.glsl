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

uniform int attrBit = 0; // Temporary uniform for indicating which attribute bit to show

#if VIEWER_MODE
uniform int viewerRenderMode;
uniform uint selectedVoxelIndex;
uniform bool randomColors;
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

struct Ray {
	vec3 o;
	vec3 d;
};

struct traversal_status {
	float t_current;
	int node_index;
        uint hdr; // node header (contains child bitmask)
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
};

//////////////////////////// STACK STUFF ////////////////////

#define MAX_STACK_DEPTH (INNER_LEVELS+1)

ivec3 stack[MAX_STACK_DEPTH];
uint stack_size = 0;

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
	const uvec4 tmp = texelFetch(nodes, idx/4);
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
	const uvec4 tmp = texelFetch(nodes, idx/4);
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
	const int realIdx = idx/4;
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
	const int realIdx = idx/4;
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

void dda_next(inout traversal_status ts) {
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

// In the attribute DAG, we need to go down the attribute root to get each bit of the attributes
void go_down_attr_root(in const Ray r, inout traversal_status ts, const int attr_bit) {
    const int childPtrPos = bitCount((ts.hdr & 0xFF) >> attr_bit);
    ts.node_index = int(myFetch(ts.node_index + childPtrPos));
}

void go_down_one_level(in const Ray r, inout traversal_status ts) {
	++ts.level;
	ts.cell_size *= 0.5;
	
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
	ts.child_linear_index =  voxel_to_linear_idx(ts.mirror_mask, ts.local_idx, ts.current_node_size);
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
//      W: node index (-1 => no intersection)

vec4 trace_ray(in Ray r, in vec2 t_min_max, const in float projection_factor, const in int attr_bit) {
	
	if (!transform_ray(r, t_min_max))
		return vec4(-4.0,0,0,-1); // out of scene Bbox
	
	const float scale = 2.0 * rootHalfSide;
	traversal_status ts;
	ts.t_current = t_min_max.x;
	init(r, ts);

	// For attr dag: start at level 1, without modifying the ray/traversal status
	go_down_attr_root(r, ts, attr_bit);

	int iteration_count = 0;
	const uint max_level = min(INNER_LEVELS - 1, drawLevel-1); // remove -1 at INNER_LEVELS when not descending attr root
	do {
		bool full_voxel = fetch_voxel_bit(ts);
	  
		if (!full_voxel) {
			dda_next(ts);
			if (!in_bounds(ts.local_idx, ts.current_node_size)) {
				if (stack_is_empty()) {
					return vec4(-1.,0,float(iteration_count),-1); // inside scene BBox, but no intersection
				}
				up_in(r, ts);
			}
		} else {
			const bool hit = (ts.level >= max_level || resolution_ok(ts.t_current, ts.cell_size, projection_factor)); 
			if (hit) {
				return vec4(ts.t_current * scale, ts.cell_size * scale, float(iteration_count), ts.node_index);  // intersection
			} else {
				down_in(r, ts);
				fetch_data(ts);
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
	
	const ivec2 p = ivec2(gl_FragCoord.xy / delta);

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

void main() {
	const vec2 screenCoords = (gl_FragCoord.xy/screenRes) * 2.0 - 1.0;
	Ray r = computeCameraRay(screenCoords);


	// Todo: Ray cast 8 times (for each bit)
	// Start at the i-th child instead of at the root

#if 0 // DEBUG for seeing the low-res first-depth pre-pass
	if(useMinDepthTex) {
		//color = texture(minDepthTex, (gl_FragCoord.xy/screenRes) ).xxx/10000.0f;
		color = vec3(getMinT(8) / length(sceneBBoxMax-sceneBBoxMin));
		return;
	}
#endif
	vec2 t_min_max = vec2(useMinDepthTex ? getMinT(8) : 0, 1e30);
	vec4 result = trace_ray(r, t_min_max, projectionFactor, randomColors ? 0 : attrBit);

	int attr = 0;
	if (result.x >= 0) {
		attr |= 1 << 0;
	}
	if (randomColors) {
		// Todo: Only count as intersection when they are intersecting the same position (the nearest one)
		// probably requires one extra attr bit, that is the original DAG containing all geometry
		for (int i = 1; i < 8; ++i) {
			// Todo: reset stack?
//			stack_size = 0;
//			r = computeCameraRay(screenCoords);
//			t_min_max = vec2(useMinDepthTex ? getMinT(8) : 0, 1e30);
			result = trace_ray(r, t_min_max, projectionFactor, i);
			if (result.x >= 0) { // Intersection!!!
				attr |= 1 << i;
			}
		}
		float t = (attr / 255.0);
		color = vec3(t, t, t);
		return;
	}


	if (attr > 0)
	{
		float t = 0;
		if(viewerRenderMode == 0) // ITERATIONS
			t =  1. - (result.z / float(maxIters));
		else if(viewerRenderMode == 1) // DEPTH
			t = result.x / length(sceneBBoxMax-sceneBBoxMin);
		else if (viewerRenderMode == 2) // VOXEL LEVELS
			t = log2(2. * rootHalfSide / result.y) / float(LEVELS);

		color = vec3(t, t, t);


		/////////////////////////////////
		// CUSTOMIZED:
		/////////////////////////////////
		// Assign random colors based on the index of a node
//		if (randomColors && selectedVoxelIndex == 0 && result.w > 0) {
////			color *= randomColor(uint(result.w));
//			uint nodeIndex = uint(result.w);
//			vec3 randomColor = normalize(vec3(
//				(nodeIndex % 100) / 100.f,
//				((3 * nodeIndex) % 200) / 200.f,
//				((2 * nodeIndex) % 300) / 300.f
//			));
//			color *= randomColor;
//		} else if (result.w == selectedVoxelIndex) {
//			// Highlight selected voxel index with blue
//			color.r *= 0.5f;
//                        color.g *= 0.5f;
//			color.b = 1.f;
//		}
	}
	else {
		if (result.x == -4.) // no intersection, out of BBox => background
		{ 
			color = 0.5 * (screenCoords.y + 1) * vec3(0.3,0.2,0.9);
		}
		else if (result.x == -1. || result.x == -2.) // inside BBox, but no intersection
		{ 
			float t = (viewerRenderMode == 0) ? 1. - (result.z / float(maxIters)) : 1.;
			color = vec3(t,t*0.5,0);
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
	vec4 result = trace_ray(r, t_min_max, projectionFactor, attrBit);
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
	vec4 result = trace_ray(r, t_min_max, projFactor, attrBit);

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

	for(int i=0; i<numAORays; i++) {
		aoRay.d = normalize(tnb *  hsSamples[(i+int(gl_FragCoord.x*gl_FragCoord.y))%N_HS_SAMPLES]);
		//aoRay.d = normalize(tnb *  hsSamples[i]);
		vec4 result = trace_ray(aoRay, t_min_max, projFactor, attrBit);
		if(result.x > 0) k += 1.0;// - (result.x/0.3);
	}
	float visibility = (numAORays>0) ? 1.0 - (k/float(numAORays)) : 1.0;
	output_t = visibility * visibility;
}

#endif
