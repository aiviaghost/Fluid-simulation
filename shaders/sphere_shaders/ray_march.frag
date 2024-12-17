#version 430

// ----------------------------------------------------------------------------
//
// uniforms
//
// ----------------------------------------------------------------------------

struct Particle {
	vec3 position;
	vec3 predicted_position;
	vec3 velocity;
};

layout(std430, binding = 0) buffer ParticleBuffer {
    Particle particles[];
};

layout(std430, binding = 3) buffer SpatialBuffer {
    ivec2 spatial[];
};

layout(std430, binding = 4) buffer StartBuffer {
    int start_inds[];
};

uniform int num_particles;
uniform float smoothing_radius;
int PRIME1 = 86183;
int PRIME2 = 7475723;
int PRIME3 = 447409;
float PI = 3.14159265358979;
float mass = 1.0;

uniform float half_width;
uniform float half_height;
uniform float half_depth;
uniform float step_size;
uniform float density_multiplier;

uniform vec3 light_position;
uniform vec3 diffuse_colour;
uniform mat4 vertex_clip_to_world;
uniform vec3 camera_position;

in VS_OUT {
	vec3 vertex;
	vec3 normal;
	vec3 colour;
	vec3 texcoords;
} fs_in;


out vec4 frag_color;

// ----------------------------------------------------------------------------
//
// functions
//
// ----------------------------------------------------------------------------


int good_mod(int a, int m){
	if (a >= 0) return a % m;
	return (-(-a % m) + m) % m;
}

int point_to_hash(vec3 point){
	int x = int(point.x / smoothing_radius);
	int y = int(point.y / smoothing_radius);
	int z = int(point.z / smoothing_radius);
	int hash = x * PRIME1 + y * PRIME2 + z * PRIME3;
	return good_mod(hash, num_particles);
}

float smoothing_kernel(float dist, float r) {
	if (dist >= r) return 0.0;
	float volume = (PI * pow(r, 4.0)) / 6.0;
	float x = r - dist;
	return x * x / volume;
}

float near_density_kernel(float dist, float r) {
	if (dist >= r) return 0.0;
	float v = r - dist;
	float volume = PI * pow(r, 5.0) / 10.0;
	return v * v * v / volume;
};

float calculate_density(vec3 point) {
	float density = 0.0;

	int hash = point_to_hash(point);

	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for(int k = -1; k <= 1; k++) {
				int cur_hash = hash + i * PRIME1 + j * PRIME2 + k * PRIME3;
				cur_hash = good_mod(cur_hash, num_particles);
				int start_spatial_ind = start_inds[cur_hash];
				if (start_spatial_ind == -1) continue;
				for (int spatial_ind = start_spatial_ind; spatial_ind < num_particles && spatial[spatial_ind].x == cur_hash; spatial_ind++) {
					float dist = length(particles[spatial[spatial_ind].y].predicted_position - point);

					density += mass * smoothing_kernel(dist, smoothing_radius);
				}
			}
		}
	}

	return density;
}


float eps = 0.000001;
vec2 intersections(vec3 orig, vec3 dir){
	float mini = 10000000.0;
	float maxi = -1.0;

	vec3 bounds = vec3(half_width, half_height, half_depth);

	float l;
	vec3 p;

	if(abs(dir.x) > eps){
		l = (bounds.x - orig.x) / dir.x;
		p = orig + l * dir;
		if(abs(p.y) < bounds.y && abs(p.z) < bounds.z){		
			mini = min(mini, l);
			maxi = max(maxi, l);
		}
		l = (-bounds.x - orig.x) / dir.x;
		p = orig + l * dir;
		if(abs(p.y) < bounds.y && abs(p.z) < bounds.z){		
			mini = min(mini, l);
			maxi = max(maxi, l);
		}
	}

	if(abs(dir.y) > eps){
		l = (bounds.y - orig.y) / dir.y;
		p = orig + l * dir;
		if(abs(p.x) < bounds.x && abs(p.z) < bounds.z){		
			mini = min(mini, l);
			maxi = max(maxi, l);
		}
		l = (-bounds.y - orig.y) / dir.y;
		p = orig + l * dir;
		if(abs(p.x) < bounds.x && abs(p.z) < bounds.z){		
			mini = min(mini, l);
			maxi = max(maxi, l);
		}
	}

	if(abs(dir.z) > eps){
		l = (bounds.z - orig.z) / dir.z;
		p = orig + l * dir;
		if(abs(p.y) < bounds.y && abs(p.x) < bounds.x){		
			mini = min(mini, l);
			maxi = max(maxi, l);
		}
		l = (-bounds.z - orig.z) / dir.z;
		p = orig + l * dir;
		if(abs(p.y) < bounds.y && abs(p.x) < bounds.x){		
			mini = min(mini, l);
			maxi = max(maxi, l);
		}
	}

	return vec2(mini, maxi);
}



void main()
{


	vec2 uv = fs_in.texcoords.xy * 2.0 - 1.0;
	vec4 clip_space = vec4(uv, -1.0, 1.0);
	vec4 world_space = vertex_clip_to_world * clip_space;
	world_space /= world_space.w;

	vec3 orig = camera_position;
	vec3 dir = normalize(world_space.xyz - orig);

	vec2 inter = intersections(orig, dir);


	frag_color = vec4(fs_in.texcoords.xy, 0.0, 1.0);
	frag_color = vec4(dir, 1.0);

	if(inter.y > 0){
		frag_color = vec4(1.0);
	}

	float tot_density = 0.0;

	if(inter.x < 0) inter.x = 0;
	for(float lambda = inter.x + eps; lambda < inter.y - eps; lambda += step_size){
		vec3 point = orig + lambda * dir;
		float density = calculate_density(point) * density_multiplier * step_size;
		tot_density += density;
	}

	vec3 water = vec3(94.0, 219.0, 228.0) / 255.0;
	frag_color = vec4(tot_density * water, 1.0);

}
