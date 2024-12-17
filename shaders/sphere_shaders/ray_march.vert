#version 430

// Remember how we enabled vertex attributes in assignment 2 and attached some
// data to each of them, here we retrieve that data. Attribute 0 pointed to the
// vertices inside the OpenGL buffer object, so if we say that our input
// variable `vertex` is at location 0, which corresponds to attribute 0 of our
// vertex array, vertex will be effectively filled with vertices from our
// buffer.
// Similarly, normal is set to location 1, which corresponds to attribute 1 of
// the vertex array, and therefore will be filled with normals taken out of our
// buffer.
layout (location = 0) in vec3 vertex;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 texcoords;
layout (location = 3) in vec3 tangent;
layout (location = 4) in vec3 binormal;

struct Particle {
	vec3 position;
	vec3 predicted_position;
	vec3 velocity;
};

struct Colour {
	vec3 colour;
};

layout(std430, binding = 0) buffer ParticleBuffer {
    Particle particles[];
};

layout(std430, binding = 6) buffer ColourBuffer {
    Colour colours[];
};

uniform mat4 vertex_model_to_world;
uniform mat4 normal_model_to_world;
uniform mat4 vertex_world_to_clip;

uniform int has_textures;
uniform int has_diffuse_texture;
uniform int has_opacity_texture;

// This is the custom output of this shader. If you want to retrieve this data
// from another shader further down the pipeline, you need to declare the exact
// same structure as in (for input), with matching name for the structure
// members and matching structure type. Have a look at
// shaders/EDAF80/diffuse.frag.
out VS_OUT {
	vec3 vertex;
	vec3 normal;
	vec3 colour;
	vec3 texcoords;
} vs_out;


void main()
{
	vec3 shift = vec3(-1.0, -1.0, -1.0);
	vs_out.texcoords = texcoords;
	gl_Position = vec4(vertex + shift, 1.0);
}
