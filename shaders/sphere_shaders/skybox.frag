#version 410

uniform sampler2D diffuse_texture;
uniform int has_diffuse_texture;
uniform vec3 camera_position;

in VS_OUT {
	vec3 vertex;
	vec3 normal;
} fs_in;

out vec4 frag_color;

uniform samplerCube cubemap;

void main()
{
	frag_color = texture(cubemap, normalize(fs_in.vertex - camera_position));
}
