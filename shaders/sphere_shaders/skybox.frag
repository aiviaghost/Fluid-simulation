#version 410

uniform sampler2D diffuse_texture;
uniform int has_diffuse_texture;

in VS_OUT {
	vec3 vertex;
	vec3 normal;
} fs_in;

out vec4 frag_color;

uniform samplerCube cubemap;

void main()
{
	frag_color = texture(cubemap, fs_in.vertex);
}
