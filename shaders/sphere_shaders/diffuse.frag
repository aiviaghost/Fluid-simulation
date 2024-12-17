#version 410

uniform vec3 light_position;
uniform vec3 diffuse_colour;

in VS_OUT {
	vec3 vertex;
	vec3 normal;
	vec3 colour;
} fs_in;

out vec4 frag_color;

void main()
{
	vec3 L = normalize(light_position - fs_in.vertex);
	frag_color = vec4(fs_in.colour, 1.0) * clamp(dot(normalize(fs_in.normal), L), 0.0, 1.0);
	frag_color = vec4(fs_in.colour, 1.0);
}
