#include "parametric_shapes.hpp"
#include "core/Log.h"

#include <glm/glm.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

bonobo::mesh_data
parametric_shapes::createSphere(float const radius,
                                unsigned int const longitude_split_count,
                                unsigned int const latitude_split_count)
{

	//! \todo Implement this function

	auto const longitude_slice_edges_count = longitude_split_count + 1u;
	auto const latitude_slice_edges_count = latitude_split_count + 1u;
	auto const latitude_slice_vertices_count = latitude_slice_edges_count + 1u;
	auto const longitude_slice_vertices_count = longitude_slice_edges_count + 1u;
	auto const vertices_nb = latitude_slice_vertices_count * longitude_slice_vertices_count;

	auto vertices = std::vector<glm::vec3>(vertices_nb);
	auto normals = std::vector<glm::vec3>(vertices_nb);
	auto texcoords = std::vector<glm::vec3>(vertices_nb);
	auto tangents = std::vector<glm::vec3>(vertices_nb);
	auto binormals = std::vector<glm::vec3>(vertices_nb);

	float const d_theta = glm::two_pi<float>() / (static_cast<float>(longitude_slice_edges_count));
	float const d_phi = glm::pi<float>() / (static_cast<float>(latitude_slice_edges_count));

	size_t index = 0u;
	float theta = 0.0f;
	float phi = 0.0f;
	float const r = radius;
	for (unsigned int i = 0u; i < longitude_slice_vertices_count; ++i) {
		float const cos_theta = std::cos(theta);
		float const sin_theta = std::sin(theta);

		phi = 0.0f;
		for (unsigned int j = 0u; j < latitude_slice_vertices_count; ++j) {
			float const cos_phi = std::cos(phi);
			float const sin_phi = std::sin(phi);
			
			// vertex
			vertices[index] = glm::vec3(r * sin_theta * sin_phi,
				-r * cos_phi,
				r * cos_theta * sin_phi);

			// texture coordinates
			texcoords[index] = glm::vec3(static_cast<float>(i) / (static_cast<float>(longitude_slice_vertices_count)),
				static_cast<float>(j) / (static_cast<float>(latitude_slice_vertices_count)),
				0.0f);

			// tangent
			auto const t = glm::normalize(glm::vec3(r * cos_theta, 0, -r * sin_theta));

			tangents[index] = t;

			// binormal
			auto const b = glm::normalize(glm::vec3(r * sin_theta * cos_phi, r * sin_phi, r * cos_theta * cos_phi));
			binormals[index] = b;

			// normal
			auto const n = glm::cross(t, b);

			normals[index] = n;

			phi += d_phi;
			++index;
		}

		theta += d_theta;
	}

	// create index array
	auto index_sets = std::vector<glm::uvec3>(2u * longitude_slice_edges_count * latitude_slice_edges_count);

	// generate indices iteratively
	index = 0u;
	for (unsigned int i = 0u; i < longitude_slice_edges_count; ++i)
	{
		for (unsigned int j = 0u; j < latitude_slice_edges_count; ++j)
		{
			index_sets[index] = glm::uvec3(latitude_slice_vertices_count * (i + 0u) + (j + 0u),
										   latitude_slice_vertices_count * (i + 1u) + (j + 1u),
										   latitude_slice_vertices_count * (i + 0u) + (j + 1u));
			++index;

			index_sets[index] = glm::uvec3(latitude_slice_vertices_count * (i + 0u) + (j + 0u),
										   latitude_slice_vertices_count * (i + 1u) + (j + 0u),
										   latitude_slice_vertices_count * (i + 1u) + (j + 1u));
			++index;
		}
	}

	bonobo::mesh_data data;
	glGenVertexArrays(1, &data.vao);
	assert(data.vao != 0u);
	glBindVertexArray(data.vao);

	auto const vertices_offset = 0u;
	auto const vertices_size = static_cast<GLsizeiptr>(vertices.size() * sizeof(glm::vec3));
	auto const normals_offset = vertices_size;
	auto const normals_size = static_cast<GLsizeiptr>(normals.size() * sizeof(glm::vec3));
	auto const texcoords_offset = normals_offset + normals_size;
	auto const texcoords_size = static_cast<GLsizeiptr>(texcoords.size() * sizeof(glm::vec3));
	auto const tangents_offset = texcoords_offset + texcoords_size;
	auto const tangents_size = static_cast<GLsizeiptr>(tangents.size() * sizeof(glm::vec3));
	auto const binormals_offset = tangents_offset + tangents_size;
	auto const binormals_size = static_cast<GLsizeiptr>(binormals.size() * sizeof(glm::vec3));
	auto const bo_size = static_cast<GLsizeiptr>(vertices_size
		+ normals_size
		+ texcoords_size
		+ tangents_size
		+ binormals_size
		);
	glGenBuffers(1, &data.bo);
	assert(data.bo != 0u);
	glBindBuffer(GL_ARRAY_BUFFER, data.bo);
	glBufferData(GL_ARRAY_BUFFER, bo_size, nullptr, GL_STATIC_DRAW);

	glBufferSubData(GL_ARRAY_BUFFER, vertices_offset, vertices_size, static_cast<GLvoid const*>(vertices.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::vertices));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::vertices), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(0x0));

	glBufferSubData(GL_ARRAY_BUFFER, normals_offset, normals_size, static_cast<GLvoid const*>(normals.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::normals));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::normals), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(normals_offset));

	glBufferSubData(GL_ARRAY_BUFFER, texcoords_offset, texcoords_size, static_cast<GLvoid const*>(texcoords.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::texcoords));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::texcoords), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(texcoords_offset));

	glBufferSubData(GL_ARRAY_BUFFER, tangents_offset, tangents_size, static_cast<GLvoid const*>(tangents.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::tangents));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::tangents), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(tangents_offset));

	glBufferSubData(GL_ARRAY_BUFFER, binormals_offset, binormals_size, static_cast<GLvoid const*>(binormals.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::binormals));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::binormals), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(binormals_offset));

	glBindBuffer(GL_ARRAY_BUFFER, 0u);

	data.indices_nb = static_cast<GLsizei>(index_sets.size() * 3u);
	glGenBuffers(1, &data.ibo);
	assert(data.ibo != 0u);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data.ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(index_sets.size() * sizeof(glm::uvec3)), reinterpret_cast<GLvoid const*>(index_sets.data()), GL_STATIC_DRAW);

	glBindVertexArray(0u);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0u);

	return data;
}


bonobo::mesh_data parametric_shapes::createHighTesselationQuad(
	float const width, float const height,
	unsigned int const horizontal_split_count,
	unsigned int const vertical_split_count) {

	int num_tris = horizontal_split_count * 2 * vertical_split_count;
	auto index_sets = std::vector<glm::uvec3>(0);

	int n = horizontal_split_count + 1;
	int m = vertical_split_count + 1;
	//std::swap(n, m);
	int vertices_nb = n * m;

	auto vertices = std::vector<glm::vec3>(vertices_nb);
	auto texcoords = std::vector<glm::vec3>(vertices_nb);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			float s = i * height;
			float t = j * width;

			int ind = j + i * m;
			vertices[ind] = glm::vec3(s, t, 0.0);
			texcoords[ind] = glm::vec3((1.0 * i) / (1.0 * vertical_split_count), (1.0 * j) / horizontal_split_count, 0.0);


			if (i < n - 1 && j < m - 1) {
				index_sets.push_back(glm::vec3(ind, ind + 1, ind + 1 + m));
				index_sets.push_back(glm::vec3(ind, ind + 1 + m, ind + m));
			}
		}
	}



	bonobo::mesh_data data;
	glGenVertexArrays(1, &data.vao);
	assert(data.vao != 0u);
	glBindVertexArray(data.vao);

	auto const vertices_offset = 0u;
	auto const vertices_size = static_cast<GLsizeiptr>(vertices.size() * sizeof(glm::vec3));

	auto const texcoords_offset = vertices_size;
	auto const texcoords_size = static_cast<GLsizeiptr>(texcoords.size() * sizeof(glm::vec3));

	auto const bo_size = static_cast<GLsizeiptr>(vertices_size + texcoords_size);


	glGenBuffers(1, &data.bo);
	assert(data.bo != 0u);
	glBindBuffer(GL_ARRAY_BUFFER, data.bo);
	glBufferData(GL_ARRAY_BUFFER, bo_size, nullptr, GL_STATIC_DRAW);

	glBufferSubData(GL_ARRAY_BUFFER, vertices_offset, vertices_size, static_cast<GLvoid const*>(vertices.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::vertices));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::vertices), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(0x0));

	glBufferSubData(GL_ARRAY_BUFFER, texcoords_offset, texcoords_size, static_cast<GLvoid const*>(texcoords.data()));
	glEnableVertexAttribArray(static_cast<unsigned int>(bonobo::shader_bindings::texcoords));
	glVertexAttribPointer(static_cast<unsigned int>(bonobo::shader_bindings::texcoords), 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<GLvoid const*>(texcoords_offset));

	glBindBuffer(GL_ARRAY_BUFFER, 0u);

	data.indices_nb = static_cast<GLsizei>(index_sets.size() * 3u);
	glGenBuffers(1, &data.ibo);
	assert(data.ibo != 0u);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data.ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(index_sets.size() * sizeof(glm::vec3)), reinterpret_cast<GLvoid const*>(index_sets.data()), GL_STATIC_DRAW);

	glBindVertexArray(0u);
	//glBindBuffer(GL_ARRAY_BUFFER, 0u);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0u);

	return data;
}
