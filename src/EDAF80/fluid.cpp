#include "fluid.hpp"
#include "parametric_shapes.hpp"

#include "config.hpp"
#include "core/Bonobo.h"
#include "core/FPSCamera.h"
#include "core/helpers.hpp"
#include "core/node.hpp"
#include "core/ShaderProgramManager.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <imgui.h>
#include <tinyfiledialogs.h>

#include <clocale>
#include <stdexcept>

#include <windows.h>
#include <ppl.h>

// particle renderer
#include "core/Log.h"
#include "core/opengl.hpp"

#include <glm/gtc/matrix_transform.hpp>

struct Particle {
	glm::vec3 position;
	float padding0;
	glm::vec3 predicted_position;
	float padding1;
	glm::vec3 velocity;
	float padding2;
	float density;
	float near_density;
	glm::vec2 padding3;
};

class ParticleRenderer {
private:
	// Geometry data
	GLuint _vao{ 0u };
	GLsizei _vertices_nb{ 0u };
	GLsizei _indices_nb{ 0u };
	GLenum _drawing_mode{ GL_TRIANGLES };
	bool _has_indices{ false };

	// Program data
	GLuint const* _program{ nullptr };
	std::function<void(GLuint)> _set_uniforms;

	// Material data
	bonobo::material_data _constants;

	// Transformation data
	TRSTransformf _transform;

	// Debug data
	std::string _name{ "Render un-named node" };

	void render(glm::mat4 const& view_projection, glm::mat4& world, GLuint program, std::function<void(GLuint)> const& set_uniforms, std::vector<glm::vec3> const& positions, std::vector<glm::vec3> const& colours) const {
		if (_vao == 0u || program == 0u)
			return;

		utils::opengl::debug::beginDebugGroup(_name);

		glUseProgram(program);

		//glm::mat4 const& world = glm::mat4(1.0);

		/*for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << world[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::exit(0);*/

		auto const normal_model_to_world = glm::transpose(glm::inverse(world));

		set_uniforms(program);

		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_model_to_world"), 1, GL_FALSE, glm::value_ptr(world));
		glUniformMatrix4fv(glGetUniformLocation(program, "normal_model_to_world"), 1, GL_FALSE, glm::value_ptr(normal_model_to_world));
		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_world_to_clip"), 1, GL_FALSE, glm::value_ptr(view_projection));

		/*for (size_t i = 0u; i < _textures.size(); ++i) {
			auto const& texture = _textures[i];
			glActiveTexture(GL_TEXTURE0 + static_cast<GLenum>(i));
			glBindTexture(std::get<2>(texture), std::get<1>(texture));
			glUniform1i(glGetUniformLocation(program, std::get<0>(texture).c_str()), static_cast<GLint>(i));

			std::string texture_presence_var_name = "has_" + std::get<0>(texture);
			glUniform1i(glGetUniformLocation(program, texture_presence_var_name.c_str()), 1);
		}*/

		glUniform3fv(glGetUniformLocation(program, "diffuse_colour"), 1, glm::value_ptr(_constants.diffuse));
		glUniform3fv(glGetUniformLocation(program, "specular_colour"), 1, glm::value_ptr(_constants.specular));
		glUniform3fv(glGetUniformLocation(program, "ambient_colour"), 1, glm::value_ptr(_constants.ambient));
		glUniform3fv(glGetUniformLocation(program, "emissive_colour"), 1, glm::value_ptr(_constants.emissive));
		glUniform1f(glGetUniformLocation(program, "shininess_value"), _constants.shininess);
		glUniform1f(glGetUniformLocation(program, "index_of_refraction_value"), _constants.indexOfRefraction);
		glUniform1f(glGetUniformLocation(program, "opacity_value"), _constants.opacity);



		glBindVertexArray(_vao);

		unsigned int instanceVBO;
		const unsigned int vertex_offset_index = 5;
		glGenBuffers(1, &instanceVBO);
		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * positions.size(), static_cast<GLvoid const*>(positions.data()), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glEnableVertexAttribArray(vertex_offset_index);
		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glVertexAttribPointer(vertex_offset_index, 3, GL_FLOAT, GL_FALSE, 0, (void*)0/*960*/);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glVertexAttribDivisor(vertex_offset_index, 1);


		const unsigned int colour_index = 6;
		glGenBuffers(1, &instanceVBO);
		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * colours.size(), static_cast<GLvoid const*>(colours.data()), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glEnableVertexAttribArray(colour_index);
		glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
		glVertexAttribPointer(colour_index, 3, GL_FLOAT, GL_FALSE, 0, (void*)0/*960*/);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glVertexAttribDivisor(colour_index, 1);



		if (_has_indices)
			glDrawElementsInstanced(_drawing_mode, _indices_nb, GL_UNSIGNED_INT, reinterpret_cast<GLvoid const*>(0x0), positions.size());
		else
			glDrawArrays(_drawing_mode, 0, _vertices_nb);
		glBindVertexArray(0u);

		/*for (auto const& texture : _textures) {
			glBindTexture(std::get<2>(texture), 0);
			glUniform1i(glGetUniformLocation(program, std::get<0>(texture).c_str()), 0);

			std::string texture_presence_var_name = "has_" + std::get<0>(texture);
			glUniform1i(glGetUniformLocation(program, texture_presence_var_name.c_str()), 0);
		}*/

		glUseProgram(0u);

		utils::opengl::debug::endDebugGroup();
	}

public:
	void render(glm::mat4 const& view_projection, std::vector<glm::vec3> const& positions, std::vector<glm::vec3> const& colours) const {
		if (_program != nullptr) {
			render(view_projection, _transform.GetMatrix(), *_program, _set_uniforms, positions, colours);
		}
	}

	void set_program(GLuint const* const program, std::function<void(GLuint)> const& set_uniforms) {
		if (program == nullptr) {
			LogError("Program can not be a null pointer; this operation will be discarded.");
			return;
		}

		_program = program;
		_set_uniforms = set_uniforms;
	}

	void set_geometry(bonobo::mesh_data const& shape) {
		_vao = shape.vao;
		_vertices_nb = static_cast<GLsizei>(shape.vertices_nb);
		_indices_nb = static_cast<GLsizei>(shape.indices_nb);
		_drawing_mode = shape.drawing_mode;
		_has_indices = shape.ibo != 0u;
		_name = std::string("Render ") + shape.name;

		_constants = shape.material;
	}

	TRSTransformf& get_transform()
	{
		return _transform;
	}

	void set_material_constants(bonobo::material_data const& constants)
	{
		_constants = constants;
	}
};

edaf80::Fluid::Fluid(WindowManager& windowManager) :
	mCamera(0.5f * glm::half_pi<float>(),
		static_cast<float>(config::resolution_x) / static_cast<float>(config::resolution_y),
		0.01f, 1000.0f),
	inputHandler(), mWindowManager(windowManager), window(nullptr)
{
	WindowManager::WindowDatum window_datum{ inputHandler, mCamera, config::resolution_x, config::resolution_y, 0, 0, 0, 0 };

	window = mWindowManager.CreateGLFWWindow("EDAF80: Fluid", window_datum, config::msaa_rate);
	if (window == nullptr) {
		throw std::runtime_error("Failed to get a window: aborting!");
	}

	bonobo::init();
}

edaf80::Fluid::~Fluid()
{
	bonobo::deinit();
}

void
edaf80::Fluid::run()
{
	// Set up the camera
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, 90.0f));
	mCamera.mMouseSensitivity = glm::vec2(0.003f);
	mCamera.mMovementSpeed = glm::vec3(3.0f); // 3 m/s => 10.8 km/h

	// Create the shader programs
	ShaderProgramManager program_manager;
	GLuint fallback_shader = 0u;
	program_manager.CreateAndRegisterProgram("Fallback",
		{ { ShaderType::vertex, "common/fallback.vert" },
		  { ShaderType::fragment, "common/fallback.frag" } },
		fallback_shader);
	if (fallback_shader == 0u) {
		LogError("Failed to load fallback shader");
		return;
	}

	GLuint diffuse_shader = 0u;
	program_manager.CreateAndRegisterProgram("Diffuse",
		{ { ShaderType::vertex, "sphere_shaders/diffuse.vert" },
		  { ShaderType::fragment, "sphere_shaders/diffuse.frag" } },
		diffuse_shader);
	if (diffuse_shader == 0u) {
		LogError("Failed to load diffuse shader");
		return;
	}


	GLuint predicted_position_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Predicted positions", "compute_shaders/predicted.comp", predicted_position_shader);
	if (predicted_position_shader == 0u) {
		LogError("Failed to load predicted position shader");
		return;
	}

	GLuint spatial1_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Spatial 1", "compute_shaders/spatial1.comp", spatial1_shader);
	if (spatial1_shader == 0u) {
		LogError("Failed to load spatial 1 shader");
		return;
	}

	GLuint bitonic_sort_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Bitonic sort", "compute_shaders/bitonic_sort.comp", bitonic_sort_shader);
	if (bitonic_sort_shader == 0u) {
		LogError("Failed to load bitonic sort shader");
		return;
	}

	//
	// Todo: Insert the creation of other shader programs.
	//       (Check how it was done in assignment 3.)
	//

	//
	// Todo: Load your geometry
	//




	auto camera_position = mCamera.mWorld.GetTranslation();
	mCamera.mWorld.LookTowards(glm::vec3(0, 0, 1));
	float elapsed_time_s = 0.0f;
	auto light_position = glm::vec3(-20.0f, 40.0f, -120.0f);

	auto const set_uniforms = [&light_position, &elapsed_time_s, &camera_position](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform1f(glGetUniformLocation(program, "elapsed_time_s"), elapsed_time_s);
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
	};







	
	float grid_sphere_radius = 0.03;
	auto grid_sphere = parametric_shapes::createSphere(grid_sphere_radius, 2u, 2u);
	if (grid_sphere.vao == 0u) {
		LogError("Failed to retrieve the mesh for the grid_sphere");
		return;
	}

	float PI = 3.14159265358979;

	int num_particles = 16384;
	int num_work_groups = num_particles / 1024;
	int sqrtN = sqrt(num_particles);
	float time_step = 1 / 60.0;
	float damping_factor = 0.95;
	float width = 16.0f, height = 9.0f;
	float half_width = width / 2, half_height = height / 2;
	float gravity_strength = 10.0;
	float smoothing_radius = 0.272;//1.2;
	float target_density = 77.5;
	float pressure_multiplier = 64.0;
	float near_pressure_multiplier = 2.0;
	float viscosity_strength = 0.1;

	int PRIME1 = 86183;
	int PRIME2 = 7475723;

	std::cout << "Max local size: " << GL_MAX_COMPUTE_WORK_GROUP_SIZE << std::endl;
	std::cout << "Max local size product: " << GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS << std::endl;

	std::cout << "X: " << glGetIntegeri_v << std::endl;
	std::cout << "X: " << GL_MAX_COMPUTE_WORK_GROUP_COUNT  << std::endl;


	volatile float avg_density = 0;

	std::vector<Node> nodes;
	std::vector<glm::vec3> positions;
	//std::vector<glm::vec3> velocities;
	std::vector<glm::vec3> colours;
	std::vector<Particle> particles;

	ParticleRenderer particle_renderer;
	particle_renderer.set_geometry(grid_sphere);
	particle_renderer.set_program(&diffuse_shader, set_uniforms);
	bonobo::material_data material;
	material.diffuse = glm::vec3(0.0, 0.6, 1.0);
	particle_renderer.set_material_constants(material);

	//float spacing = 3 * grid_sphere_radius;
	glm::vec3 square_center = glm::vec3(width / 10, height / 10, 0.0);
	float spacing = half_height / ((sqrtN - 1));
	for (int i = 0; i < sqrtN; i++) {
		for (int j = 0; j < sqrtN; j++) {
			auto pos = glm::vec3(-half_height / 2 + i * spacing, -half_height / 2 + j * spacing, 100.0f) + square_center;
			//auto pos = glm::vec3((rand() * 1.0f / RAND_MAX) * width - half_width, (rand() * 1.0f / RAND_MAX) * height - half_height, 100.0f);
			positions.push_back(pos);
			particles.push_back(Particle());
			particles.back().position = pos;

			auto velocity = glm::vec3(0, 0, 0);
			//velocities.push_back(velocity);
			colours.push_back(glm::vec3(0.0, 0.0, 1.0));

			Node node;
			node.set_geometry(grid_sphere);
			node.set_program(&diffuse_shader, set_uniforms);
			node.get_transform().SetTranslate(pos);
			bonobo::material_data material;
			material.diffuse = glm::vec3(0.0, 0.6, 1.0);
			node.set_material_constants(material);
			nodes.push_back(node);
		}
	}


	num_particles = particles.size();
	std::vector<std::pair<float, float>> densities(num_particles);
	std::vector<glm::vec3> predicted_positions(num_particles);

	auto handle_collision = [&](int idx) -> void {
		if (particles[idx].position.x + grid_sphere_radius > half_width) {
			particles[idx].position.x = half_width - grid_sphere_radius;
			particles[idx].velocity.x *= -1 * damping_factor;
		}
		if (particles[idx].position.x - grid_sphere_radius < -half_width) {
			particles[idx].position.x = -half_width + grid_sphere_radius;
			particles[idx].velocity.x *= -1 * damping_factor;
		}
		if (particles[idx].position.y + grid_sphere_radius > half_height) {
			particles[idx].position.y = half_height - grid_sphere_radius;
			particles[idx].velocity.y *= -1 * damping_factor;
		}
		if (particles[idx].position.y - grid_sphere_radius < -half_height) {
			particles[idx].position.y = -half_height + grid_sphere_radius;
			particles[idx].velocity.y *= -1 * damping_factor;
		}
	};

	std::vector<glm::ivec2> spatial(num_particles);
	std::vector<int> start_inds(num_particles);

	auto good_mod = [](int a, int m) -> int {
		if (a >= 0) return a % m;
		return (-(-a % m) + m) % m;
	};

	auto point_to_hash = [&](glm::vec3 point) -> int {
		int x = (point.x / smoothing_radius);
		int y = (point.y / smoothing_radius);
		int hash = x * PRIME1 + y * PRIME2;
		return good_mod(hash, num_particles);
	};

	auto update_spatial = [&]() -> void {
		concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			spatial[i] = { point_to_hash(particles[i].predicted_position), i };
		});

		sort(spatial.begin(), spatial.end(), [](glm::ivec2 a, glm::ivec2 b) {return a.x < b.x || (a.x == b.x && a.y < b.y);});


		std::fill(start_inds.begin(), start_inds.end(), -1);
		start_inds[spatial[0].x] = 0;
		concurrency::parallel_for(size_t(1), size_t(num_particles), [&](size_t i) {
			if (spatial[i].x != spatial[i - 1].x) {
				start_inds[spatial[i].x] = i;
			}
		});
	};


	auto smoothing_kernel = [&](float dist, float r) -> float {
		/*float volume = PI * std::pow(r, 8.0) / 4.0;
		float value = std::max(0.0f, r * r - dist * dist);
		return value * value * value / volume;*/

		if (dist >= r) return 0;
		float volume = (PI * std::pow(r, 4)) / 6;
		float x = r - dist;
		return x * x / volume;
	};

	auto smoothing_kernel_derivative = [&](float dist, float r) -> float {
		/*if (dist >= r) return 0;
		float f = r * r - dist * dist;
		float scale = -24 / (PI * std::pow(r, 8.0));
		return scale * dist * f * f;*/
		if (dist >= r) return 0;
		float scale = 12 / (std::pow(r, 4) * PI);
		return (dist - r) * scale;
	};

	auto viscosity_smoothing_kernel = [&](float dist, float r) -> float {
		float volume = PI * std::pow(r, 8.0) / 4.0;
		float value = std::max(0.0f, r * r - dist * dist);
		return value * value * value / volume;
	};

	auto near_density_kernel = [&](float dist, float r) -> float {
		if (dist >= r) return 0;
		float v = r - dist;
		float volume = PI * std::pow(r, 5) / 10;
		return v * v * v / volume;
	};

	auto near_density_derivative = [&](float dist, float r) -> float {
		if (dist >= r) return 0;
		float v = r - dist;
		float volume = PI * std::pow(r, 5) / 30;
		return -v * v / volume;
	};


	const float mass = 1.0;
	auto calculate_density = [&](glm::vec3 point) -> std::pair<float, float> {

		float density = 0.0;
		float near_density = 0.0;

		int hash = point_to_hash(point);

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int cur_hash = hash + i * PRIME1 + j * PRIME2;
				cur_hash = good_mod(cur_hash, num_particles);
				int start_spatial_ind = start_inds[cur_hash];
				if (start_spatial_ind == -1) continue;
				for (int spatial_ind = start_spatial_ind; spatial_ind < num_particles && spatial[spatial_ind].x == cur_hash; spatial_ind++) {
					//float dist = glm::l2Norm(positions[spatial[spatial_ind].second] - point);
					float dist = glm::l2Norm(particles[spatial[spatial_ind].y].predicted_position - point);
					float influence = smoothing_kernel(dist, smoothing_radius);
					density += mass * influence;

					near_density += mass * near_density_kernel(dist, smoothing_radius);
				}
			}
		}

		/*for (int i = 0; i < num_particles; i++) {
			float dist = glm::l2Norm(positions[i] - point);
			float influence = smoothing_kernel(dist, smoothing_radius);
			density += mass * influence;
		}*/

		return { density, near_density };
	};


	auto density_to_pressure = [&](float density) -> float {
		float density_error = density - target_density;
		float pressure = density_error * pressure_multiplier;
		return pressure;
	};

	auto near_density_to_pressure = [&](float near_density) -> float {
		return near_density * near_pressure_multiplier;
	};

	auto get_random_dir = [&]() -> glm::vec3 {
		float theta = (rand() * 1.0f / RAND_MAX) * 2.0f * PI;
		return glm::vec3(glm::cos(theta), glm::sin(theta), 0.0f);
	};

	auto calculcate_shared_pressure = [&](float a, float b) -> float {
		return (density_to_pressure(a) + density_to_pressure(b)) * 0.5;
	};

	auto calculate_pressure_force = [&](int idx) -> glm::vec3 {
		glm::vec3 force = glm::vec3(0.0);
		//glm::vec3 point = positions[idx];
		glm::vec3 point = particles[idx].predicted_position;
		float pressure = density_to_pressure(particles[idx].density);
		float near_pressure = near_density_to_pressure(particles[idx].near_density);

		int hash = point_to_hash(point);

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int cur_hash = hash + i * PRIME1 + j * PRIME2;
				cur_hash = good_mod(cur_hash, num_particles);
				int start_spatial_ind = start_inds[cur_hash];
				if (start_spatial_ind == -1) continue;

				for (int spatial_ind = start_spatial_ind; spatial_ind < num_particles && spatial[spatial_ind].x == cur_hash; spatial_ind++) {
					int other_idx = spatial[spatial_ind].y;
					if (other_idx == idx) continue;

					//glm::vec3 offset = positions[other_idx] - point;
					glm::vec3 offset = particles[other_idx].predicted_position - point;
					float dist = glm::l2Norm(offset);
					glm::vec3 dir = dist == 0 ? get_random_dir() : offset / dist;
					float slope = smoothing_kernel_derivative(dist, smoothing_radius);
					float other_density = particles[other_idx].density;
					float other_near_density = particles[other_idx].near_density;

					float other_pressure = density_to_pressure(other_density);
					float other_near_pressure = near_density_to_pressure(other_near_density);

					float shared_pressure = (pressure + other_pressure) * 0.5;
					float shared_near_pressure = (near_pressure + other_near_pressure) * 0.5;
					force += shared_pressure * dir * smoothing_kernel_derivative(dist, smoothing_radius) * mass / other_density;
					force += shared_near_pressure * dir * near_density_derivative(dist, smoothing_radius) * mass / other_near_density;
				}
			}
		}



		/*for (int i = 0; i < num_particles; i++) {
			if (i == idx) continue;

			glm::vec3 offset = positions[i] - point;
			float dist = glm::l2Norm(offset);
			glm::vec3 dir = dist == 0 ? get_random_dir() : offset / dist;
			float slope = smoothing_kernel_derivative(dist, smoothing_radius);
			float density = densities[i];
			float shared_pressure = calculcate_shared_pressure(density, densities[idx]);
			force += shared_pressure * dir * slope * mass / density;
		}*/

		return force;
	};


	auto calculate_viscosity_force = [&](int idx) -> glm::vec3 {
		glm::vec3 force = glm::vec3(0.0);
		glm::vec3 point = particles[idx].predicted_position;

		int hash = point_to_hash(point);

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int cur_hash = hash + i * PRIME1 + j * PRIME2;
				cur_hash = good_mod(cur_hash, num_particles);
				int start_spatial_ind = start_inds[cur_hash];
				if (start_spatial_ind == -1) continue;

				for (int spatial_ind = start_spatial_ind; spatial_ind < num_particles && spatial[spatial_ind].x == cur_hash; spatial_ind++) {
					int other_idx = spatial[spatial_ind].y;

					//glm::vec3 offset = positions[other_idx] - point;
					glm::vec3 offset = particles[other_idx].predicted_position - point;
					float dist = glm::l2Norm(offset);
					float influence = viscosity_smoothing_kernel(dist, smoothing_radius);
					force += (particles[other_idx].velocity - particles[idx].velocity) * influence;
				}
			}
		}

		return force * viscosity_strength;
	};


	auto simulation_step = [&](float delta_time) -> void {

		GLuint ssbo_particles;
		glGenBuffers(1, &ssbo_particles);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_particles);
		glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo_particles);

		GLuint ssbo_spatial;
		glGenBuffers(1, &ssbo_spatial);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_spatial);
		glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(glm::vec2), spatial.data(), GL_DYNAMIC_DRAW);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo_spatial);

		GLuint ssbo_start_inds;
		glGenBuffers(1, &ssbo_start_inds);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_start_inds);
		glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(int), start_inds.data(), GL_DYNAMIC_DRAW);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssbo_start_inds);


		
		/*concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			particles[i].velocity += glm::vec3(0.0, -1.0, 0.0) * gravity_strength * delta_time;
			particles[i].predicted_position = particles[i].position + particles[i].velocity * delta_time;
		});*/
		glUseProgram(predicted_position_shader);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
		////////////////////////////////////////////////////////////



		/*concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			spatial[i] = { point_to_hash(particles[i].predicted_position), i };
		});*/
		glUseProgram(spatial1_shader);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
		////////////////////////////////////////////////////////////



		/*for (int k = 2; k <= num_particles; k *= 2) {
			for (int j = k / 2; j > 0; j /= 2) {
				for (int i = 0; i < num_particles; i++) {
					int l = i ^ j;
					if (l > i) {
						if ((i & k) == 0 && spatial[i].x > spatial[l].x || (i & k) != 0 && spatial[i].x < spatial[l].x) {
							std::swap(spatial[i], spatial[l]);
						}
					}
				}
			}
		}*/
		//sort(spatial.begin(), spatial.end(), [](glm::ivec2 a, glm::ivec2 b) {return a.x < b.x || (a.x == b.x && a.y < b.y);});
		glUseProgram(bitonic_sort_shader);
		for (int k = 2; k <= num_particles; k *= 2) {
			for (int j = k / 2; j > 0; j /= 2) {
				GLint kk = k;
				GLint jj = j;
				glUniform1i(glGetUniformLocation(bitonic_sort_shader, "k"), kk);
				glUniform1i(glGetUniformLocation(bitonic_sort_shader, "j"), jj);
				glDispatchCompute(num_work_groups, 1, 1);
				glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
			}
		}

		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_particles);
		glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_particles * sizeof(Particle), particles.data());
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_start_inds);
		glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_particles * sizeof(int), start_inds.data());
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_spatial);
		glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_particles * sizeof(glm::vec2), spatial.data());



		std::fill(start_inds.begin(), start_inds.end(), -1);
		start_inds[spatial[0].x] = 0;
		concurrency::parallel_for(size_t(1), size_t(num_particles), [&](size_t i) {
		for(int i = 1; i < num_particles; i++)
			if (spatial[i].x != spatial[i - 1].x) {
				start_inds[spatial[i].x] = i;
			}
		});


		//update_spatial();





		concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			std::pair<float, float> dens = calculate_density(particles[i].predicted_position);
			particles[i].density = dens.first;
			particles[i].near_density = dens.second;
		});

		concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			glm::vec3 pressure_force = calculate_pressure_force(i);
			glm::vec3 pressure_acceleration = pressure_force / particles[i].density;
			particles[i].velocity += pressure_acceleration * delta_time;
		});
		
		concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			glm::vec3 viscosity_force = calculate_viscosity_force(i);
			particles[i].velocity += viscosity_force * delta_time;
		});

		concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
			particles[i].position += particles[i].velocity * delta_time;
			handle_collision(i);
		});
	};



	auto speed_to_color = [](float speed) -> glm::vec3 {
		bonobo::material_data mat;
		if (speed > 17) {
			return glm::vec3(1.0, 0.0, 0.0);
		}
		else if (speed > 12) {
			return glm::vec3(0.9, 0.3, 0.0);
		}
		else if (speed > 8) {
			return glm::vec3(0.7, 0.4, 0.0);
		}
		else if (speed > 5) {
			return glm::vec3(0.5, 0.5, 0.0);
		}
		else if (speed > 2) {
			return glm::vec3(0.0, 0.5, 0.5);
		}
		else {
			return glm::vec3(0.0, 0.6, 1.0);
		}
	};

	glClearDepthf(1.0f);
	glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	glEnable(GL_DEPTH_TEST);


	auto lastTime = std::chrono::high_resolution_clock::now();

	bool show_logs = true;
	bool show_gui = true;
	bool shader_reload_failed = false;
	bool show_basis = false;
	float basis_thickness_scale = 1.0f;
	float basis_length_scale = 1.0f;

	while (!glfwWindowShouldClose(window)) {
		half_width = width / 2;
		half_height = height / 2;

		auto const nowTime = std::chrono::high_resolution_clock::now();
		auto const deltaTimeUs = std::chrono::duration_cast<std::chrono::microseconds>(nowTime - lastTime);
		lastTime = nowTime;

		auto& io = ImGui::GetIO();
		inputHandler.SetUICapture(io.WantCaptureMouse, io.WantCaptureKeyboard);

		glfwPollEvents();
		inputHandler.Advance();
		mCamera.Update(deltaTimeUs, inputHandler);

		if (inputHandler.GetKeycodeState(GLFW_KEY_R) & JUST_PRESSED) {
			shader_reload_failed = !program_manager.ReloadAllPrograms();
			if (shader_reload_failed)
				tinyfd_notifyPopup("Shader Program Reload Error",
					"An error occurred while reloading shader programs; see the logs for details.\n"
					"Rendering is suspended until the issue is solved. Once fixed, just reload the shaders again.",
					"error");
		}
		if (inputHandler.GetKeycodeState(GLFW_KEY_F3) & JUST_RELEASED)
			show_logs = !show_logs;
		if (inputHandler.GetKeycodeState(GLFW_KEY_F2) & JUST_RELEASED)
			show_gui = !show_gui;
		if (inputHandler.GetKeycodeState(GLFW_KEY_F11) & JUST_RELEASED)
			mWindowManager.ToggleFullscreenStatusForWindow(window);


		// Retrieve the actual framebuffer size: for HiDPI monitors,
		// you might end up with a framebuffer larger than what you
		// actually asked for. For example, if you ask for a 1920x1080
		// framebuffer, you might get a 3840x2160 one instead.
		// Also it might change as the user drags the window between
		// monitors with different DPIs, or if the fullscreen status is
		// being toggled.
		int framebuffer_width, framebuffer_height;
		glfwGetFramebufferSize(window, &framebuffer_width, &framebuffer_height);
		glViewport(0, 0, framebuffer_width, framebuffer_height);


		//
		// Todo: If you need to handle inputs, you can do it here
		//


		mWindowManager.NewImGuiFrame();

		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

		float time_to_step;

		// time_step = std::chrono::duration<float>(deltaTimeUs).count();
		if (!shader_reload_failed) {
			//
			// Todo: Render all your geometry here.
			//
			/*for (int i = 0; i < num_particles; i++) {
				velocities[i] += glm::vec3(0.0, -gravity_strength, 0.0) * time_step;
				predicted_positions[i] = positions[i] + velocities[i] * (1 / 120.0f);
			}*/

			auto t0 = std::chrono::high_resolution_clock::now();
			simulation_step(time_step);
			auto t1 = std::chrono::high_resolution_clock::now();
			time_to_step = std::chrono::duration<float>(t1 - t0).count();


			concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
				colours[i] = speed_to_color(glm::l2Norm(particles[i].velocity));
				positions[i] = particles[i].position;
				/*nodes[i].get_transform().SetTranslate(positions[i]);
				nodes[i].render(mCamera.GetWorldToClipMatrix());
				nodes[i].set_material_constants(speed_to_color(glm::l2Norm(velocities[i])));*/
			});

			particle_renderer.render(mCamera.GetWorldToClipMatrix(), positions, colours);
		}


		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//
		// Todo: If you want a custom ImGUI window, you can set it up
		//       here
		//
		bool const opened = ImGui::Begin("Scene Controls", nullptr, ImGuiWindowFlags_None);
		if (opened) {
			ImGui::Text("Avg density: %.4f", avg_density / num_particles);
			ImGui::Text("Time per frame: %.10f", std::chrono::duration<float>(deltaTimeUs).count());
			ImGui::Text("Time per update: %.10f", time_to_step);

			ImGui::Checkbox("Show basis", &show_basis);
			ImGui::SliderFloat("Basis thickness scale", &basis_thickness_scale, 0.0f, 100.0f);
			ImGui::SliderFloat("Basis length scale", &basis_length_scale, 0.0f, 100.0f);

			ImGui::SliderFloat("Bounding box width", &width, 10, 100);
			ImGui::SliderFloat("Bounding box height", &height, 10, 100);

			ImGui::SliderFloat("Damping factor", &damping_factor, 0.0f, 1.0f);
			ImGui::SliderFloat("Gravity strength", &gravity_strength, 0.0f, 100.0f);
			ImGui::SliderFloat("Smoothing radius", &smoothing_radius, grid_sphere_radius, 2.0f);
			//ImGui::SliderFloat("Time step", &time_step, 0.0f, 0.5f);
			ImGui::SliderFloat("Target_density", &target_density, 0.0f, 700.0f);
			ImGui::SliderFloat("Pressure multiplier", &pressure_multiplier, 0.0f, 1000.0f);
			ImGui::SliderFloat("Viscosity strength", &viscosity_strength, 0.0f, 1.0f);
			ImGui::SliderFloat("Near pressure multiplier", &near_pressure_multiplier, 0.0f, 10.0f);
		}
		ImGui::End();

		if (show_basis)
			bonobo::renderBasis(basis_thickness_scale, basis_length_scale, mCamera.GetWorldToClipMatrix());
		if (show_logs)
			Log::View::Render();
		mWindowManager.RenderImGuiFrame(show_gui);

		glfwSwapBuffers(window);
	}
}

int main()
{
	std::setlocale(LC_ALL, "");

	Bonobo framework;

	try {
		edaf80::Fluid assignment5(framework.GetWindowManager());
		assignment5.run();
	}
	catch (std::runtime_error const& e) {
		LogError(e.what());
	}
}
