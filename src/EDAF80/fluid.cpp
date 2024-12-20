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
};

struct ViscosityForce {
	glm::vec3 force;
	float padding;
};

struct Colour {
	glm::vec3 colour;
	float padding;
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
	std::vector<std::tuple<std::string, GLuint, GLenum>> _textures;
	bonobo::material_data _constants;

	// Transformation data
	TRSTransformf _transform;

	// Debug data
	std::string _name{ "Render un-named node" };

	void render(glm::mat4 const& view_projection, glm::mat4 const& view_projection_inv, glm::mat4& world, GLuint program, std::function<void(GLuint)> const& set_uniforms, int num_particles) const {
		if (_vao == 0u || program == 0u)
			return;

		utils::opengl::debug::beginDebugGroup(_name);

		glUseProgram(program);

		auto const normal_model_to_world = glm::transpose(glm::inverse(world));

		set_uniforms(program);

		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_model_to_world"), 1, GL_FALSE, glm::value_ptr(world));
		glUniformMatrix4fv(glGetUniformLocation(program, "normal_model_to_world"), 1, GL_FALSE, glm::value_ptr(normal_model_to_world));
		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_world_to_clip"), 1, GL_FALSE, glm::value_ptr(view_projection));
		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_clip_to_world"), 1, GL_FALSE, glm::value_ptr(view_projection_inv));
		
		for (size_t i = 0u; i < _textures.size(); ++i) {
			auto const& texture = _textures[i];
			glActiveTexture(GL_TEXTURE0 + static_cast<GLenum>(i));
			glBindTexture(std::get<2>(texture), std::get<1>(texture));
			glUniform1i(glGetUniformLocation(program, std::get<0>(texture).c_str()), static_cast<GLint>(i));

			std::string texture_presence_var_name = "has_" + std::get<0>(texture);
			glUniform1i(glGetUniformLocation(program, texture_presence_var_name.c_str()), 1);
		}

		glUniform3fv(glGetUniformLocation(program, "diffuse_colour"), 1, glm::value_ptr(_constants.diffuse));
		glUniform3fv(glGetUniformLocation(program, "specular_colour"), 1, glm::value_ptr(_constants.specular));
		glUniform3fv(glGetUniformLocation(program, "ambient_colour"), 1, glm::value_ptr(_constants.ambient));
		glUniform3fv(glGetUniformLocation(program, "emissive_colour"), 1, glm::value_ptr(_constants.emissive));
		glUniform1f(glGetUniformLocation(program, "shininess_value"), _constants.shininess);
		glUniform1f(glGetUniformLocation(program, "index_of_refraction_value"), _constants.indexOfRefraction);
		glUniform1f(glGetUniformLocation(program, "opacity_value"), _constants.opacity);



		glBindVertexArray(_vao);

		if (_has_indices)
			glDrawElementsInstanced(_drawing_mode, _indices_nb, GL_UNSIGNED_INT, reinterpret_cast<GLvoid const*>(0x0), num_particles);
		else
			glDrawArrays(_drawing_mode, 0, _vertices_nb);
		glBindVertexArray(0u);

		for (auto const& texture : _textures) {
			glBindTexture(std::get<2>(texture), 0);
			glUniform1i(glGetUniformLocation(program, std::get<0>(texture).c_str()), 0);

			std::string texture_presence_var_name = "has_" + std::get<0>(texture);
			glUniform1i(glGetUniformLocation(program, texture_presence_var_name.c_str()), 0);
		}

		glUseProgram(0u);

		utils::opengl::debug::endDebugGroup();
	}

	void render2(glm::mat4 const& view_projection, glm::mat4 const& view_projection_inv, glm::mat4& world, GLuint program, std::function<void(GLuint)> const& set_uniforms) const {
		if (_vao == 0u || program == 0u)
			return;

		utils::opengl::debug::beginDebugGroup(_name);

		glUseProgram(program);

		auto const normal_model_to_world = glm::transpose(glm::inverse(world));

		set_uniforms(program);

		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_model_to_world"), 1, GL_FALSE, glm::value_ptr(world));
		glUniformMatrix4fv(glGetUniformLocation(program, "normal_model_to_world"), 1, GL_FALSE, glm::value_ptr(normal_model_to_world));
		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_world_to_clip"), 1, GL_FALSE, glm::value_ptr(view_projection));
		glUniformMatrix4fv(glGetUniformLocation(program, "vertex_clip_to_world"), 1, GL_FALSE, glm::value_ptr(view_projection_inv));

		for (size_t i = 0u; i < _textures.size(); ++i) {
			auto const& texture = _textures[i];
			glActiveTexture(GL_TEXTURE0 + static_cast<GLenum>(i));
			glBindTexture(std::get<2>(texture), std::get<1>(texture));
			glUniform1i(glGetUniformLocation(program, std::get<0>(texture).c_str()), static_cast<GLint>(i));

			std::string texture_presence_var_name = "has_" + std::get<0>(texture);
			glUniform1i(glGetUniformLocation(program, texture_presence_var_name.c_str()), 1);
		}


		glUniform3fv(glGetUniformLocation(program, "diffuse_colour"), 1, glm::value_ptr(_constants.diffuse));
		glUniform3fv(glGetUniformLocation(program, "specular_colour"), 1, glm::value_ptr(_constants.specular));
		glUniform3fv(glGetUniformLocation(program, "ambient_colour"), 1, glm::value_ptr(_constants.ambient));
		glUniform3fv(glGetUniformLocation(program, "emissive_colour"), 1, glm::value_ptr(_constants.emissive));
		glUniform1f(glGetUniformLocation(program, "shininess_value"), _constants.shininess);
		glUniform1f(glGetUniformLocation(program, "index_of_refraction_value"), _constants.indexOfRefraction);
		glUniform1f(glGetUniformLocation(program, "opacity_value"), _constants.opacity);



		glBindVertexArray(_vao);

		if (_has_indices)
			glDrawElements(_drawing_mode, _indices_nb, GL_UNSIGNED_INT, reinterpret_cast<GLvoid const*>(0x0));
		else
			glDrawArrays(_drawing_mode, 0, _vertices_nb);
		glBindVertexArray(0u);


		for (auto const& texture : _textures) {
			glBindTexture(std::get<2>(texture), 0);
			glUniform1i(glGetUniformLocation(program, std::get<0>(texture).c_str()), 0);

			std::string texture_presence_var_name = "has_" + std::get<0>(texture);
			glUniform1i(glGetUniformLocation(program, texture_presence_var_name.c_str()), 0);
		}

		glUseProgram(0u);

		utils::opengl::debug::endDebugGroup();
	}

public:
	void render(glm::mat4 const& view_projection, glm::mat4 const& view_projection_inv, int num_particles) const {
		if (_program != nullptr) {
			render(view_projection, view_projection_inv, _transform.GetMatrix(), *_program, _set_uniforms, num_particles);
		}
	}

	void render2(glm::mat4 const& view_projection, glm::mat4 const& view_projection_inv) const {
		if (_program != nullptr) {
			render2(view_projection, view_projection_inv, _transform.GetMatrix(), *_program, _set_uniforms);
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

	void add_texture(std::string const& name, GLuint tex_id, GLenum type)
	{
		GLint max_combined_texture_image_units{ -1 };
		glGetIntegerv(GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS, &max_combined_texture_image_units);
		std::size_t const max_active_texture_count
			= (max_combined_texture_image_units > 0) ? static_cast<std::size_t>(max_combined_texture_image_units)
			: 80; // OpenGL 4.x guarantees at least 80.

		if (_textures.size() >= max_active_texture_count) {
			LogWarning("Trying to add more textures to an object than supported (%llu); the texture %s with ID %u will **not** be added. If you really need that many textures, do not use the `Node` class and roll your own solution instead.",
				max_active_texture_count, name.c_str(), tex_id);
			return;
		}
		if (tex_id == 0u) {
			LogWarning("0 is not a valid texture ID; the texture %s (with ID %u) will **not** be added.",
				name.c_str(), tex_id);
			return;
		}

		_textures.emplace_back(name, tex_id, type);
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
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, -25.0f));
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

	GLuint ray_march_shader = 0u;
	program_manager.CreateAndRegisterProgram("Ray March",
		{ { ShaderType::vertex, "sphere_shaders/ray_march.vert" },
		  { ShaderType::fragment, "sphere_shaders/ray_march.frag" } },
		ray_march_shader);
	if (ray_march_shader == 0u) {
		LogError("Failed to load ray_march shader");
		return;
	}


	GLuint skybox_shader = 0u;
	program_manager.CreateAndRegisterProgram("Skybox",
		{ { ShaderType::vertex, "sphere_shaders/skybox.vert" },
		  { ShaderType::fragment, "sphere_shaders/skybox.frag" } },
		skybox_shader);
	if (skybox_shader == 0u)
		LogError("Failed to load skybox shader");


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

	GLuint spatial2_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Spatial 2", "compute_shaders/spatial2.comp", spatial2_shader);
	if (spatial2_shader == 0u) {
		LogError("Failed to load spatial 2 shader");
		return;
	}

	GLuint density_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Density", "compute_shaders/density.comp", density_shader);
	if (density_shader == 0u) {
		LogError("Failed to load density shader");
		return;
	}

	GLuint pressure_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Pressure", "compute_shaders/pressure.comp", pressure_shader);
	if (pressure_shader == 0u) {
		LogError("Failed to load pressure shader");
		return;
	}

	GLuint viscosity_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Viscosity", "compute_shaders/viscosity.comp", viscosity_shader);
	if (viscosity_shader == 0u) {
		LogError("Failed to load viscosity shader");
		return;
	}

	GLuint viscosity_velocity_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Viscosity velocity", "compute_shaders/viscosity_velocity.comp", viscosity_velocity_shader);
	if (viscosity_velocity_shader == 0u) {
		LogError("Failed to load viscosity velocity shader");
		return;
	}

	GLuint update_position_shader = 0u;
	program_manager.CreateAndRegisterComputeProgram("Update position", "compute_shaders/update_position.comp", update_position_shader);
	if (update_position_shader == 0u) {
		LogError("Failed to load update_position shader");
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
	mCamera.mWorld.LookAt(glm::vec3(0, 5, 0), glm::vec3(0, 1, 0));
	float elapsed_time_s = 0.0f;
	auto light_position = glm::vec3(-20.0f, 40.0f, -20.0f);

	auto const set_uniforms = [&light_position, &elapsed_time_s, &camera_position](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform1f(glGetUniformLocation(program, "elapsed_time_s"), elapsed_time_s);
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
	};


	float grid_sphere_radius = 0.07;
	auto grid_sphere = parametric_shapes::createSphere(grid_sphere_radius, 15u, 15u);
	if (grid_sphere.vao == 0u) {
		LogError("Failed to retrieve the mesh for the grid_sphere");
		return;
	}


	auto grid_quad = parametric_shapes::createHighTesselationQuad(2.0, 2.0, 1, 1);
	if (grid_quad.vao == 0u) {
		LogError("Failed to retrieve the mesh for the grid quad");
		return;
	}


	auto skybox_shape = parametric_shapes::createSphere(500.0f, 100u, 100u);
	if (skybox_shape.vao == 0u) {
		LogError("Failed to retrieve the mesh for the skybox");
		return;
	}



	GLuint cubemap = bonobo::loadTextureCubeMap(
		config::resources_path("cubemaps/NissiBeach2/posx.jpg"),
		config::resources_path("cubemaps/NissiBeach2/negx.jpg"),
		config::resources_path("cubemaps/NissiBeach2/posy.jpg"),
		config::resources_path("cubemaps/NissiBeach2/negy.jpg"),
		config::resources_path("cubemaps/NissiBeach2/posz.jpg"),
		config::resources_path("cubemaps/NissiBeach2/negz.jpg"));


	GLuint cubemap2 = bonobo::loadTextureCubeMap(
		config::resources_path("cubemaps/NissiBeach2/posx.jpg"),
		config::resources_path("cubemaps/NissiBeach2/negx.jpg"),
		config::resources_path("cubemaps/NissiBeach2/posy.jpg"),
		config::resources_path("cubemaps/NissiBeach2/negy.jpg"),
		config::resources_path("cubemaps/NissiBeach2/posz.jpg"),
		config::resources_path("cubemaps/NissiBeach2/negz.jpg"));


	Node skybox;
	skybox.set_geometry(skybox_shape);
	skybox.set_program(&skybox_shader, set_uniforms);
	skybox.add_texture("cubemap", cubemap, GL_TEXTURE_CUBE_MAP);



	float PI = 3.14159265358979;

	int num_particles = 16384 * 16;
	int num_work_groups = num_particles / 1024;
	int sqrtN = sqrt(num_particles);
	int cbrtN = pow(num_particles, 0.35);
	float time_step = 1 / 60.0;
	float damping_factor = 0.95;
	float width = 16.0f, height = 9.0f, depth = 16.0f;
	float half_width = width / 2, half_height = height / 2, half_depth = depth / 2;;
	float gravity_strength = 10.0;
	float smoothing_radius = 0.272;//1.2;
	float target_density = 77.5;
	float pressure_multiplier = 64.0;
	float near_pressure_multiplier = 2.0;
	float viscosity_strength = 0.03;
	float step_size = 0.05;
	float density_multiplier = 0.01;
	int num_refractions = 1;
	glm::vec3 scattering_coefficients = glm::vec3(0.173, 0.069, 0.0);

	int PRIME1 = 86183;
	int PRIME2 = 7475723;
	int PRIME3 = 447409;

	//std::vector<glm::vec3> velocities;
	std::vector<Colour> colours;
	std::vector<Particle> particles;
	std::vector<ViscosityForce> viscosity_forces;

	int particles_created = 0;

	//float spacing = 3 * grid_sphere_radius;
	glm::vec3 square_center = glm::vec3(width / 10, height / 10, 0.0);
	//float spacing = half_height / ((sqrtN - 0));
	//for (int i = 0; i < sqrtN + 1; i++) {
	//	for (int j = 0; j < sqrtN + 1 && particles_created < num_particles; j++) {
	//		particles_created++;
	//		auto pos = glm::vec3(-half_height / 2 + i * spacing, -half_height / 2 + j * spacing, rand() * 1.0 / RAND_MAX) + square_center;
	//		//auto pos = glm::vec3((rand() * 1.0f / RAND_MAX) * width - half_width, (rand() * 1.0f / RAND_MAX) * height - half_height, 100.0f);
	//		particles.push_back(Particle());
	//		particles.back().position = pos;

	//		auto velocity = glm::vec3(0, 0, 0);
	//		//velocities.push_back(velocity);
	//		colours.push_back(Colour());
	//		viscosity_forces.push_back(ViscosityForce());
	//	}
	//}


	float spacing = half_height / cbrtN;
	for (int i = 0; i < cbrtN + 1; i++) {
		for (int j = 0; j < cbrtN + 1; j++) {
			for (int k = 0; k < cbrtN + 1 && particles_created < num_particles; k++) {
				particles_created++;
				auto pos = glm::vec3(-half_height / 2 + i * spacing, -half_height / 2 + j * spacing, -half_depth / 2 + k * spacing) + square_center;
				//auto pos = glm::vec3((rand() * 1.0f / RAND_MAX) * width - half_width, (rand() * 1.0f / RAND_MAX) * height - half_height, 100.0f);
				particles.push_back(Particle());
				particles.back().position = pos;

				auto velocity = glm::vec3(0, 0, 0);
				//velocities.push_back(velocity);
				colours.push_back(Colour());
				viscosity_forces.push_back(ViscosityForce());
			}
		}
	}


	num_particles = particles.size();
	std::cout << num_particles << std::endl;
	std::vector<float> densities(num_particles);
	std::vector<float> near_densities(num_particles);
	std::vector<glm::vec3> predicted_positions(num_particles);

	std::vector<glm::ivec2> spatial(num_particles);
	std::vector<int> start_inds(num_particles);

	auto const set_uniforms_quad = [&](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform1f(glGetUniformLocation(program, "elapsed_time_s"), elapsed_time_s);
		glUniform1f(glGetUniformLocation(program, "step_size"), step_size);
		glUniform1f(glGetUniformLocation(program, "smoothing_radius"), smoothing_radius);
		glUniform1i(glGetUniformLocation(program, "num_particles"), num_particles);
		glUniform1i(glGetUniformLocation(program, "num_refractions"), num_refractions);
		glUniform1f(glGetUniformLocation(program, "density_multiplier"), density_multiplier);
		glUniform1f(glGetUniformLocation(program, "half_width"), half_width);
		glUniform1f(glGetUniformLocation(program, "half_height"), half_height);
		glUniform1f(glGetUniformLocation(program, "half_depth"), half_depth);
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "scattering_coefficients"), 1, glm::value_ptr(scattering_coefficients));
	};


	GLuint ssbo_particles;
	glGenBuffers(1, &ssbo_particles);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_particles);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo_particles);

	GLuint ssbo_densities;
	glGenBuffers(1, &ssbo_densities);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_densities);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(float), densities.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo_densities);

	GLuint ssbo_near_densities;
	glGenBuffers(1, &ssbo_near_densities);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_near_densities);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(float), near_densities.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssbo_near_densities);

	GLuint ssbo_spatial;
	glGenBuffers(1, &ssbo_spatial);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_spatial);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(glm::vec2), spatial.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ssbo_spatial);

	GLuint ssbo_start_inds;
	glGenBuffers(1, &ssbo_start_inds);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_start_inds);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(int), start_inds.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, ssbo_start_inds);

	GLuint ssbo_viscosity_forces;
	glGenBuffers(1, &ssbo_viscosity_forces);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_viscosity_forces);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(ViscosityForce), viscosity_forces.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, ssbo_viscosity_forces);

	GLuint ssbo_colours;
	glGenBuffers(1, &ssbo_colours);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_colours);
	glBufferData(GL_SHADER_STORAGE_BUFFER, num_particles * sizeof(Colour), colours.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, ssbo_colours);

	std::vector<Colour> speedmap(3);
	speedmap[0].colour = glm::vec3(0.0, 0.0, 1.0);
	speedmap[1].colour = glm::vec3(1.0, 1.0, 0.0);
	speedmap[2].colour = glm::vec3(1.0, 0.0, 0.0);

	GLuint ssbo_speedmap;
	glGenBuffers(1, &ssbo_speedmap);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_speedmap);
	glBufferData(GL_SHADER_STORAGE_BUFFER, speedmap.size() * sizeof(Colour), speedmap.data(), GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, ssbo_speedmap);





	ParticleRenderer particles_renderer;
	particles_renderer.set_geometry(grid_sphere);
	particles_renderer.set_program(&diffuse_shader, set_uniforms);

	ParticleRenderer ray_marcher;
	ray_marcher.set_geometry(grid_quad);
	ray_marcher.set_program(&ray_march_shader, set_uniforms_quad);
	ray_marcher.add_texture("cubemap", cubemap2, GL_TEXTURE_CUBE_MAP);




	auto simulation_step = [&](float delta_time) -> void {
		// Gravity and update predicted positions
		glUseProgram(predicted_position_shader);
		glUniform1f(glGetUniformLocation(predicted_position_shader, "gravity_strength"), gravity_strength);
		glUniform1f(glGetUniformLocation(predicted_position_shader, "delta_time"), delta_time);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);



		// Spatial hashing step 1
		glUseProgram(spatial1_shader);
		glUniform1i(glGetUniformLocation(spatial1_shader, "num_particles"), num_particles);
		glUniform1f(glGetUniformLocation(spatial1_shader, "smoothing_radius"), smoothing_radius);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);


		// Bitonic sort
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


		// Spatial hashing last step
		glUseProgram(spatial2_shader);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);



		// Density
		glUseProgram(density_shader);
		glUniform1i(glGetUniformLocation(density_shader, "num_particles"), num_particles);
		glUniform1f(glGetUniformLocation(density_shader, "smoothing_radius"), smoothing_radius);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);




		// Pressure
		glUseProgram(pressure_shader);
		glUniform1i(glGetUniformLocation(pressure_shader, "num_particles"), num_particles);
		glUniform1f(glGetUniformLocation(pressure_shader, "smoothing_radius"), smoothing_radius);
		glUniform1f(glGetUniformLocation(pressure_shader, "target_density"), target_density);
		glUniform1f(glGetUniformLocation(pressure_shader, "pressure_multiplier"), pressure_multiplier);
		glUniform1f(glGetUniformLocation(pressure_shader, "near_pressure_multiplier"), near_pressure_multiplier);
		glUniform1f(glGetUniformLocation(pressure_shader, "delta_time"), delta_time);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);



		// Viscosity
		glUseProgram(viscosity_shader);
		glUniform1i(glGetUniformLocation(viscosity_shader, "num_particles"), num_particles);
		glUniform1f(glGetUniformLocation(viscosity_shader, "smoothing_radius"), smoothing_radius);
		glUniform1f(glGetUniformLocation(viscosity_shader, "target_density"), target_density);
		glUniform1f(glGetUniformLocation(viscosity_shader, "viscosity_strength"), viscosity_strength);
		glUniform1f(glGetUniformLocation(viscosity_shader, "delta_time"), delta_time);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

		glUseProgram(viscosity_velocity_shader);
		glUniform1f(glGetUniformLocation(viscosity_velocity_shader, "delta_time"), delta_time);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);



		// Update position and check collision with wall
		glUseProgram(update_position_shader);
		glUniform1f(glGetUniformLocation(update_position_shader, "delta_time"), delta_time);
		glUniform1f(glGetUniformLocation(update_position_shader, "half_width"), half_width);
		glUniform1f(glGetUniformLocation(update_position_shader, "half_height"), half_height);
		glUniform1f(glGetUniformLocation(update_position_shader, "half_depth"), half_depth);
		glUniform1f(glGetUniformLocation(update_position_shader, "grid_sphere_radius"), grid_sphere_radius);
		glUniform1f(glGetUniformLocation(update_position_shader, "damping_factor"), damping_factor);
		glDispatchCompute(num_work_groups, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
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

	bool pause = true;
	bool particle = true;
	bool last_render_val = false;

	GLuint64 timer, timer2;
	GLint64 timer1;
	unsigned int query;
	int done = 0;

	simulation_step(time_step);

	while (!glfwWindowShouldClose(window)) {
		half_width = width / 2;
		half_height = height / 2;
		half_depth = depth / 2;

		

		auto const nowTime = std::chrono::high_resolution_clock::now();
		auto const deltaTimeUs = std::chrono::duration_cast<std::chrono::microseconds>(nowTime - lastTime);
		lastTime = nowTime;

		auto& io = ImGui::GetIO();
		inputHandler.SetUICapture(io.WantCaptureMouse, io.WantCaptureKeyboard);

		glfwPollEvents();
		inputHandler.Advance();
		mCamera.Update(deltaTimeUs, inputHandler);
		camera_position = mCamera.mWorld.GetTranslation();

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
			
			auto t0 = std::chrono::high_resolution_clock::now();
			if (!pause) {
				//glGenQueries(1, &query);
				//glBeginQuery(GL_TIME_ELAPSED, query);
				simulation_step(time_step);
				// if (!particle) simulation_step(time_step);
				//glEndQuery(GL_TIME_ELAPSED);
				//while (!done) {
				//	glGetQueryObjectiv(query,
				//		GL_QUERY_RESULT_AVAILABLE,
				//		&done);
				//}
				//glGetQueryObjectui64v(query, GL_QUERY_RESULT, &timer);
			}
			auto t1 = std::chrono::high_resolution_clock::now();
			time_to_step = std::chrono::duration<float>(t1 - t0).count();



			//if (particle != last_render_val) {
			//	if (particle) {
			//		//particle_renderer.add_texture("cubemap", cubemap, GL_TEXTURE_CUBE_MAP);
			//		particle_renderer.set_geometry(grid_sphere);
			//		particle_renderer.set_program(&diffuse_shader, set_uniforms);
			//	}
			//	else {
			//		//particle_renderer.add_texture("cubemap", cubemap, GL_TEXTURE_CUBE_MAP);
			//		particle_renderer.set_geometry(grid_quad);
			//		particle_renderer.set_program(&ray_march_shader, set_uniforms_quad);
			//	}
			//	last_render_val = particle;
			//}


			if (particle) {
				particles_renderer.render(mCamera.GetWorldToClipMatrix(), mCamera.GetClipToWorldMatrix(), num_particles);
				skybox.render(mCamera.GetWorldToClipMatrix());
			}
			else {
				ray_marcher.render2(mCamera.GetWorldToClipMatrix(), mCamera.GetClipToWorldMatrix());
			}

		}

		


		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//
		// Todo: If you want a custom ImGUI window, you can set it up
		//       here
		//
		bool const opened = ImGui::Begin("Scene Controls", nullptr, ImGuiWindowFlags_None);
		if (opened) {
			ImGui::Text("Time per frame:  %.10f", std::chrono::duration<float>(deltaTimeUs).count());
			// ImGui::Text("Time per update: %.10f", timer / 1e9);

			ImGui::Checkbox("Pause simulation", &pause);
			ImGui::Checkbox("Render particles", &particle);


			ImGui::SliderFloat("Coeff r", &scattering_coefficients.x, 0.0f, 10.0f);
			ImGui::SliderFloat("Coeff g", &scattering_coefficients.y, 0.0f, 10.0f);
			ImGui::SliderFloat("Coeff b", &scattering_coefficients.z, 0.0f, 10.0f);
			// ImGui::SliderInt("Number of refractions", &num_refractions, 0, 5);

			ImGui::Checkbox("Show basis", &show_basis);
			ImGui::SliderFloat("Basis thickness scale", &basis_thickness_scale, 0.0f, 100.0f);
			ImGui::SliderFloat("Basis length scale", &basis_length_scale, 0.0f, 100.0f);

			ImGui::SliderFloat("Bounding box width", &width, 1, 16);
			ImGui::SliderFloat("Bounding box height", &height, 1, 10);
			ImGui::SliderFloat("Bounding box depth", &depth, 1, 16);

			ImGui::SliderFloat("Damping factor", &damping_factor, 0.0f, 1.0f);
			ImGui::SliderFloat("Gravity strength", &gravity_strength, 0.0f, 100.0f);
			ImGui::SliderFloat("Smoothing radius", &smoothing_radius, grid_sphere_radius, 2.0f);
			//ImGui::SliderFloat("Time step", &time_step, 0.0f, 0.5f);
			ImGui::SliderFloat("Target_density", &target_density, 0.0f, 700.0f);
			ImGui::SliderFloat("Pressure multiplier", &pressure_multiplier, 0.0f, 1000.0f);
			ImGui::SliderFloat("Viscosity strength", &viscosity_strength, 0.0f, 0.1f);
			ImGui::SliderFloat("Near pressure multiplier", &near_pressure_multiplier, 0.0f, 10.0f);
			ImGui::SliderFloat("Ray march step size", &step_size, 0.01f, 0.5f);
			ImGui::SliderFloat("Ray march density multiplier", &density_multiplier, 0.001f, 0.1f);
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
