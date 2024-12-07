#include "assignment5.hpp"
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

edaf80::Assignment5::Assignment5(WindowManager& windowManager) :
	mCamera(0.5f* glm::half_pi<float>(),
		static_cast<float>(config::resolution_x) / static_cast<float>(config::resolution_y),
		0.01f, 1000.0f),
	inputHandler(), mWindowManager(windowManager), window(nullptr)
{
	WindowManager::WindowDatum window_datum{ inputHandler, mCamera, config::resolution_x, config::resolution_y, 0, 0, 0, 0 };

	window = mWindowManager.CreateGLFWWindow("EDAF80: Assignment 5", window_datum, config::msaa_rate);
	if (window == nullptr) {
		throw std::runtime_error("Failed to get a window: aborting!");
	}

	bonobo::init();
}

edaf80::Assignment5::~Assignment5()
{
	bonobo::deinit();
}

void
edaf80::Assignment5::run()
{
	// Set up the camera
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, 6.0f));
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








	float grid_sphere_radius = 0.5;
	auto grid_sphere = parametric_shapes::createSphere(grid_sphere_radius, 8u, 8u);
	if (grid_sphere.vao == 0u) {
		LogError("Failed to retrieve the mesh for the grid_sphere");
		return;
	}

	float PI = 3.14159265358979;

	int num_particles = 200;
	float time_step = 1 / 120.0;
	float damping_factor = 1.0;
	int width = 80.0f, height = 45.0f;
	int half_width = width / 2, half_height = height / 2;
	float gravity_strength = 0.0;
	float smoothing_radius = 12.0;
	float target_density = 2.75;
	float pressure_multiplier = 0.8;


	std::vector<Node> nodes;
	std::vector<glm::vec3> positions, velocities, predicted_positions(num_particles);
	std::vector<float> densities(num_particles);

	float spacing = 3 * grid_sphere_radius;
	for (int i = 0; i * i <= num_particles; i++) {
		for (int j = 0; j * j <= num_particles; j++) {
			auto pos = glm::vec3(-half_width / 2 + i * spacing, -half_height / 2 + j * spacing, 100.0f);
			positions.push_back(pos);

			auto velocity = glm::vec3(0, 0, 0);
			velocities.push_back(velocity);

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

	auto handle_collision = [&](int idx) -> void {
			if (positions[idx].x + grid_sphere_radius > half_width) {
				positions[idx].x = half_width - grid_sphere_radius;
				velocities[idx].x *= -1 * damping_factor;
			}
			if (positions[idx].x - grid_sphere_radius < -half_width) {
				positions[idx].x = -half_width + grid_sphere_radius;
				velocities[idx].x *= -1 * damping_factor;
			}
			if (positions[idx].y + grid_sphere_radius > half_height) {
				positions[idx].y = half_height - grid_sphere_radius;
				velocities[idx].y *= -1 * damping_factor;
			}
			if (positions[idx].y - grid_sphere_radius < -half_height) {
				positions[idx].y = -half_height + grid_sphere_radius;
				velocities[idx].y *= -1 * damping_factor;
			}
		};

	auto smoothing_kernel = [&](float dist) -> float {
			if (dist >= smoothing_radius) {
				return 0;
			}
			float volume = (PI * std::pow(smoothing_radius, 4)) / 6;
			return (smoothing_radius - dist) * (smoothing_radius - dist) / volume;
		};

	auto smoothing_kernel_derivative = [&](float dist) -> float {
			if (dist >= smoothing_radius) {
				 return 0;
			}
			float scale = 12 / (std::pow(smoothing_radius, 4) * PI);
			return (dist - smoothing_radius) * scale;
		};

	auto get_density = [&](glm::vec3 point) -> float {
			float density = 0;
			for (int i = 0; i < num_particles; i++) {
				float dist = glm::length(positions[i] - point);
				density += smoothing_kernel(dist);
			}
			return density;
		};

	auto convert_density_to_pressure = [&](float density) -> float {
			float diff = density - target_density;
			return density * pressure_multiplier;
		};


	auto get_random_dir = []() -> glm::vec3 {
			return glm::vec3(rand() % 2 ? -1 : 1, rand() % 2 ? -1 : 1, 0);
		};

	auto get_pressure_force = [&](int idx) -> glm::vec3 {
			glm::vec3 pressure_force = glm::vec3(0.0f, 0.0f, 0.0f);
			float p1 = convert_density_to_pressure(densities[idx]);
			for (int i = 0; i < num_particles; i++) {
				if (i == idx) continue;
				glm::vec3 diff = positions[i] - positions[idx];
				float dist = glm::length(diff);
				glm::vec3 dir = dist > 0 ? diff / dist : get_random_dir();
				float slope = smoothing_kernel_derivative(dist);
				float p2 = convert_density_to_pressure(densities[i]);
				pressure_force += ((p1 + p2) / 2) * dir * slope / densities[i];
			}
			return pressure_force;
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

		auto & io = ImGui::GetIO();
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

		time_step = std::chrono::duration<float>(deltaTimeUs).count();
		if (!shader_reload_failed) {
			//
			// Todo: Render all your geometry here.
			//
			concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
				velocities[i] += glm::vec3(0.0, -gravity_strength, 0.0) * time_step;
				predicted_positions[i] = positions[i] + velocities[i] * (1 / 120.0f);
				
			});

			concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
				densities[i] = get_density(predicted_positions[i]);
				//std::cout << target_density << " " << densities[i] << std::endl;
			});

			concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
				glm::vec3 pressure_force = get_pressure_force(i);
				glm::vec3 pressure_acceleration = pressure_force / densities[i];
				velocities[i] += pressure_acceleration * time_step;
			});

			concurrency::parallel_for(size_t(0), size_t(num_particles), [&](size_t i) {
				positions[i] += velocities[i] * time_step;
				handle_collision(i);
			});

			for (int i = 0; i < num_particles; i++) {
				nodes[i].get_transform().SetTranslate(positions[i]);
				nodes[i].render(mCamera.GetWorldToClipMatrix());
			}
		}


		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//
		// Todo: If you want a custom ImGUI window, you can set it up
		//       here
		//
		bool const opened = ImGui::Begin("Scene Controls", nullptr, ImGuiWindowFlags_None);
		if (opened) {
			ImGui::Checkbox("Show basis", &show_basis);
			ImGui::SliderFloat("Basis thickness scale", &basis_thickness_scale, 0.0f, 100.0f);
			ImGui::SliderFloat("Basis length scale", &basis_length_scale, 0.0f, 100.0f);

			ImGui::SliderInt("Bounding box width", &width, 10, 1000);
			ImGui::SliderInt("Bounding box height", &height, 10, 1000);

			ImGui::SliderFloat("Damping factor", &damping_factor, 0.0f, 1.0f);
			ImGui::SliderFloat("Gravity strength", &gravity_strength, 0.0f, 100.0f);
			ImGui::SliderFloat("Smoothing radius", &smoothing_radius, 0.01f, 100.0f);
			//ImGui::SliderFloat("Time step", &time_step, 0.0f, 0.5f);
			ImGui::SliderFloat("Target_density", &target_density, 0.0f, 1.0f);
			ImGui::SliderFloat("Pressure multiplier", &pressure_multiplier, 0.0f, 1000.0f);
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
		edaf80::Assignment5 assignment5(framework.GetWindowManager());
		assignment5.run();
	}
	catch (std::runtime_error const& e) {
		LogError(e.what());
	}
}
