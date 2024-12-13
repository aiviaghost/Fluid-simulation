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

	int num_particles = 3248;
	int sqrtN = sqrt(num_particles);
	float time_step = 1 / 120.0;
	float damping_factor = 0.95;
	float width = 16.0f, height = 9.0f;
	float half_width = width / 2, half_height = height / 2;
	float gravity_strength = 0.0;
	float smoothing_radius = 0.3;//1.2;
	float target_density = 2.75;
	float pressure_multiplier = 40.0;

	int PRIME1 = 86183;
	int PRIME2 = 7475723;

	volatile float avg_density = 0;

	std::vector<Node> nodes;
	std::vector<glm::vec3> positions, velocities;

	

	//float spacing = 3 * grid_sphere_radius;
	glm::vec3 square_center = glm::vec3(width / 10, height / 10, 0.0);
	float spacing = half_height / ((sqrtN - 1));
	for (int i = 0; i < sqrtN; i++) {
		for (int j = 0; j < sqrtN; j++) {
			auto pos = glm::vec3(-half_height / 2 + i * spacing, -half_height / 2 + j * spacing, 100.0f) + square_center;
			//auto pos = glm::vec3((rand() * 1.0f / RAND_MAX) * width - half_width, (rand() * 1.0f / RAND_MAX) * height - half_height, 100.0f);
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


	num_particles = positions.size();
	std::vector<float> densities(num_particles);
	std::vector<glm::vec3> predicted_positions(num_particles);

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

	std::vector<std::pair<int, int>> spatial(num_particles);
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

	auto update_spatial = [&](std::vector<glm::vec3> &vec) -> void {
		for (int i = 0; i < num_particles; i++) {
			spatial[i] = {point_to_hash(vec[i]), i};
		}

		sort(spatial.begin(), spatial.end());


		std::fill(start_inds.begin(), start_inds.end(), -1);
		start_inds[spatial[0].first] = 0;
		for (int i = 1; i < num_particles; i++) {
			if (spatial[i].first != spatial[i - 1].first) {
				start_inds[spatial[i].first] = i;
			}
		}
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


	const float mass = 1.0;
	auto calculate_density = [&](glm::vec3 point) -> float {
		float density = 0.0;

		int hash = point_to_hash(point);

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int cur_hash = hash + i * PRIME1 + j * PRIME2;
				cur_hash = good_mod(cur_hash, num_particles);
				int start_spatial_ind = start_inds[cur_hash];
				if (start_spatial_ind == -1) continue;
				for (int spatial_ind = start_spatial_ind; spatial_ind < num_particles && spatial[spatial_ind].first == cur_hash; spatial_ind++) {
					//float dist = glm::l2Norm(positions[spatial[spatial_ind].second] - point);
					float dist = glm::l2Norm(predicted_positions[spatial[spatial_ind].second] - point);
					float influence = smoothing_kernel(dist, smoothing_radius);
					density += mass * influence;
				}
			}
		}

		/*for (int i = 0; i < num_particles; i++) {
			float dist = glm::l2Norm(positions[i] - point);
			float influence = smoothing_kernel(dist, smoothing_radius);
			density += mass * influence;
		}*/

		return density;
	};


	auto update_densities = [&]() -> void {
		for (int i = 0; i < num_particles; i++) {
			densities[i] = calculate_density(positions[i]);
		}
	};

	auto density_to_pressure = [&](float density) -> float {
		float density_error = density - target_density;
		return density_error * pressure_multiplier;
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
		glm::vec3 point = predicted_positions[idx];

		int hash = point_to_hash(point);

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				int cur_hash = hash + i * PRIME1 + j * PRIME2;
				cur_hash = good_mod(cur_hash, num_particles);
				int start_spatial_ind = start_inds[cur_hash];
				if (start_spatial_ind == -1) continue;

				for (int spatial_ind = start_spatial_ind; spatial_ind < num_particles && spatial[spatial_ind].first == cur_hash; spatial_ind++) {
					int other_idx = spatial[spatial_ind].second;
					if (other_idx == idx) continue;

					//glm::vec3 offset = positions[other_idx] - point;
					glm::vec3 offset = predicted_positions[other_idx] - point;
					float dist = glm::l2Norm(offset);
					glm::vec3 dir = dist == 0 ? get_random_dir() : offset / dist;
					float slope = smoothing_kernel_derivative(dist, smoothing_radius);
					float density = densities[other_idx];
					float shared_pressure = calculcate_shared_pressure(density, densities[idx]);
					force += shared_pressure * dir * slope * mass / density;
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


	auto simulation_step = [&](float delta_time) -> void {
		for (int i = 0; i < num_particles; i++) {
			velocities[i] += glm::vec3(0.0, -1.0, 0.0) * gravity_strength * delta_time;
			predicted_positions[i] = positions[i] + velocities[i] * delta_time;
		}

		update_spatial(predicted_positions);

		for (int i = 0; i < num_particles; i++) {
			densities[i] = calculate_density(predicted_positions[i]);
		}

		for (int i = 0; i < num_particles; i++) {
			glm::vec3 pressure_force = calculate_pressure_force(i);
			glm::vec3 pressure_acceleration = pressure_force / densities[i];
			velocities[i] += pressure_acceleration * delta_time;
		}

		for (int i = 0; i < num_particles; i++) {
			positions[i] += velocities[i] * delta_time;
			handle_collision(i);
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

		//time_step = std::chrono::duration<float>(deltaTimeUs).count();
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
			ImGui::SliderFloat("Smoothing radius", &smoothing_radius, grid_sphere_radius, 100.0f);
			ImGui::SliderFloat("Time step", &time_step, 0.0f, 0.5f);
			ImGui::SliderFloat("Target_density", &target_density, 0.0f, 700.0f);
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
		edaf80::Fluid assignment5(framework.GetWindowManager());
		assignment5.run();
	}
	catch (std::runtime_error const& e) {
		LogError(e.what());
	}
}
