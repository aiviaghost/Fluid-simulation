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

edaf80::Assignment5::Assignment5(WindowManager& windowManager) :
	mCamera(0.5f * glm::half_pi<float>(),
	        static_cast<float>(config::resolution_x) / static_cast<float>(config::resolution_y),
	        0.01f, 1000.0f),
	inputHandler(), mWindowManager(windowManager), window(nullptr)
{
	WindowManager::WindowDatum window_datum{ inputHandler, mCamera, config::resolution_x, config::resolution_y, 0, 0, 0, 0};

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

#define rep(i, a, b) for(int i = a; i < (b); i++)
#define sz(x) (int)(x.size())

typedef long long ll;
typedef std::vector<ll> vi;


typedef double T; // long double, Rational, double + mod<P>...
typedef std::vector<T> vd;
typedef std::vector<vd> vvd;

const T eps = 1e-8, inf = 1e20;
#define MP std::make_pair
#define ltj(X) if(s == -1 || MP(X[j],N[j]) < MP(X[s],N[s])) s=j

struct LPSolver {
	int m, n;
	vi N, B;
	vvd D;

	LPSolver(const vvd& A, const vd& b, const vd& c) :
		m(sz(b)), n(sz(c)), N(n + 1), B(m), D(m + 2, vd(n + 2)) {
		rep(i, 0, m) rep(j, 0, n) D[i][j] = A[i][j];
		rep(i, 0, m) { B[i] = n + i; D[i][n] = -1; D[i][n + 1] = b[i]; }
		rep(j, 0, n) { N[j] = j; D[m][j] = -c[j]; }
		N[n] = -1; D[m + 1][n] = 1;
	}

	void pivot(int r, int s) {
		T* a = D[r].data(), inv = 1 / a[s];
		rep(i, 0, m + 2) if (i != r && abs(D[i][s]) > eps) {
			T* b = D[i].data(), inv2 = b[s] * inv;
			rep(j, 0, n + 2) b[j] -= a[j] * inv2;
			b[s] = a[s] * inv2;
		}
		rep(j, 0, n + 2) if (j != s) D[r][j] *= inv;
		rep(i, 0, m + 2) if (i != r) D[i][s] *= -inv;
		D[r][s] = inv;
		std::swap(B[r], N[s]);
	}

	bool simplex(int phase) {
		int x = m + phase - 1;
		for (;;) {
			int s = -1;
			rep(j, 0, n + 1) if (N[j] != -phase) ltj(D[x]);
			if (D[x][s] >= -eps) return true;
			int r = -1;
			rep(i, 0, m) {
				if (D[i][s] <= eps) continue;
				if (r == -1 || MP(D[i][n + 1] / D[i][s], B[i])
					< MP(D[r][n + 1] / D[r][s], B[r])) r = i;
			}
			if (r == -1) return false;
			pivot(r, s);
		}
	}

	T solve(vd & x) {
		int r = 0;
		rep(i, 1, m) if (D[i][n + 1] < D[r][n + 1]) r = i;
		if (D[r][n + 1] < -eps) {
			pivot(r, n);
			if (!simplex(2) || D[m + 1][n + 1] < -eps) return -inf;
			rep(i, 0, m) if (B[i] == -1) {
				int s = 0;
				rep(j, 1, n + 1) ltj(D[i]);
				pivot(i, s);
			}
		}
		bool ok = simplex(1); x = vd(n);
		rep(i, 0, m) if (B[i] < n) x[B[i]] = D[i][n + 1];
		return ok ? D[m][n + 1] : inf;
	}
};

void
edaf80::Assignment5::run()
{
	// Set up the camera
	float pi = 3.1415;
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, -110.0f));
	mCamera.mMouseSensitivity = glm::vec2(0.003f);
	mCamera.mMovementSpeed = glm::vec3(3.0f); // 3 m/s => 10.8 km/h
	mCamera.mWorld.SetRotateY(pi);

	auto camera_position = mCamera.mWorld.GetTranslation();

	float elapsed_time_s = 0.0f;

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

	GLuint number_sphere_shader = 0u;
	program_manager.CreateAndRegisterProgram("Fallback",
		{ { ShaderType::vertex, "EDAF80/number_sphere.vert" },
		  { ShaderType::fragment, "EDAF80/number_sphere.frag" } },
		number_sphere_shader);
	if (number_sphere_shader == 0u) {
		LogError("Failed to load number_sphere shader");
		return;
	}

	GLuint eye_shader = 0u;
	program_manager.CreateAndRegisterProgram("Fallback",
		{ { ShaderType::vertex, "EDAF80/eye.vert" },
		  { ShaderType::fragment, "EDAF80/eye.frag" } },
		eye_shader);
	if (eye_shader == 0u) {
		LogError("Failed to load eye shader");
		return;
	}

	//
	// Todo: Insert the creation of other shader programs.
	//       (Check how it was done in assignment 3.)
	//

	//
	// Todo: Load your geometry
	//

	auto light_position = glm::vec3(-20.0f, 40.0f, -120.0f);

	auto const set_uniforms = [&light_position, &elapsed_time_s, &camera_position](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform1f(glGetUniformLocation(program, "elapsed_time_s"), elapsed_time_s);
		glUniform3fv(glGetUniformLocation(program, "camera_position"), 1, glm::value_ptr(camera_position));
	};

	float grid_sphere_radius = 5.0;
	auto grid_sphere = parametric_shapes::createSphere(grid_sphere_radius, 100u, 100u);
	if (grid_sphere.vao == 0u) {
		LogError("Failed to retrieve the mesh for the grid_sphere");
		return;
	}

	float eye_radius = 50.0;
	auto eye_sphere = parametric_shapes::createSphere(eye_radius, 100u, 100u);
	if (eye_sphere.vao == 0u) {
		LogError("Failed to retrieve the mesh for the eye_sphere");
		return;
	}


	auto set_number_sphere = [](int num, glm::vec3 color, Node& node) {
		bonobo::material_data material;
		material.diffuse = color;
		material.emissive.x = 1.0 * num;

		node.set_material_constants(material);
	};



	srand(1234);
	int rows = 5, cols = 5;
	std::vector<std::vector<int>> num_grid(rows, std::vector<int>(cols));
	std::vector<std::vector<Node>> grid(rows, std::vector<Node>(cols));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			grid[i][j].set_geometry(grid_sphere);
			grid[i][j].get_transform().SetTranslate(glm::vec3(j * 15, (rows - 1) * 15 - i * 15, -50));
			grid[i][j].set_program(&number_sphere_shader, set_uniforms);
			num_grid[i][j] = rand() % 10;
			set_number_sphere(num_grid[i][j], glm::vec3(1.0, 0.0, 1.0), grid[i][j]);
		}
	}

	glm::vec3 grid_middle = 0.5f * (grid[0][0].get_transform().GetTranslation() + grid[rows - 1][cols - 1].get_transform().GetTranslation());

	double minA = 101;
	int t = rows;
	int m = cols;
	vvd A(t, vd(m));
	rep(i, 0, t) rep(j, 0, m) {
		A[i][j] = num_grid[i][j];
		minA = std::min(minA, A[i][j]);
	}

	double alpha = -minA + 1;
	rep(i, 0, t) rep(j, 0, m) A[i][j] += alpha;

	// Solve v, Maj's strategy
	vd bMaj(t, 1);
	vd cMaj(m, 1);
	vd maj;
	T valMaj = LPSolver(A, bMaj, cMaj).solve(maj);
	double omegaModMaj = 1 / valMaj;
	rep(i, 0, m) maj[i] *= omegaModMaj;
	double omegaMaj = omegaModMaj - alpha;

	vd pref(m + 1);
	rep(i, 1, m + 1) {
		pref[i] = pref[i - 1] + maj[i - 1];
	}
	int n = 100000;
	vi maj_play(n);
	rep(i, 0, n) {
		double r = ((double)rand() / (double)RAND_MAX);

		int out = 1;
		while (out < m && r > pref[out]) out++;
		maj_play[i] = out - 1;
	}


	int cur_round = 0;

	Node eye;
	eye.set_geometry(eye_sphere);
	eye.set_program(&eye_shader);
	eye.get_transform().SetTranslate(glm::vec3(grid_middle.x, grid_middle.y, 50));

	auto project_to_line = [](glm::vec3 u, glm::vec3 v) -> glm::vec3 {
		return u - v * glm::dot(u, v);
	};

	int chosen_row = -1;
	int eye_chosen = -1;
	float chosen_time = -5;
	float eye_animation_started;
	glm::vec3 eye_rot_axis;
	float eye_rot_angle;

	bool points_added = false;
	int points = 0;

	float choice_duration = 3.0;

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
		auto const nowTime = std::chrono::high_resolution_clock::now();
		auto const deltaTimeUs = std::chrono::duration_cast<std::chrono::microseconds>(nowTime - lastTime);
		lastTime = nowTime;
		elapsed_time_s += deltaTimeUs.count() / 1e6;

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

		auto cam_center = mCamera.mWorld.GetTranslation();
		auto v = mCamera.mWorld.GetFront();

		if (elapsed_time_s - chosen_time > choice_duration) {
			chosen_row = -1;
			eye_chosen = -1;
		}

		int looking_row = -1;

		for (int i = 0; i < rows && chosen_row == -1; i++) {			for (int j = 0; j < cols && looking_row == -1; j++) {				glm::vec3 sphere_center = grid[i][j].get_transform().GetTranslation();
				if (glm::length(project_to_line(sphere_center - cam_center, v)) < grid_sphere_radius) {
					looking_row = i;
				}			}		}

		glm::vec3 sphere_standard_color = glm::vec3(1.0, 0.0, 1.0);
		glm::vec3 sphere_look_color = glm::vec3(0.0, 1.0, 0.0);
		glm::vec3 sphere_chosen_color = glm::vec3(1.0, 0.0, 0.0);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				glm::vec3 color;
				if (i == chosen_row && j == eye_chosen) {
					color = sphere_chosen_color;
				}
				else if (i == looking_row && chosen_row == -1) {
					color = sphere_look_color;
				}
				else {
					color = sphere_standard_color;
				}
				set_number_sphere(num_grid[i][j], color, grid[i][j]);			}
		}


		if ((inputHandler.GetKeycodeState(GLFW_KEY_SPACE) & JUST_PRESSED) && looking_row != -1 && chosen_row == -1) {
			chosen_time = elapsed_time_s;
			chosen_row = looking_row;
			eye_animation_started = elapsed_time_s;
			eye_chosen = maj_play[cur_round];
			cur_round++;


			glm::vec3 eye_should_look = grid[chosen_row][eye_chosen].get_transform().GetTranslation() - eye.get_transform().GetTranslation();
			eye_should_look = glm::normalize(eye_should_look);
			glm::vec3 current_eye_dir = glm::vec3(0.0, 0.0, -1.0);

			if (eye_chosen == cols / 2 && chosen_row == rows / 2 && cols % 2 == 1 && rows % 2 == 1) {
				eye_rot_angle = 2 * pi;
				eye_rot_axis = glm::vec3(0.0, 1.0, 0.0);
			}
			else {
				eye_rot_angle = glm::acos(glm::dot(current_eye_dir, eye_should_look));
				eye_rot_axis = glm::normalize(glm::cross(current_eye_dir, eye_should_look));
			}

			points_added = false;

			// eye.get_transform().SetRotate(rot_angle, rot_axis);
		}

		if (eye_chosen == -1) {			eye.get_transform().LookAt(glm::vec3(0.0, 0.0, -10000.0));		}


		if (eye_chosen != -1) {
			// Animation
			float from_angle;
			float to_angle;
			if (elapsed_time_s - eye_animation_started < choice_duration / 3.0) {				from_angle = 0.0;
				to_angle = eye_rot_angle;			}			else if (elapsed_time_s - eye_animation_started < 2.0 * choice_duration / 3.0) {				from_angle = eye_rot_angle;
				to_angle = eye_rot_angle;

				if (!points_added) {					points += num_grid[chosen_row][eye_chosen];
					points_added = true;				}			}			else {
				from_angle = eye_rot_angle;				to_angle = 0.0;			}

			float interp = glm::mod(elapsed_time_s - eye_animation_started, choice_duration / 3.0f) / (choice_duration / 3.0);
			float rot_angle = (1 - interp) * from_angle + interp * to_angle;

			eye.get_transform().SetRotate(rot_angle, eye_rot_axis);		}


		mWindowManager.NewImGuiFrame();

		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);


		if (!shader_reload_failed) {
			//
			// Todo: Render all your geometry here.
			//
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					grid[i][j].render(mCamera.GetWorldToClipMatrix());
				}
			}

			eye.render(mCamera.GetWorldToClipMatrix());
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
			ImGui::Text("Points: %i", points);
			ImGui::Text("If you played optimally: %f", omegaMaj * cur_round);
			ImGui::Text("Eye probabilities:");
			for (int i = 0; i < m; i++) ImGui::Text("%f", maj[m - i - 1]);
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
	} catch (std::runtime_error const& e) {
		LogError(e.what());
	}
}
