//------- Ignore this ----------
#include<filesystem>
namespace fs = std::filesystem;
//------------------------------

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <random>
#include<iostream>
#include<glad/glad.h>
#include<GLFW/glfw3.h>
#include<stb/stb_image.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>

#include"Texture.h"
#include"shaderClass.h"
#include"VAO.h"
#include"VBO.h"
#include"EBO.h"
#include"Camera.h"
#include "FluidSimulation.h"

// Define precision for fluid simulation

const unsigned int width = 800;
const unsigned int height = 800;
float size = 1;

int s = (N + 2) * (N + 2) * 8;
int m = 0;
int t = 20 * 20 * 20 * 8;
GLfloat* vertices = new GLfloat[s];


// Vertices coordinates
/*GLfloat vertices[] =
{ //     COORDINATES     /        COLORS      /   TexCoord  //
		0.0f, 0.0f, 0.0f,     0.83f, 0.70f, 0.44f,	0.0f, 0.0f,
		0.1f, 0.1f, 0.1f,     0.83f, 0.70f, 0.44f,	5.0f, 0.0f,
		0.5f, 0.0f, -0.5f,     0.83f, 0.70f, 0.44f,	0.0f, 0.0f,
		 (size * 0.5f), (size * 0.0f),  (size * 0.5f),     0.83f, 0.70f, 0.44f,	5.0f, 0.0f,
		 (size * 0.0f), (size * 0.8f),  (size * 0.0f),     0.92f, 0.86f, 0.76f,	2.5f, 5.0f
};*/

// Indices for vertices order
GLuint indices[] =
{
	0, 1, 2,
	0, 2, 3,
	0, 1, 4,
	1, 2, 4,
	2, 3, 4,
	3, 0, 4
};

float generateRandomNumber() {
	// Create a random device and a generator
	std::random_device rd; // Get a random number from hardware
	std::mt19937 generator(rd()); // Seed the generator
	std::uniform_real_distribution<float> distribution(-1.0f, 1.0f); // Define the range

	return distribution(generator); // Generate the random number
}

/*void generateRandom(Shader& shaderProgram, int n)
{
	VAO vao;
	vao.Bind();

	GLfloat* va = new GLfloat[n];

	va[m++] = generateRandomNumber();
	va[m++] = generateRandomNumber();
	va[m++] = generateRandomNumber();
	va[m++] = 0.83f;
	va[m++] = 0.70f;
	va[m++] = 0.44f;
	va[m++] = 0.0f;
	va[m++] = 0.0f;

	VBO vbo(vertices, n * sizeof(GLfloat));
	vao.Bind(); // Bind the VAO
	vao.LinkAttrib(vbo, 0, 3, GL_FLOAT, sizeof(GLfloat) * 8, (void*)0);

	glPointSize(3);

	glDrawArrays(GL_POINTS, 0, n / 8);

	delete[] va;
}*/

void draw_dens(FluidSimulation& fluidSim, Shader& shaderProgram) 
{
	VAO vao;
	vao.Bind();

	int size = (N + 2) * (N + 2) * (N + 2);

	// iterate over all densities
	// assign a point to said density
	// render point on screen
	
	for (int i = 0; i < N + 2; i++) {
		for (int j = 0; j < N + 2; j++) {
			float density = fluidSim.dens[IX(i, j)];
				
			if (density > 0)
			{
				float x = (i / (float)(N + 2)) * 2.0f - 1.0f;
				float y = (j / (float)(N + 2)) * 2.0f - 1.0f; // Normalized y
					
				vertices[m++] = x;
				vertices[m++] = y;
				vertices[m++] = 0.0f;
				vertices[m++] = 0.83f;
				vertices[m++] = 0.70f;
				vertices[m++] = 0.44f;
				vertices[m++] = 0.0f;
				vertices[m++] = 0.0f;
			}
		}
	}

	GLint opacityLoc = glGetUniformLocation(shaderProgram.ID, "opacity");
	float newOpacity = 1.0f;
	glUniform1f(opacityLoc, newOpacity);
	// Create a VBO with the vertex positions
	VBO vbo(vertices, s * sizeof(GLfloat));
	vao.Bind(); // Bind the VAO
	vao.LinkAttrib(vbo, 0, 3, GL_FLOAT, sizeof(GLfloat) * 8, (void*)0);

	glPointSize(5);

	glDrawArrays(GL_POINTS, 0, s / 8);

	m = 0;

}

int main()
{
	srand(time(0)); // Seed random number generator
	// Initialize GLFW
	glfwInit();

	// Tell GLFW what version of OpenGL we are using 
	// In this case we are using OpenGL 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	// Tell GLFW we are using the CORE profile
	// So that means we only have the modern functions
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Create a GLFWwindow object of 800 by 800 pixels, naming it "YoutubeOpenGL"
	GLFWwindow* window = glfwCreateWindow(width, height, "YoutubeOpenGL", NULL, NULL);
	// Error check if the window fails to create
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	// Introduce the window into the current context
	glfwMakeContextCurrent(window);

	//Load GLAD so it configures OpenGL
	gladLoadGL();
	// Specify the viewport of OpenGL in the Window
	// In this case the viewport goes from x = 0, y = 0, to x = 800, y = 800
	glViewport(0, 0, width, height);


	// Generates Shader object using shaders default.vert and default.frag
	Shader shaderProgram("default.vert", "default.frag");


	// Generates Vertex Array Object and binds it
	VAO VAO1;
	VAO1.Bind();

	// Generates Vertex Buffer Object and links it to vertices
	VBO VBO1(vertices, sizeof(vertices));
	// Generates Element Buffer Object and links it to indices
	EBO EBO1(indices, sizeof(indices));

	// Links VBO attributes such as coordinates and colors to VAO
	VAO1.LinkAttrib(VBO1, 0, 3, GL_FLOAT, 8 * sizeof(float), (void*)0);
	VAO1.LinkAttrib(VBO1, 1, 3, GL_FLOAT, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	VAO1.LinkAttrib(VBO1, 2, 2, GL_FLOAT, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	// Unbind all to prevent accidentally modifying them
	VAO1.Unbind();
	VBO1.Unbind();
	EBO1.Unbind();



	/*
	* I'm doing this relative path thing in order to centralize all the resources into one folder and not
	* duplicate them between tutorial folders. You can just copy paste the resources from the 'Resources'
	* folder and then give a relative path from this folder to whatever resource you want to get to.
	* Also note that this requires C++17, so go to Project Properties, C/C++, Language, and select C++17
	*/
	std::string parentDir = (fs::current_path().fs::path::parent_path()).string();
	std::string texPath = "/Resources/YoutubeOpenGL 7 - Going 3D/";

	// Texture
	Texture brickTex((parentDir + texPath + "brick.png").c_str(), GL_TEXTURE_2D, GL_TEXTURE0, GL_RGBA, GL_UNSIGNED_BYTE);
	brickTex.texUnit(shaderProgram, "tex0", 0);

	// Original code from the tutorial
	/*Texture brickTex("brick.png", GL_TEXTURE_2D, GL_TEXTURE0, GL_RGBA, GL_UNSIGNED_BYTE);
	brickTex.texUnit(shaderProgram, "tex0", 0);*/



	// Enables the Depth Buffer
	glEnable(GL_DEPTH_TEST);

	// Creates camera object
	Camera camera(width, height, glm::vec3(0.0f, 0.0f, 2.0f));

	shaderProgram.Activate();

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 330");

	bool drawTriangle = false;

	float size = 1.0f;
	float color[4] = { 0.8f, 0.3f, 0.02f, 1.0f };
	glUseProgram(shaderProgram.ID);
	glUniform1f(glGetUniformLocation(shaderProgram.ID, "size"), size);
	glUniform4f(glGetUniformLocation(shaderProgram.ID, "color"), color[0], color[1], color[2], color[3]);
	float lastFrame = 0.0f;
	float deltaTime = 0.0f;
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Main while loop
	while (!glfwWindowShouldClose(window))
	{
		FluidSimulation fluidSim;
		float currentFrame = glfwGetTime();
		deltaTime = (currentFrame - lastFrame) * 50;
		float dt = currentFrame - lastFrame;
		// Start a new ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// Render the OpenGL content
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shaderProgram.Activate();

		ImGuiIO& io = ImGui::GetIO();
		if (!io.WantCaptureKeyboard && !io.WantCaptureMouse) {
			camera.Inputs(window, deltaTime);
		}

		camera.Matrix(45.0f, 0.1f, 100.0f, shaderProgram, "camMatrix");


		// Bind texture and VAO, then draw
		brickTex.Bind();
		VAO1.Bind();

		for (int i = 0; i < N + 2; i++) {
			for (int j = 0; j < N + 2; j++) {
				for (int k = 0; k < N + 2; k++) {
					// Randomize velocity directions in 3D
					fluidSim.vel_x[IX(i, j, k)] = ((rand() % 100) / 100.0f - 0.5f) * 0.1f;
					fluidSim.vel_y[IX(i, j, k)] = ((rand() % 100) / 100.0f - 0.5f) * 0.1f;

					// Randomize density
					fluidSim.dens[IX(i, j, k)] = (rand() % 100) / 100.0f * 0.5f;
				}
			}
		}
		
		//fluidSim.vel_step(N, fluidSim.vel_x, fluidSim.vel_y, fluidSim.vel_x_prev, fluidSim.vel_y_prev, 0.6f, 0.6f);
		//fluidSim.dens_step(N, fluidSim.dens, fluidSim.dens_prev, fluidSim.vel_x, fluidSim.vel_y, 0.5, 0.6f);
		draw_dens(fluidSim, shaderProgram);
		// generateRandom(shaderProgram, 100);

		if (drawTriangle) {
			glDrawElements(GL_TRIANGLES, sizeof(indices) / sizeof(int), GL_UNSIGNED_INT, 0);
		}

		ImGui::Begin("Test");
		ImGui::Text("hello world");
		ImGui::Checkbox("Draw Triangle", &drawTriangle);
		ImGui::SliderFloat("Size", &size, 0.5f, 2.0f);
		ImGui::ColorEdit4("Color", color);
		ImGui::End();

		glUseProgram(shaderProgram.ID);
		glUniform1f(glGetUniformLocation(shaderProgram.ID, "size"), size);
		glUniform4f(glGetUniformLocation(shaderProgram.ID, "color"), color[0], color[1], color[2], color[3]);
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		// Swap buffers and poll events
		glfwSwapBuffers(window);
		glfwPollEvents();
		lastFrame = currentFrame;
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();


	// Delete all the objects we've created
	VAO1.Delete();
	VBO1.Delete();
	EBO1.Delete();
	brickTex.Delete();
	shaderProgram.Delete();
	// Delete window before ending the program
	glfwDestroyWindow(window);
	// Terminate GLFW before ending the program
	glfwTerminate();
	return 0;
}