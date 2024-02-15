#include "Application.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <stdexcept>

using namespace std;

Application* currentApplication = nullptr;

Application& Application::getInstance() {
  if (currentApplication) {
    return *currentApplication;
  }
  throw std::runtime_error("There is no current Application");
}

Application::Application() : state(stateReady), width(640), height(480), title("Application") {
  currentApplication = this;

  cout << "[Info] GLFW initialisation" << endl;

  if (!glfwInit()) {
    throw std::runtime_error("Couldn't init GLFW");
  }

  int major = 3;
  int minor = 2;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, major);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, minor);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
  if (!window) {
    glfwTerminate();
    throw std::runtime_error("Couldn't create a window");
  }

  glfwMakeContextCurrent(window);

  glewExperimental = GL_TRUE;
  const GLenum err = glewInit();

  if (err != GLEW_OK) {
    glfwTerminate();
    throw std::runtime_error(string("Could initialize GLEW, error = ") + (const char*)glewGetErrorString(err));
  }

  // get version info
  const GLubyte* renderer = glGetString(GL_RENDERER);
  const GLubyte* version = glGetString(GL_VERSION);
  cout << "Renderer: " << renderer << endl;
  cout << "OpenGL version supported " << version << endl;

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
}

GLFWwindow* Application::getWindow() const {
  return window;
}

void Application::exit() {
  state = stateExit;
}

float Application::getFrameDeltaTime() const {
  return deltaTime;
}

float Application::getTime() const {
  return time;
}

void Application::run() {
  state = stateRun;

  glfwMakeContextCurrent(window);
  glfwWindowHint(GLFW_SAMPLES, 8);

  time = glfwGetTime();

  while (state == stateRun) {
    float t = glfwGetTime();
    deltaTime = t - time;
    time = t;

    detectWindowDimensionChange();

    loop();

    glfwSwapBuffers(window);

    glfwPollEvents();
  }

  glfwTerminate();
}

void Application::detectWindowDimensionChange() {
  int w, h;
  glfwGetWindowSize(getWindow(), &w, &h);
  dimensionChanged = (w != width || h != height);
  if (dimensionChanged) {
    width = w;
    height = h;
    glViewport(0, 0, width, height);
  }
}

void Application::loop() {
  cout << "[INFO] : loop" << endl;
}

int Application::getWidth() const {
  return width;
}

int Application::getHeight() const {
  return height;
}

float Application::getWindowRatio() const {
  return static_cast<float>(width) / static_cast<float>(height);
}

bool Application::windowDimensionChanged() const {
  return dimensionChanged;
}
