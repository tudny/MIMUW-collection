#ifndef GK_MYAPPLICATION
#define GK_MYAPPLICATION

#include <vector>

#include "Application.hpp"
#include "Shader.hpp"

#define MAX_LEAFS 100

class MyApplication : public Application {
 public:
  MyApplication();

 protected:
  virtual void loop();

 private:
  const int size = 1000;

  // shader
  Shader vertexShader;
  Shader fragmentShader;
  Shader geometryShader;
  ShaderProgram shaderProgram;

  // shader matrix uniform
  glm::mat4 projection = glm::mat4(1.0);
  glm::mat4 view = glm::mat4(1.0);
  glm::mat4 model = glm::mat4(1.0);

  // leafs
  Shader leafVertexShader;
  ShaderProgram leafShaderProgram;
  std::vector<glm::mat4> leafs;
  std::vector<glm::mat4> branches;
  const int branch_granularity = 100;
  std::vector<glm::mat4> forest;

  ShaderProgram branchShaderProgram;

  // VBO/VAO/ibo
  GLuint vao, vbo, ibo;
  GLuint leaf_vao, leaf_vbo, leaf_ibo;
  GLuint branch_vao, branch_vbo, branch_ibo;
};

#endif
