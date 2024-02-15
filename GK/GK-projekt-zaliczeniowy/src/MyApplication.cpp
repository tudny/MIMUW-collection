#include "MyApplication.hpp"

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include <iostream>
#include <vector>
#include "PerlinNoise.hpp"

#include "LSystems.h"
#include "asset.hpp"
#include "glError.hpp"

#define MAX_INSTANCES 255

struct VertexType {
  glm::vec3 color;
  glm::vec3 position;
};

struct LSPar {
  glm::vec4 position;
  glm::vec3 ambient;
  glm::vec3 direct;
  glm::vec3 attenuation;
};

// float bellCurve(float x, float sigma, float mu) {
//   return exp(-pow(x - mu, 2) / (2 * pow(sigma, 2))) / (sigma * sqrt(2 *
//   M_PI));
// }

float heightMap(const glm::vec2 position) {
  static const siv::PerlinNoise::seed_type seed = 123456u;
  static const siv::PerlinNoise perlin{seed};
  // float s = 1. / 4.;
  // return bellCurve(position.x, s, 0) * bellCurve(position.y, s, 0);

  // return sin(position.x) * (sin(position.y) / (2.2 + cos(10 * position.x) +
  //                                              cos(10 * position.y)));
  // return abs(sin(position.x) * sin(position.y)) /
  //        (0.1 + abs(position.x * position.y));
  return perlin.noise2D(position.x, position.y) * perlin.noise2D(position.x / 3., position.y / 3.) / 6.f * glm::length(position);
  // return 0.0;
}

VertexType getHeightMap(const glm::vec2 position) {
  VertexType v;
  static const siv::PerlinNoise::seed_type seed = 123456u;
  static const siv::PerlinNoise perlin{seed};
  float h = heightMap(position);

  v.position = glm::vec3(position, h);

  float c = sin(h * 5.f) * 0.1 + 0.4 + 0.05 * (sin((perlin.noise2D(position.x, position.y)) * 50.f));
  // 67, 92, 25
  v.color = glm::vec4(67.0 / 256.0, c, 25.0 / 256.0, 1.0);
  return v;
}

MyApplication::MyApplication()
    : Application(),
      vertexShader(SHADER_DIR "/shader.vert", GL_VERTEX_SHADER),
      fragmentShader(SHADER_DIR "/shader.frag", GL_FRAGMENT_SHADER),
      geometryShader(SHADER_DIR "/shader.geom", GL_GEOMETRY_SHADER),
      leafVertexShader(SHADER_DIR "/leaf.vert", GL_VERTEX_SHADER),
      shaderProgram({vertexShader, geometryShader, fragmentShader}),
      leafShaderProgram({leafVertexShader, geometryShader, fragmentShader}),
      branchShaderProgram({leafVertexShader, geometryShader, fragmentShader}) {
  glCheckError(__FILE__, __LINE__);

  // creation of the mesh ------------------------------------------------------
  std::vector<VertexType> vertices;
  std::vector<GLuint> index;

  for (int y = 0; y <= size; ++y) {
    for (int x = 0; x <= size; ++x) {
      float xx = (x - size / 2) * 0.1f;
      float yy = (y - size / 2) * 0.1f;
      vertices.push_back(getHeightMap({xx, yy}));
    }
  }

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      index.push_back((x + 0) + (size + 1) * (y + 0));
      index.push_back((x + 1) + (size + 1) * (y + 0));
      index.push_back((x + 1) + (size + 1) * (y + 1));

      index.push_back((x + 1) + (size + 1) * (y + 1));
      index.push_back((x + 0) + (size + 1) * (y + 1));
      index.push_back((x + 0) + (size + 1) * (y + 0));
    }
  }

  std::cout << "vertices=" << vertices.size() << std::endl;
  std::cout << "index=" << index.size() << std::endl;

  // auto leaf_appender = [&](float x, float y, float z) {
  //   auto first_tri = vertices.size();
  //   vertices.push_back({
  //       .color = glm::vec3(0.0, 1.0, 0.0),
  //       .position = glm::vec3(x, y, z),
  //   });
  //   vertices.push_back({
  //       .color = glm::vec3(0.0, 1.0, 0.0),
  //       .position = glm::vec3(x + 0.2, y, z + 0.2),
  //   });
  //   vertices.push_back({
  //       .color = glm::vec3(0.0, 1.0, 0.0),
  //       .position = glm::vec3(x - 0.2, y, z + 0.2),
  //   });
  //   index.push_back(first_tri + 0);
  //   index.push_back(first_tri + 1);
  //   index.push_back(first_tri + 2);
  // };
  //
  // leaf_appender(0.0, 0.0, 2.0);

  // creation of the vertex array buffer----------------------------------------

  // vbo
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(VertexType), vertices.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // ibo
  glGenBuffers(1, &ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, index.size() * sizeof(GLuint), index.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  // vao
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  // bind vbo
  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  // map vbo to shader attributes
  shaderProgram.setAttribute("in_Position", 3, sizeof(VertexType), offsetof(VertexType, position));
  shaderProgram.setAttribute("in_Colour", 3, sizeof(VertexType), offsetof(VertexType, color));

  // bind the ibo
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

  // vao end
  glBindVertexArray(0);

  // -- leafs -----------------------------------------------------------------

  std::vector<VertexType> leafs_vertices;
  std::vector<GLuint> leafs_index;

  const auto color = glm::vec3(0.0, 1.0, 0.0);
  const auto color_middle = glm::vec3(0.0, 0.3, 0.0);

  VertexType middle = {color_middle, {0.0, 0.0, 1. / 3.}};
  VertexType b = {color, {0.0, 0.0, 0.0}};
  VertexType b1 = {color, {sqrt(3.0) / 6.0, 0.0, 1. / 6.}};
  VertexType b2 = {color, {-sqrt(3.0) / 6.0, 0.0, 1. / 6.}};
  VertexType m1 = {color, {0.35, 0.0, 0.45}};
  VertexType m2 = {color, {-0.35, 0.0, 0.45}};
  VertexType t1 = {color, {0.25, 0.0, 0.8}};
  VertexType t2 = {color, {-0.25, 0.0, 0.8}};
  VertexType t = {color, {0.0, 0.0, 1.0}};

  leafs_vertices = {middle, b, b1, b2, m1, m2, t1, t2, t};
  leafs_index = {0, 1, 2, 4, 6, 8, 7, 5, 3, 1};

  glGenBuffers(1, &leaf_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, leaf_vbo);
  glBufferData(GL_ARRAY_BUFFER, leafs_vertices.size() * sizeof(VertexType), leafs_vertices.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &leaf_ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, leaf_ibo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, leafs_index.size() * sizeof(GLuint), leafs_index.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glGenVertexArrays(1, &leaf_vao);
  glBindVertexArray(leaf_vao);

  glBindBuffer(GL_ARRAY_BUFFER, leaf_vbo);

  leafShaderProgram.setAttribute("position", 3, sizeof(VertexType), offsetof(VertexType, position));
  leafShaderProgram.setAttribute("color", 3, sizeof(VertexType), offsetof(VertexType, color));

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, leaf_ibo);
  glBindVertexArray(0);

  // -- branch ----------------------------------------------------------------

  const auto branch_color = glm::vec3(128.0 / 256.0, 79.0 / 256.0, 42.0 / 256.0);
  const auto branch_diff_color = glm::vec3(92.0 / 256.0, 54.0 / 256.0, 25.0 / 256.0);

  const float delta_phi = 2.0 * M_PI / (float)branch_granularity;
  float phi = 0.;
  std::vector<VertexType> branch_vertices;
  std::vector<GLuint> branch_index;

  const float r = 0.5;
  for (int i = 0; i < branch_granularity; i++, phi += delta_phi) {
    branch_vertices.push_back({.color = branch_color, .position = glm::vec3(r * std::cos(phi), r * std::sin(phi), 0)});
  }
  const auto branch_scalar = 0.75;
  for (int i = 0; i < branch_granularity; i++, phi += delta_phi) {
    branch_vertices.push_back(
        {.color = branch_diff_color, .position = glm::vec3(r * std::cos(phi) * branch_scalar, r * std::sin(phi) * branch_scalar, 1.)});
  }
  branch_vertices.push_back({.color = branch_color, .position = glm::vec3(0, 0, 1.07)});

  for (int i = 0; i < branch_granularity; ++i) {
    branch_index.push_back(i);
    branch_index.push_back((i + 1) % branch_granularity);
    branch_index.push_back((i + 1) % branch_granularity + branch_granularity);
    branch_index.push_back(i + branch_granularity);
  }

  branch_index.push_back(2 * branch_granularity);
  for (int i = 0; i < branch_granularity; ++i) {
    branch_index.push_back(i + branch_granularity);
  }

  glGenBuffers(1, &branch_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, branch_vbo);
  glBufferData(GL_ARRAY_BUFFER, branch_vertices.size() * sizeof(VertexType), branch_vertices.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &branch_ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, branch_ibo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, branch_index.size() * sizeof(GLuint), branch_index.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glGenVertexArrays(1, &branch_vao);
  glBindVertexArray(branch_vao);

  glBindBuffer(GL_ARRAY_BUFFER, branch_vbo);

  branchShaderProgram.setAttribute("position", 3, sizeof(VertexType), offsetof(VertexType, position));
  branchShaderProgram.setAttribute("color", 3, sizeof(VertexType), offsetof(VertexType, color));

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, branch_ibo);
  glBindVertexArray(0);

  generateTree(leafs, branches);

  std::cout << "leafs=" << leafs.size() << std::endl;
  std::cout << "branches=" << branches.size() << std::endl;

  const auto num_trees = 128;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution dis(-1.0, 1.0);
  const float ssize = 20.f;

  for (int i = 0; i < num_trees; ++i) {
    float x = (float) dis(gen) * ssize;
    float y = (float) dis(gen) * ssize;
    float rot = (float) dis(gen) * 2 * M_PI;
    float sc = 1. + dis(gen) * 0.5;
    forest.push_back(glm::scale(glm::rotate(glm::translate(glm::mat4(1.0), glm::vec3(x, y, 0.0)), rot, glm::vec3(0, 0, 1)), glm::vec3(sc, sc, sc)));
  }
}

void MyApplication::loop() {
  // exit on window close button pressed
  if (glfwWindowShouldClose(getWindow())) {
    exit();
  }

  float t = getTime();
  // set matrix : projection + view
  glm::vec4 eyepos = glm::vec4(10.0 * sin(t / 10), 10.0 * cos(t / 10), 5.0, 1.0);
  projection =
      glm::perspective(static_cast<float>(2.0 * std::atan(getHeight() / 1920.f)), getWindowRatio(), 0.1f, 100.f);
  view = glm::lookAt(glm::vec3(eyepos), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0, 0, 1));
  model = glm::mat4(1.0);

  const glm::mat4 view_projection = projection * view;

  // clear
  glClear(GL_COLOR_BUFFER_BIT);
  glClearColor(135. / 256., 206. / 256., 235. / 256., 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  shaderProgram.use();

  // send uniforms
  shaderProgram.setUniform("pm", projection);
  shaderProgram.setUniform("vm", view);
  shaderProgram.setUniform("mm", model);
  shaderProgram.setUniform("eyepos", eyepos);
  shaderProgram.setUniform("vpm", view_projection);
  // shaderProgram.setUniform("nls", 1);
  // shaderProgram.setUniform("ls_mask", 0x1);
  const auto light_source = glm::vec4(cos(t) * 10.0, sin(t) * 10.0, 10.0, 0.0);
  LSPar ls[] = {{
      light_source,
      glm::vec3(0.2, 0.2, 0.25),
      glm::vec3(0.8, 0.8, 0.8),
      glm::vec3(1.0, 1.0, 1.0),
  }};
  shaderProgram.setUniform("ls_position", ls[0].position);
  shaderProgram.setUniform("ls_ambient", ls[0].ambient);
  shaderProgram.setUniform("ls_direct", ls[0].direct);
  shaderProgram.setUniform("ls_attenuation", ls[0].attenuation);
  glBindVertexArray(vao);

  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

  glCheckError(__FILE__, __LINE__);
  glDrawElements(GL_TRIANGLES,             // mode
                 size * size * 2 * 3 + 3,  // count
                 GL_UNSIGNED_INT,          // type
                 NULL                      // element array buffer offset
  );

  glBindVertexArray(0);

  shaderProgram.unuse();

  for (int k = 0; k < forest.size(); k++) {
    leafShaderProgram.use();

    leafShaderProgram.setUniform("vpm", view_projection);
    leafShaderProgram.setUniform("eyepos", eyepos);
    leafShaderProgram.setUniform("ls_position", ls[0].position);
    leafShaderProgram.setUniform("ls_ambient", ls[0].ambient);
    leafShaderProgram.setUniform("ls_direct", ls[0].direct);
    leafShaderProgram.setUniform("ls_attenuation", ls[0].attenuation);

    glBindVertexArray(leaf_vao);

    int idx = 0;
    for (auto leaf : leafs) {
      if (idx >= MAX_INSTANCES)
        break;
      leaf = forest[k] * leaf;
      leafShaderProgram.setUniform("modelMatrix[" + std::to_string(idx) + "]", leaf);
      idx++;
    }
    glBindBuffer(GL_ARRAY_BUFFER, leaf_vbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, leaf_ibo);

    glCheckError(__FILE__, __LINE__);

    glDrawElementsInstanced(GL_TRIANGLE_FAN, 10, GL_UNSIGNED_INT, NULL, leafs.size());
    glBindVertexArray(0);

    leafShaderProgram.unuse();
    branchShaderProgram.use();

    branchShaderProgram.setUniform("vpm", view_projection);
    branchShaderProgram.setUniform("eyepos", eyepos);
    branchShaderProgram.setUniform("ls_position", ls[0].position);
    branchShaderProgram.setUniform("ls_ambient", ls[0].ambient);
    branchShaderProgram.setUniform("ls_direct", ls[0].direct);
    branchShaderProgram.setUniform("ls_attenuation", ls[0].attenuation);

    glBindVertexArray(branch_vao);

    idx = 0;
    for (auto branch : branches) {
      if (idx >= MAX_INSTANCES)
        break;
      branch = forest[k] * branch;
      branchShaderProgram.setUniform("modelMatrix[" + std::to_string(idx) + "]", branch);
      idx++;
    }
    glBindBuffer(GL_ARRAY_BUFFER, branch_vbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, branch_ibo);

    glCheckError(__FILE__, __LINE__);

    for (int i = 0; i < branch_granularity; ++i) {
      glDrawElementsInstanced(GL_TRIANGLE_FAN, 4, GL_UNSIGNED_INT, (GLvoid*)((i * 4) * sizeof(GLuint)),
                              branches.size());
    }

    glDrawElementsInstanced(GL_TRIANGLE_FAN, branch_granularity + 1, GL_UNSIGNED_INT,
                            (GLvoid*)((4 * branch_granularity) * sizeof(GLuint)), branches.size());
    glBindVertexArray(0);

    branchShaderProgram.unuse();
  }

  for (int i = 0; i < leafs.size(); ++i) {
    auto leafs_movment = glm::rotate(glm::mat4(1.0f), 0.008f * (float)cos(t) * (float)sin(i), glm::vec3(1, 0, 0));
    leafs[i] = leafs[i] * leafs_movment;
  }
}
