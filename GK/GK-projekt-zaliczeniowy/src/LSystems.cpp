//
// Created by tudny on 1/17/24.
//

#include "LSystems.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include <map>
#include <string>
#include <vector>

const float phi = 15. / 360. * 2. * M_PI;

class LSystem {
 public:
  explicit LSystem(std::vector<glm::mat4>& leafs, std::vector<glm::mat4>& branches)
      : _leafs(leafs), _branches(branches) {}

  void make(glm::mat4 pos, float x, int n, bool is_first = false) {
    if (n == 0) {
      static int counter = 0;
      if (counter++ % 2 == 0) return;
      _leafs.push_back(glm::scale(pos, glm::vec3(x, x, x)));
      return;
    }

    const auto b_main = is_first ? x : x * 0.5f;
    auto branch_matrix = glm::scale(pos, glm::vec3(x * 0.1, x * 0.1, b_main));
    _branches.push_back(branch_matrix);

    const auto next = glm::translate(pos, glm::vec3(0, 0, b_main));

    const auto next1 = glm::rotate(next, 3 * phi, glm::vec3(0, 1, 0));
    const auto next11 = glm::rotate(next1, phi, glm::vec3(1, 0, 0));
    const auto next12 = glm::rotate(next1, -phi * 2.f, glm::vec3(1, 0, 0));
    // const auto next13 = glm::rotate(next1, phi * 4.f, glm::vec3(sqrt(0.5), sqrt(.5), 0));
    make(next11, x * 0.45f, n - 1);
    make(next12, x * 0.5f, n - 1);
    // make(next13, x * 0.6f, n - 1);

    const auto next2 = glm::rotate(next, -phi * 1.12312f, glm::vec3(0, 1, 0));

    const auto x_prim = x * 0.4f;
    auto branch_matrix2 = glm::scale(next2, 0.6f * glm::vec3(x * 0.1, x * 0.1, x_prim));
    _branches.push_back(branch_matrix2);

    const auto next_next = glm::translate(next2, 0.6f * glm::vec3(0, 0, x_prim));

    const auto next3 = glm::rotate(next_next, 2 * phi * 0.998998f, glm::vec3(sqrt(0.5), sqrt(0.5), 0));
    make(next3, x * 0.35f, n - 1);

    const auto next4 = glm::rotate(next_next, -2 * phi * .08136f, glm::vec3(sqrt(0.5), sqrt(0.5), 0));
    make(next4, x * 0.35f, n - 1);

    const auto next5 = glm::rotate(next_next, -6 * phi * 1.12436f, glm::vec3(sqrt(0.5), sqrt(0.5), 0));
    make(next5, x * 0.30f, n - 1);
  }

 private:
  std::vector<glm::mat4>& _leafs;
  std::vector<glm::mat4>& _branches;
};

void generateTree(std::vector<glm::mat4>& leafs, std::vector<glm::mat4>& branches) {
  LSystem l(leafs, branches);
  l.make(glm::mat4(1.0), 4, 4);
}
