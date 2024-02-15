#ifndef GK_APPLICATION_HPP
#define GK_APPLICATION_HPP

#include <string>

struct GLFWwindow;

class Application {
 public:
  Application();

  static Application& getInstance();

  [[nodiscard]] GLFWwindow* getWindow() const;

  void exit();

  [[nodiscard]] float getFrameDeltaTime() const;
  [[nodiscard]] float getTime() const;

  void run();

  int getWidth() const;
  int getHeight() const;
  float getWindowRatio() const;
  bool windowDimensionChanged() const;

 private:
  enum State { stateReady, stateRun, stateExit };

  State state;

  Application& operator=(const Application&) { return *this; }

  GLFWwindow* window;

  float time;
  float deltaTime;

  int width;
  int height;
  bool dimensionChanged;
  void detectWindowDimensionChange();

 protected:
  Application(const Application&) {}

  std::string title;

  virtual void loop();
};

#endif
