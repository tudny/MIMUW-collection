cmake_minimum_required(VERSION 3.22)
project(grafika-komputerowa-projekt-zaliczeniowy)

add_executable(${CMAKE_PROJECT_NAME}
        src/Application.cpp
        src/Application.cpp
        src/MyApplication.cpp
        src/MyApplication.cpp
        src/glError.hpp
        src/glError.cpp
        src/main.cpp
        src/Shader.hpp
        src/Shader.cpp
        src/LSystems.cpp
        src/LSystems.h
)

set_property(TARGET ${CMAKE_PROJECT_NAME} PROPERTY CXX_STANDARD 20)
target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE -Wall -Wextra -O3)

add_definitions(-DGLEW_STATIC)
add_subdirectory(lib/glfw EXCLUDE_FROM_ALL)
add_subdirectory(lib/glew EXCLUDE_FROM_ALL)
add_subdirectory(lib/glm EXCLUDE_FROM_ALL)
include_directories(lib/PerlinNoise)

target_link_libraries(${CMAKE_PROJECT_NAME}
        PRIVATE glfw
        PRIVATE libglew_static
        PRIVATE glm
)

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/asset.hpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/src/asset.hpp
)
target_include_directories(${CMAKE_PROJECT_NAME}
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/src
)
