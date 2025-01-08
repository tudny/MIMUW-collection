#pragma once

#include <filesystem>

namespace fs = std::filesystem;

bool file_exists(const fs::path &path);
