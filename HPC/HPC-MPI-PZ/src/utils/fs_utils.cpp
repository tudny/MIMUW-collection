#include "fs_utils.hpp"

bool file_exists(const fs::path &path) {
    return fs::exists(path) && fs::is_regular_file(path);
}
