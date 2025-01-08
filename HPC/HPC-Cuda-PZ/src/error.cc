#include "error.h"

void print_usage() {
    std::cerr << (
            "Usage: gpugen [options]\n"
            "Options:\n"
            "  -i <database> <index>      Index a database\n"
            "  <query> <index> <output>   Query an index\n"
    );
    exit(1);
}
