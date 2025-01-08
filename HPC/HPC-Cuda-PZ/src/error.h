#ifndef GPUGENV_ERROR_H
#define GPUGENV_ERROR_H

#include "logger.h"

#define PRINT_ERROR(msg) do { LOGERROR(msg); exit(1); } while (false)

void print_usage();

#endif //GPUGENV_ERROR_H
