#pragma once

#include <cassert>

#ifndef NDEBUG
#define DO_IN_DEBUG(code) do { code } while (0)
#else
#define DO_IN_DEBUG(code)
#endif
