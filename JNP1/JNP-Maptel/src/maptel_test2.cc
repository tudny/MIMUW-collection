#include <cassert>
#include <cstddef>
#include <cstring>
#include "maptel.h"

namespace {
  unsigned long testuj() {
    unsigned long id;

    id = ::jnp1::maptel_create();
    ::jnp1::maptel_insert(id, "997", "112");

    return id;
  }

  unsigned long id = testuj();
} // anonymous namespace

int main() {
  char tel[::jnp1::TEL_NUM_MAX_LEN + 1];

  ::jnp1::maptel_transform(id, "997", tel, ::jnp1::TEL_NUM_MAX_LEN + 1);
  assert(::std::strcmp(tel, "112") == 0);
  ::jnp1::maptel_delete(id);
}
