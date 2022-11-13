#include <assert.h>
#include <string.h>
#include "maptel.h"

static const char t112[] = "112";
static const char t997[] = "997";
static const char t416[] = "225544416";
static const char t401[] = "225544401";
static const char t333[] = "223333333";

int main() {
  unsigned long id;
  char tel[TEL_NUM_MAX_LEN + 1]; /* +1 na terminalne zero */

  id = maptel_create();
  maptel_insert(id, t112, t997);
  maptel_insert(id, t416, t401);
  maptel_insert(id, t401, t333);
  maptel_transform(id, t112, tel, TEL_NUM_MAX_LEN + 1);
  assert(strcmp(tel, t997) == 0);
  maptel_transform(id, t333, tel, TEL_NUM_MAX_LEN + 1);
  assert(strcmp(tel, t333) == 0);
  maptel_erase(id, t333);
  maptel_erase(id, t112);
  maptel_transform(id, t112, tel, TEL_NUM_MAX_LEN + 1);
  assert(strcmp(tel, t112) == 0);
  maptel_transform(id, t416, tel, TEL_NUM_MAX_LEN + 1);
  assert(strcmp(tel, t333) == 0);
  maptel_insert(id, t333, t333);
  maptel_transform(id, t333, tel, TEL_NUM_MAX_LEN + 1);
  assert(strcmp(tel, t333) == 0);
  maptel_delete(id);
}
