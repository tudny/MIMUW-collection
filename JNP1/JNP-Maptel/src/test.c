#include "maptel.h"
#include <stdlib.h>
#include <string.h>

int main() {
    unsigned long id = maptel_create();

    char *str = "1234567890123456789012";
    maptel_insert(id, str, "123");

    char *dst = malloc(sizeof(char) * 3);
    maptel_transform(id, str, dst, 4);

    char *not_zero = malloc(sizeof(char) * 3);
    for (size_t i = 0; i < TEL_NUM_MAX_LEN; ++i) {
        not_zero[i] = '1';
    }

    not_zero[TEL_NUM_MAX_LEN] = '1';
    not_zero[TEL_NUM_MAX_LEN + 1] = '\0';
    maptel_insert(id, dst, not_zero);
}
