
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#ifndef __DEBUG__
#define __DEBUG__ 0
#endif

#define DEBUG_PRINT(fmt, ...) \
    do { if (__DEBUG__) fprintf(stderr, "DEBUG: " fmt, ##__VA_ARGS__); } while (0)

#define BOOL_TRUE 1
#define BOOL_FALSE 0
#define BOOL_T int32_t

// LLVM equivalent:
// %struct.String = type { i32, i32, i8* }
typedef struct __attribute__((packed)) {
    int32_t len;
    BOOL_T is_heap_allocated;
    int8_t *buff;
} String;

void _exit_with_message(char *message) {
    DEBUG_PRINT("calling exit_with_message: %s\n", message);
    fprintf(stderr, "%s\n", message);
    exit(1);
}

void *_safe_malloc(size_t num, size_t size) {
    DEBUG_PRINT("calling safe_malloc: %lu, %lu\n", num, size);
    void *res = calloc(num, size);
    if (res == NULL) {
        _exit_with_message("ERROR WHILE ALLOCATING SPACE");
    }
    return res;
}

int32_t _get_str_len(String *s) {
    DEBUG_PRINT("calling get_str_len: %p\n", s);
    return s->len;
}

int8_t *_get_str_data(String *s) {
    DEBUG_PRINT("calling get_str_data: %p\n", s);
    return s->buff;
}

String *_new_str(int8_t *data, int32_t len, BOOL_T is_heap_allocated) {
    DEBUG_PRINT("calling new_str: %p, %d, %d\n", data, len, is_heap_allocated);
    String *res = _safe_malloc(sizeof(String), 1);
    res->len = len;
    res->buff = data;
    res->is_heap_allocated = is_heap_allocated;
    return res;
}

String *_new_str_from_literal(int8_t *data, int32_t len) {
    DEBUG_PRINT("calling new_str_from_literal: %p, %d\n", data, len);
    return _new_str(data, len, BOOL_FALSE);
}

void _remove_str(String *s) {
    DEBUG_PRINT("calling remove_str: %p\n", s);
    if (s->is_heap_allocated) {
        free(s->buff);
        s->buff = NULL;
    }
    free(s);
}

String *_concat_str(String *a, String *b) {
    DEBUG_PRINT("calling concat_str: %p, %p\n", a, b);
    int32_t len = a->len + b->len;
    int8_t *data = _safe_malloc(len, sizeof(int8_t));
    String *res = _new_str(data, len, BOOL_TRUE);

    int8_t *it = res->buff;
    for (int i = 0; i < a->len; i++) {
        *(it++) = a->buff[i];
    }
    for (int i = 0; i < b->len; i++) {
        *(it++) = b->buff[i];
    }
    return res;
}

BOOL_T _str_eq(String *a, String *b) {
    DEBUG_PRINT("calling str_eq: %p, %p\n", a, b);
    if (a->len != b->len) return BOOL_FALSE;
    int8_t *ai = a->buff;
    int8_t *bi = b->buff;

    for (int i = 0; i < a->len; ++i) {
        if (*(ai++) != *(bi++)) return BOOL_FALSE;
    }
    return BOOL_TRUE;
}

BOOL_T _str_neq(String *a, String *b) {
    DEBUG_PRINT("calling str_neq: %p, %p\n", a, b);
    return 1 - _str_eq(a, b);
}

void printString(String *s) {
    DEBUG_PRINT("calling print_string: %p\n", s);
    for (int i = 0; i < s->len; ++i) {
        putchar(s->buff[i]);
    }
    putchar('\n');
    fflush(stdout);
}

String *readString() {
    DEBUG_PRINT("calling read_string\n");
    int8_t *buffer = NULL;
    unsigned long len = 0;
    errno = 0;
    size_t read = getline((char **) &buffer, &len, stdin);
    if (errno == EINVAL) _exit_with_message("INVALID GETLINE ARGUMENTS");
    if (errno == ENOMEM) _exit_with_message("ERROR WHILE ALLOCATING SPACE");
    if (read > 0 && buffer[read - 1] == '\n') {
        buffer[read - 1] = 0;
        --read;
    }
    String *res = _new_str(buffer, read, BOOL_TRUE);
    res->is_heap_allocated = 1;
    return res;
}

void printInt(int32_t x) {
    DEBUG_PRINT("calling print_int: %d\n", x);
    printf("%d\n", x);
    fflush(stdout);
}

int32_t readInt() {
    DEBUG_PRINT("calling read_int\n");
    int x;
    scanf(" %d ", &x);
    return x;
}

void error() {
    DEBUG_PRINT("calling error\n");
    _exit_with_message("runtime error");
}

void *_memory_allocate(int32_t size) {
    DEBUG_PRINT("calling memory_allocate: sizeof=%d\n", size);
    return _safe_malloc(1, size);
}

void *_memory_allocate_many(int32_t count, int32_t size) {
    DEBUG_PRINT("calling memory_allocate_many: sizeof=%d, num=%d\n", size, count);
    return _safe_malloc(count, size);
}

void _memory_free(void *ptr) {
    DEBUG_PRINT("calling memory_free: ptr=%p\n", ptr);
    free(ptr);
}

void *_memory_allocate_many_strings(int32_t count, int8_t *empty_str) {
    String **res = _safe_malloc(count, sizeof(String *));
    // String are immutable, so we can reuse the same empty string
    String *empty = _new_str_from_literal(empty_str, 0);
    for (int i = 0; i < count; ++i) {
        res[i] = empty;
    }
    return res;
}
