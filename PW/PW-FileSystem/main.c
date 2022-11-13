#include "HashMap.h"
#include "Tree.h"
#include "read_write_lock.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>


int test2(void) {
    Tree *tree = tree_new();

    int code = tree_create(tree, "/a/");



    code = tree_create(tree, "/a/b/");



    code = tree_create(tree, "/z/");



    code = tree_create(tree, "/a/b/c/");



    code = tree_create(tree, "/a/b/e/");



    code = tree_create(tree, "/a/b/e/");



    const char *path = "/a/b/";
    const char *list = tree_list(tree, path);



    free((void *) list);

    path = "/a/b/c/";
    code = tree_remove(tree, path);



    path = "/a/b/";
    list = tree_list(tree, path);



    free((void *) list);



    code = tree_move(tree, "/a/", "/z/x/");




    tree_free(tree);

    return 0;
}

int main(void) {
    test();

    return 0;
}
