#include <errno.h>
#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>

#include "Tree.h"
#include "HashMap.h"
#include "path_utils.h"
#include "read_write_lock.h"

#define NEW_COMPONENT_ARRAY(X) \
char X[MAX_FOLDER_NAME_LENGTH + 1]; \
memset(X, 0, sizeof(X));

typedef struct Node_ {
    HashMap *dirs;
    char *name;
    rw_lock_t *lock;
    struct Node_ *parent;
} Node;

// Let's assume that we copy given name, as it can be from stack.
Node *node_new(const char *name) {
    Node *node = malloc(sizeof(Node));
    if (node == NULL) {
        return NULL;
    }

    memset(node, 0, sizeof(Node));

    node->dirs = hmap_new();
    if (node->dirs == NULL) {
        free(node);
        return NULL;
    }

    node->name = strdup(name);
    if (node->name == NULL) {
        hmap_free(node->dirs);
        free(node);
        return NULL;
    }

    node->lock = rw_lock_new();
    if (node->lock == NULL) {
        hmap_free(node->dirs);
        free(node->name);
        free(node);
        return NULL;
    }

    node->parent = NULL;

    return node;
}

void node_free(Node *node) {
    hmap_free(node->dirs);
    rw_lock_free(node->lock);
    free(node->name);
    free(node);
}

void node_free_recursive(Node *node) {
    const char *key;
    void *value;
    HashMapIterator it = hmap_iterator(node->dirs);
    while (hmap_next(node->dirs, &it, &key, &value)) {
        Node *child = value;
        node_free_recursive(child);
    }

    node_free(node);
}

bool has_child(Node *node, const char *name) {
    return hmap_get(node->dirs, name);
}

int add_child(Node *node, const char *name, Node **child) {
    assert(node);
    assert(name);

    if (has_child(node, name)) {
        return EEXIST;
    }

    *child = node_new(name);
    bool inserted = hmap_insert(node->dirs, name, *child);
    assert(inserted);
    (*child)->parent = node;

    return 0;
}

struct Tree {
    Node *root;
};

bool starts_with(const char *target, const char *source);

Tree *tree_new() {
    Node *root = node_new("");
    if (root == NULL)
        return NULL;

    Tree *tree = malloc(sizeof(Tree));
    if (!tree) {
        node_free(root);
        return NULL;
    }

    memset(tree, 0, sizeof(Tree));

    tree->root = root;

    return tree;
}

void tree_free(Tree *tree) {
    node_free_recursive(tree->root);
    free(tree);
}

// Assume the path is correct TODO delete
int path_length(const char *path) {
    int counter = 0;
    while (*path != '\0') {
        if (*path == '/')
            ++counter;
        ++path;
    }

    return counter;
}

typedef enum LockType_ {
    WRITE,
    READ
} LockType;

int unlock_all_node_at_path(Node *, LockType);

// Assume the path is correct
int get_node_with_lock(Tree *tree, const char *path, Node **end_node, LockType lockType) {
    int err;
    char component[MAX_FOLDER_NAME_LENGTH + 1];
    memset(component, 0, sizeof(component));

    Node *node = tree->root;
    Node *parent = NULL;
    if ((err = read_lock(node->lock)) != 0) {
        return err;
    }

    const char *subpath = path;
    while ((subpath = split_path(subpath, component))) {
        parent = node;
        if (!(node = hmap_get(node->dirs, component))) {
            unlock_all_node_at_path(parent, READ);
            return ENOENT;
        }

        if ((err = read_lock(node->lock)) != 0) {
            unlock_all_node_at_path(parent, READ);
            return err;
        }
    }

    // Na całą ścieżkę nakładamy READ locki.

    // Jeżeli lockType == WRITE to chcemy zmienić read lock na write lock
    if (lockType == WRITE) {
//        printf("Zmieniam typ locka.");
        if ((err = read_unlock(node->lock)) != 0) {
            unlock_all_node_at_path(parent, READ);
            return err;
        }
        if ((err = write_lock(node->lock)) != 0) {
            unlock_all_node_at_path(parent, READ);
            return err;
        }
    }

    *end_node = node;
    return 0;
}

int unlock_all_node_at_path(Node *node, LockType lockType) {
    bool first = true;
    int err;

    while (node) {
        if (first && lockType == WRITE) {
            first = false;
            if ((err = writer_unlock(node->lock)) != 0) {
                return err;
            }
        } else {
            if ((err = read_unlock(node->lock)) != 0) {
                return err;
            }
        }

        node = node->parent;
    }

    return 0;
}

int tree_create(Tree *tree, const char *full_path) {
    // Ścieżka jest niepoprawna.
    if (!is_path_valid(full_path))
        return EINVAL;

    // Nie możemy stworzyć roota.
    if (strcmp(full_path, "/") == 0)
        return EEXIST;

    NEW_COMPONENT_ARRAY(component_to_be_added);
    char *path = make_path_to_parent(full_path, component_to_be_added);

    Node *parent_node;
    int err;
    if ((err = get_node_with_lock(tree, path, &parent_node, WRITE)) != 0) {
        free(path); // Nawet jeżeli nie jest to ENOENT to robimy free, bo czemu nie.
        return err;
    }

    // Tutaj mamy już zapewniony write_lock na wierzchołku, do którego mamy dodać nowy.
    Node *child;
    if ((err = add_child(parent_node, component_to_be_added, &child)) != 0) {
        free(path);
        unlock_all_node_at_path(parent_node, WRITE);
        return err;
    }

    free(path);
    unlock_all_node_at_path(parent_node, WRITE);
    return 0;
}

char *tree_list(Tree *tree, const char *path) {
    if (!is_path_valid(path))
        return NULL;

    NEW_COMPONENT_ARRAY(component);

    Node *node;
    if (get_node_with_lock(tree, path, &node, READ) != 0) {
        return NULL;
    }

    char *list = make_map_contents_string(node->dirs);
    if (unlock_all_node_at_path(node, READ) != 0) {
        return NULL;
    }

    return list;
}

int tree_remove(Tree *tree, const char *full_path) {
    // Ścieżka jest niepoprawna.
    if (!is_path_valid(full_path))
        return EINVAL;

    // Nie możemy usunąć roota.
    if (strcmp(full_path, "/") == 0) {
        return EBUSY;
    }

    NEW_COMPONENT_ARRAY(deleted_component);
    char *path = make_path_to_parent(full_path, deleted_component);

    Node *parent_node;
    int err;
    if ((err = get_node_with_lock(tree, path, &parent_node, WRITE)) != 0) {
        free(path);
        return err;
    }

    free(path);

    // Tutaj mamy już WRITER_LOCK na parent'cie wierzchołka do usunięcia.
    Node *node;
    if ((node = hmap_get(parent_node->dirs, deleted_component)) == NULL) {
        unlock_all_node_at_path(parent_node, WRITE);
        return ENOENT;
    }

    if (hmap_size(node->dirs) > 0) {
        unlock_all_node_at_path(parent_node, WRITE);
        return ENOTEMPTY;
    }

    hmap_remove(parent_node->dirs, node->name);
    node_free(node);

    unlock_all_node_at_path(parent_node, WRITE);
    return 0;
}

int tree_move(Tree *tree, const char *source, const char *target) {
    const char *original_source = source,
               *original_target = target;

    if (!is_path_valid(source) || !is_path_valid(target))
        return EINVAL;

    if (strcmp(source, "/") == 0) {
        return EBUSY;
    }

    if (strcmp(target, "/") == 0) {
        return EEXIST;
    }

    if (starts_with(target, source)) {
        return -1; // TODO cannot move from source to source
    }

    Node *lca;
    int err;
    const char *lca_path = path_intersection(&source, &target);
    if ((err = get_node_with_lock(tree, lca_path, &lca, WRITE)) != 0) {
        free((char *) lca_path);
        return err;
    }

    free((char *) lca_path);

    NEW_COMPONENT_ARRAY(component);

    Node *source_node = lca;
    const char *source_subpath = source;
    while ((source_subpath = split_path(source_subpath, component))) {
        if (!(source_node = hmap_get(source_node->dirs, component))) {
            unlock_all_node_at_path(lca, WRITE);
            return ENOENT;
        }
    }

    NEW_COMPONENT_ARRAY(component_to_be_replaced);
    target = make_path_to_parent(target, component_to_be_replaced);

    if (target == NULL) {
        unlock_all_node_at_path(lca, WRITE);
        if (strcmp(original_source, original_target) == 0) {
            return 0;
        }
        return EEXIST;
    }

    Node *target_parent_node = lca;
    const char *target_parent_subpath = target;
    while ((target_parent_subpath = split_path(target_parent_subpath, component))) {
        if (!(target_parent_node = hmap_get(target_parent_node->dirs, component))) {
            free((char *) target);
            unlock_all_node_at_path(lca, WRITE);
            return ENOENT;
        }
    }

    free((char *) target);

    if (hmap_get(target_parent_node->dirs, component_to_be_replaced)) {
        unlock_all_node_at_path(lca, WRITE);
        return EEXIST;
    }

    hmap_remove(source_node->parent->dirs, source_node->name);

    free(source_node->name);
    source_node->name = strdup(component_to_be_replaced);

    hmap_insert(target_parent_node->dirs, source_node->name, source_node);
    source_node->parent = target_parent_node;

    unlock_all_node_at_path(lca, WRITE);
    return 0;
}

bool starts_with(const char *target, const char *source) {
    if (strlen(target) <= strlen(source))
        return false;

    while (*target && *source) {
        if (*target != *source)
            return false;

        ++target;
        ++source;
    }

    return true;
}

// TODO delete
void test() {

    const char *a = "/abc/";
    const char *b = "/abc/def/ghj/ijk/";
    const char *inter = path_intersection(&a, &b);
    printf("%s\n", inter);
    printf("%s\n", a);
    printf("%s\n", b);
    free((char *) inter);

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

//    path = "/a/b/c/";
//    code = tree_remove(tree, path);



    path = "/a/b/";
    list = tree_list(tree, path);

    printf("List %s\n", list);

    free((void *) list);

    // code = tree_move(tree, "/a/", "/z/x/");

    Node *node;
    get_node_with_lock(tree, "/a/b/c/", &node, WRITE);
    unlock_all_node_at_path(node, WRITE);

    tree_free(tree);
}
