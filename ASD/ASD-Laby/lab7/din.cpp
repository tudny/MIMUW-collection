#include <iostream>
#include <string>
#define assert(x) do { \
} while(false)

using namespace std;

template<typename T>
T maxx(T a) {
    return a;
}

template<typename T>
T maxx(T a, T b) {
    return a > b ? a : b;
}

template<typename T, typename... Args>
T maxx(T t, Args... args) {
    return maxx(t, maxx(args...));
}

namespace Splay {
    bool isDNA(char c) {
        switch (c) {
            case 'A': case 'T': case 'C': case 'G':
                return true;
            default:
                return false;
        }
    }

    enum CHILD {
        RIGHT = true,
        LEFT = false
    };

    struct Node {
        Node *child[2];
        Node *parent;
        char data, leftC, rightC;

        bool flipped;
        bool hasToBeUpdated;
        int maks, maksL, maksP, ile;
    };

    void updateInfo(Node *);

    Node *generateNode(char c) {
        assert(isDNA(c));
        Node *node = new Node();

        node->child[0] = node->child[1] = node->parent = nullptr;
        node->flipped = false;
        node->hasToBeUpdated = false;
        node->data = node->leftC = node->rightC = c;
        node->ile = node->maks = node->maksL = node->maksP = 1;

        return node;
    }

    void destroyNode(Node *node) {
        if (node == nullptr) return;
        destroyNode(node->child[0]);
        destroyNode(node->child[1]);
        free(node);
    }

    Node *getChild(Node *node, CHILD ch) {
        assert(node);
        updateInfo(node);
        return node->child[ch];
    }

    Node *getChildUNSAFE(Node *node, CHILD ch) {
        assert(node);
        return node->child[ch];
    }

    void setChild(Node *node, CHILD ch, Node *child) {
        assert(node);
        updateInfo(node);

        assert(node->child[ch] == nullptr);

        node->child[ch] = child;
        if (child != nullptr) {
            child->parent = node;
        }

        node->hasToBeUpdated = true;

        updateInfo(node);
    }

    Node *detach(Node *node, CHILD ch) {
        assert(node);
        Node *child = getChild(node, ch);
        if (child != nullptr) child->parent = nullptr;
        node->child[ch] = nullptr;
        node->hasToBeUpdated = true;
        updateInfo(node);
        return child;
    }

    void attach(Node *node, CHILD ch, Node *child) {
        assert(node);
        setChild(node, ch, child);
        updateInfo(node);
    }

    void flip(Node *node);

    void updateInfo(Node *node) {
        assert(node != nullptr);

        if (!node->hasToBeUpdated) return;
        node->hasToBeUpdated = false;

        if (node->flipped) {
            swap(node->child[0], node->child[1]);
            node->flipped = false;
            flip(node->child[0]);
            flip(node->child[1]);
        }

        Node *left = getChildUNSAFE(node, LEFT);
        Node *right = getChildUNSAFE(node, RIGHT);

        if (left != nullptr) left->hasToBeUpdated = true;
        if (right != nullptr) right->hasToBeUpdated = true;

        if (left == nullptr && right == nullptr) {
            node->leftC = node->rightC = node->data;
            node->ile = node->maks = node->maksL = node->maksP = 1;
        }
        else if (left != nullptr && right == nullptr) {
            node->maks = maxx(
                    left->maks,
                    (left->maksP + 1) * (node->data == left->rightC)
                    );
            node->maksP = maxx(
                    (left->maksP + 1) * (node->data == left->rightC),
                    1
                    );
            node->maksL = maxx(
                    (left->ile == left->maksL && left->rightC == node->data)
                        * (left->maksL + 1),
                    left->maksL
                    );
            node->ile = left->ile + 1;
            node->rightC = node->data;
            node->leftC = left->leftC;
        }
        else if (left == nullptr && right != nullptr) {
            node->maks = maxx(
                    right->maks,
                    (right->maksL + 1) * (node->data == right->leftC)
                    );
            node->maksP = maxx(
                    right->maksP,
                    (right->ile == right->maksP && right->leftC == node->data)
                        * (right->maksP + 1)
                    );
            node->maksL = maxx(
                    (right->maksL + 1) * (node->data == right->leftC),
                    1
                    );
            node->ile = right->ile + 1;
            node->rightC = right->rightC;
            node->leftC = node->data;
        }
        else {
            node->maks = maxx(
                    left->maks,
                    right->maks,
                    (node->data == left->rightC) * (left->maksP + 1),
                    (node->data == right->leftC) * (right->maksL + 1),
                    (node->data == left->rightC) * (node->data == right->leftC)
                        * (left->maksP + 1 + right->maksL)
                    );

            node->maksP = maxx(
                    right->maksP,
                    (right->maksP == right->ile && node->data == right->leftC)
                        * (right->maksP + 1),
                    (right->maksP == right->ile && node->data == right->leftC && node->data == left->rightC)
                        * (right->maksP + 1 + left->maksP)
                    );

            node->maksL = maxx(
                    left->maksL,
                    (left->maksL == left->ile && node->data == left->rightC)
                        * (left->maksL + 1),
                    (left->maksL == left->ile && node->data == left->rightC && node->data == right->leftC)
                        * (left->maksL + 1 + right->maksL)
                    );

            node->ile = left->ile + 1 + right->ile;
            node->leftC = left->leftC;
            node->rightC = right->rightC;
        }
    }

    void flip(Node *node) {
        if (node) {
            node->flipped ^= true;
            node->hasToBeUpdated = true;
            swap(node->maksL, node->maksP);
            swap(node->rightC, node->leftC);
        }
    }

    void printNode(Node *node) {
        if (!node) return;

        std::cout << "{ " << node->data
                  << ", maks=" << node->maks
                  << ", maksL=" << node->maksL
                  << ", maksP=" << node->maksP
                  << ", ile=" << node->ile
                  << ", leftC=" << node->leftC
                  << ", rightC=" << node->rightC
                  << " }\n" << flush;
    }

    bool hasParent(Node *node) {
        assert(node);
        return node->parent != nullptr;
    }

    CHILD getType(Node *node) {
        assert(node);
        assert(hasParent(node));

        Node *parent = node->parent;
        Node *right = getChild(parent, RIGHT);
        Node *left = getChild(parent, LEFT);

        if (right == node)
            return RIGHT;
        if (left == node)
            return LEFT;

        assert(node != left && node != right);
        return LEFT;
    }

    int getIle(Node *node) {
        if (node == nullptr) return 0;
        return node->ile;
    }

    int calculateKey(Node *node, int parentKey) {
        assert(node);
        updateInfo(node);
        if (parentKey == -1) {
            return getIle(getChild(node, LEFT)) + 1;
        }

        bool isLeft = getType(node) == LEFT;
        if (isLeft)
            return getIle(getChild(node, LEFT)) + 1;
        else
            return parentKey + getIle(getChild(node, LEFT)) + 1;
    }

    Node *rotateRight(Node *node) {
        assert(node);
        Node *c = detach(node, LEFT);
        Node *l2 = detach(c, RIGHT);
        attach(node, LEFT, l2);
        attach(c, RIGHT, node);
        return c;
    }

    Node *rotateLeft(Node *node) {
        assert(node);
        Node *c = detach(node, RIGHT);
        Node *r1 = detach(c, LEFT);
        attach(node, RIGHT, r1);
        attach(c, LEFT, node);
        return c;
    }

    void printSingle(Node *node, ostream &os) {
        if (node != nullptr) {
            printSingle(getChild(node, LEFT), os);
            os << node->data << " ";
            printSingle(getChild(node, RIGHT), os);
        }
    }

    void printSequence(Node *root) {
        printSingle(root, cout);
        cout << std::endl;
    }

    Node *splay(Node *root, int key, int parentKey) {
        if (root == nullptr) return root;

        int rootKey = calculateKey(root, parentKey);

        if (rootKey == key)
            return root;

        if (rootKey > key) {

            Node *left = getChild(root, LEFT);

            if (left == nullptr) return root;

            int leftKey = calculateKey(left, rootKey);

            if (leftKey > key) {
                Node *x = detach(left, LEFT);
                x = splay(x, key, -1);
                attach(left, LEFT, x);

                root = rotateRight(root);
            }
            else if (leftKey < key) {
                Node *x = detach(left, RIGHT);
                x = splay(x,  key - leftKey, -1);
                attach(left, RIGHT, x);

                if (getChild(left, RIGHT) != nullptr) {
                    Node *y = detach(root, LEFT);
                    y = rotateLeft(y);
                    attach(root, LEFT, y);
                }
            }

            if (getChild(root, LEFT) != nullptr) {
                root = rotateRight(root);
            }

            return root;
        }
        else {
            Node *right = getChild(root, RIGHT);

            if (right == nullptr) return root;

            int rightKey = calculateKey(right, rootKey);

            if (rightKey > key) {
                Node *x = detach(right, LEFT);
                x = splay(x, key - rootKey, -1);
                attach(right, LEFT, x);

                if (getChild(right, LEFT) != nullptr) {
                    Node *y = detach(root, RIGHT);
                    y = rotateRight(y);
                    attach(root, RIGHT, y);
                }
            }
            else if (rightKey < key) {
                Node *x = detach(right, RIGHT);
                x = splay(x, key - rightKey, -1);
                attach(right, RIGHT, x);

                root = rotateLeft(root);
            }

            if (getChild(root, RIGHT) != nullptr) {
                root = rotateLeft(root);
            }

            return root;
        }
    }

    Node *insert(Node *root, char data, int key) {
        if (root == nullptr) return generateNode(data);

        root = splay(root, key, -1);
        int rootKey = calculateKey(root, -1);

        if (rootKey == key) return root;

        Node *child = generateNode(data);

        if (rootKey > key) {
            Node *x = detach(root, LEFT);
            attach(child, RIGHT, root);
            attach(child, LEFT, x);
        } else {
            Node *x = detach(root, RIGHT);
            attach(child, RIGHT, x);
            attach(child, LEFT, root);
        }

        return child;
    }
}

using Splay::Node;
using Splay::CHILD::LEFT;
using Splay::CHILD::RIGHT;

struct Data {
    int N;
    int Q;
    string s;
};

void readInput(Data &input) {
    cin >> input.N;
    cin >> input.Q;
    cin >> input.s;
    input.s = '$' + input.s;
}

struct DNA {
    string s;
    Splay::Node *root = nullptr;

    DNA(const string &s) : s(s) {
        for (int i = 1; i < (int) s.length(); ++i) {
            root = Splay::insert(root, s[i], i);
        }
    }

    size_t getN() {
        return s.length() - 1;
    }

    void performO(int j, int k) {
        root = Splay::splay(root, k, -1);
        Node *r = Splay::detach(root, RIGHT);
        root = Splay::splay(root, j, -1);
        Node *l = Splay::detach(root, LEFT);

        Splay::flip(root);
        root = Splay::splay(root, 1, -1);

        Splay::attach(root, LEFT, l);
        root = Splay::splay(root, k, -1);
        Splay::attach(root, RIGHT, r);
    }

    void performP(int j, int k, int l) {
        root = Splay::splay(root, k, -1);
        Node *right = Splay::detach(root, RIGHT);
        root = Splay::splay(root, j, -1);
        Node *left = Splay::detach(root, LEFT);

        right = Splay::splay(right, 1, -1);

        if (right == nullptr) {
            right = left;
        } else {
            Splay::attach(right, LEFT, left);
        }

        Node *base = right;

        if (base != nullptr) {
            base = Splay::splay(base, l, -1);

            if ((int) getN() + j < l + k + 1) {
                root = Splay::splay(root, 1, -1);
                Splay::attach(root, LEFT, base);
            }
            else {
                Node *leftLeftover = Splay::detach(base, LEFT);

                root = Splay::splay(root, k - j + 1, -1);

                attach(base, LEFT, root);

                base = Splay::splay(base, 1, -1);

                attach(base, LEFT, leftLeftover);

                root = base;
            }
        }
    }

    int performN(int j, int k) {
        root = Splay::splay(root, k, -1);
        Node *r = Splay::detach(root, RIGHT);
        root = Splay::splay(root, j, -1);
        Node *l = Splay::detach(root, LEFT);

        int toReturn = root->maks;

        Splay::attach(root, LEFT, l);
        root = Splay::splay(root, k, -1);
        Splay::attach(root, RIGHT, r);

        return toReturn;
    }
};


int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    Data input;
    readInput(input);

    DNA dna = input.s;
    int Q = input.Q;

    int j, k, l;

    while (Q --> 0) {
        char c;
        cin >> c;

        switch (c) {
            case 'O':
                cin >> j >> k;
                dna.performO(j, k);
                break;
            case 'P':
                cin >> j >> k >> l;
                dna.performP(j, k, l);
                break;
            case 'N':
                cin >> j >> k;
                cout << dna.performN(j, k) << "\n";
                break;
        }
    }

    return 0;
}
