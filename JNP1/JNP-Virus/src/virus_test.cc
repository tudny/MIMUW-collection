#include <new>
#include <cstdio>
#include <cstdlib>

const bool print_debug = false;

const int THROW_FREQ = 50;
bool enable_throw = false;
int to_throw = THROW_FREQ;
int throw_cnt = 0;

bool should_fail() {
    if (enable_throw && --to_throw == 0) {
        to_throw = THROW_FREQ;
        throw_cnt++;
        return true;
    }
    return false;
}

void *operator new(std::size_t size) {
    if (should_fail()) {
        if (print_debug) {
            printf("bad alloc!\n");
            fflush(stdout);
        }
        throw std::bad_alloc();
    }
    if (size == 0)
        size = 1;
    void* p;
    while ((p = ::malloc(size)) == nullptr) {
        // If malloc fails and there is a new_handler,
        // call it to try free up memory.
        std::new_handler nh = std::get_new_handler();
        if (nh)
            nh();
        else
            throw std::bad_alloc();
    }
    return p;
}

#include <iostream>
#include <set>
#include <cassert>
#include <random>
#include "virus_genealogy.h"

using namespace std;

struct MyVeryOwnException {};

void bad() {
    if (should_fail()) {
        if (print_debug)
            cout << "exception!" << endl;
        throw MyVeryOwnException{};
    }
}

class IdType {
private:
    int v;

public:
    IdType() {
        bad();
    }

    IdType(int v) {
        bad();
        this->v = v;
    };

    IdType &operator=(IdType const &other) {
        bad();
        v = other.v;
        return *this;
    }

    int operator()() const {
        return v;
    }

    auto operator<=>(IdType const &other) const {
        bad();
        return v <=> other.v;
    }

    bool operator==(IdType const &other) const {
        bad();
        return v == other.v;
    }
};

class Virus {
public:
    using id_type = IdType;

    Virus(Virus::id_type _id) : id(_id) {
        bad();
    }

    Virus::id_type get_id() const {
        bad();
        return id;
    }

private:
    id_type id;
};

template<class Ex>
bool check_throws_up(auto f) {
    try {
        f();
    }
    catch(const Ex &e) {
        return true;
    }
    return false;
}

template<class T> bool has_same_content(std::vector<T> a, std::vector<T> b) {
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    return a == b;
}

using Iter = VirusGenealogy<Virus>::children_iterator;

bool has_same_addresses(Iter begin, Iter end, std::vector<const Virus*> addresses) {
    std::vector<const Virus*> v;
    for(; begin != end; ++begin)
        v.emplace_back(&(*begin));
    return has_same_content(v, addresses);
}

const int MAX_N = 50;
const int OPS = 50000;
int root = 0;
VirusGenealogy<Virus> G(root);
set<int> ch[MAX_N];
set<int> p[MAX_N];
mt19937 rng(10);

int rand(int a, int b) {
    return uniform_int_distribution<int>(a, b)(rng);
}

// jego rodzicom usunelismy juz dzieci
void remove_dfs(int w) {
    p[w].clear();
    for (int c : ch[w]) {
        if (p[c].size() == 1)
            remove_dfs(c);
        else
            p[c].erase(w);
    }
    ch[w].clear();
}

int main(int argc, char *argv[]) {
    check_throws_up<TriedToRemoveStemVirus>([&] {G.remove(root); });
    for (int op = 0; op < OPS; op++) {
        vector<IdType> free_nums;
        vector<IdType> occup_nums;
        for (int i = 0; i < MAX_N; i++) {
            if (i != 0 && p[i].empty())
                free_nums.push_back(i);
            else
                occup_nums.push_back(i);
        }

        int type = rand(0, 7);
        try {
            if (type == 0) { // ok create
                if (free_nums.empty())
                    continue;

                IdType id = free_nums[rand(0, free_nums.size() - 1)];
                while(occup_nums.back() > id)
                    occup_nums.pop_back();
                shuffle(occup_nums.begin(), occup_nums.end(), rng);

                int to_remove = rand(0, occup_nums.size());
                while (to_remove--)
                    occup_nums.pop_back();

                if (print_debug) {
                    cout << op << "/" << OPS << ": create(" << id() << ", {";
                    for (int i = 0 ; i < occup_nums.size() ; i++) {
                        cout << occup_nums[i]();
                        if (i + 1 < occup_nums.size())
                            cout << ", ";
                    }
                    cout << "});\n";
                }

                enable_throw = true;
                G.create(id, occup_nums);
                enable_throw = false;

                for (IdType w : occup_nums) {
                    ch[w()].insert(id());
                    p[id()].insert(w());
                }
            } else if (type == 1) { // ok remove
                if (occup_nums.size() == 1)
                    continue;

                IdType id = occup_nums[rand(1, occup_nums.size() - 1)];
                if (print_debug)
                    cout << op << "/" << OPS << ": remove(" << id() << ");\n";

                enable_throw = true;
                G.remove(id);
                enable_throw = false;

                for (int w : p[id()])
                    ch[w].erase(id());
                remove_dfs(id());
            } else if (type == 2) { // ok connect
                if (occup_nums.size() == 1)
                    continue;

                int a_ind = rand(0, occup_nums.size() - 2);
                IdType a = occup_nums[a_ind];
                IdType b = occup_nums[rand(a_ind + 1, occup_nums.size() - 1)];

                if (print_debug)
                    cout << op << "/" << OPS << ": connect(" << b() << ", " << a() << ");\n";
                enable_throw = true;
                G.connect(b, a);
                enable_throw = false;

                ch[a()].insert(b());
                p[b()].insert(a());
            } else if (type == 3) { // remove nonexistent
                if (free_nums.empty())
                    continue;

                IdType id = free_nums[rand(0, free_nums.size() - 1)];
                if (print_debug)
                    cout << op << "/" << OPS << ": remove(" << id() << ");\n";

                enable_throw = true;
                check_throws_up<VirusNotFound>([&] {G.remove(id);});
                enable_throw = false;
            } else if (type == 4) { // connect parent nonexistent
                if (free_nums.size() == 0)
                    continue;

                IdType c = occup_nums[rand(0, occup_nums.size() - 1)];
                IdType pp = free_nums[rand(0, free_nums.size() - 1)];
                if (print_debug)
                    cout << op << "/" << OPS << ": connect(" << c() << ", " << pp() << ");\n";

                enable_throw = true;
                check_throws_up<VirusNotFound>([&] {G.connect(c, pp);});
                enable_throw = false;
            } else if (type == 5) { // connect child nonexistent
                if (free_nums.size() == 0)
                    continue;

                IdType pp = occup_nums[rand(0, occup_nums.size() - 1)];
                IdType c = free_nums[rand(0, free_nums.size() - 1)];
                if (print_debug)
                    cout << op << "/" << OPS << ": connect(" << c() << ", " << pp() << ");\n";

                enable_throw = true;
                check_throws_up<VirusNotFound>([&] {G.connect(c, pp);});
                enable_throw = false;
            } else if (type == 6) { // ok create single parent
                if (free_nums.empty())
                    continue;

                IdType id = free_nums[rand(0, free_nums.size() - 1)];
                while(occup_nums.back() > id)
                    occup_nums.pop_back();
                shuffle(occup_nums.begin(), occup_nums.end(), rng);

                IdType par = occup_nums.back();

                if (print_debug)
                    cout << op << "/" << OPS << ": create(" << id() << ", " << par() << ");\n";
                enable_throw = true;
                G.create(id, par);
                enable_throw = false;

                ch[par()].insert(id());
                p[id()].insert(par());
            } else if (type == 7) { // create existent
                IdType id = occup_nums[rand(0, occup_nums.size() - 1)];
                IdType par = occup_nums[rand(0, occup_nums.size() - 1)];

                if (print_debug)
                    cout << op << "/" << OPS << ": create(" << id() << ", " << par() << ");\n";
                enable_throw = true;
                check_throws_up<VirusAlreadyCreated>([&] {G.create(id, par);});
                enable_throw = false;
            }
        } catch (MyVeryOwnException const &) {
            enable_throw = false;
            throw_cnt--;
        } catch (bad_alloc const &) {
            enable_throw = false;
            throw_cnt--;
        }

        assert(G.get_stem_id()() == root);
        for (int w = 0 ; w < MAX_N ; w++) {
            if (w != 0 && p[w].empty())
                assert(!G.exists(w));
            else {
                assert(G.exists(w));
                assert(has_same_content(G.get_parents(w),
                                        std::vector<IdType>(p[w].begin(), p[w].end())));
                vector<const Virus*> addresses;
                for (int u : ch[w])
                    addresses.push_back(&G[u]);
                assert(has_same_addresses(G.get_children_begin(w),
                                          G.get_children_end(w), addresses));
            }
        }
    }

    assert(throw_cnt == 0);

    return 0;
}
