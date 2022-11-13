#include "virus_genealogy.h"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

class Virus {
public:
    using id_type = std::string;

    Virus(id_type const &_id) : id(_id) {
    }

    id_type get_id() const {
        return id;
    }

private:
    id_type id;
};

int main() {
    VirusGenealogy<Virus> gen("A1H1");
    Virus::id_type const id1 = gen.get_stem_id();
    assert(gen.exists(id1));
    assert(gen.get_parents(id1).size() == 0);
    gen.create("A", id1);
    gen.create("B", id1);
    size_t size;
    VirusGenealogy<Virus>::children_iterator ch_it;
    for (size = 0, ch_it = gen.get_children_begin(id1);
         ch_it != gen.get_children_end(id1);
         ++size, ++ch_it)
        assert((*ch_it).get_id() != id1);
    assert(size == 2);
    gen.create("C", "A");
    gen.connect("C", "B");
    assert(gen.get_parents("C").size() == 2);
    for (size = 0, ch_it = gen.get_children_begin("A");
         ch_it != gen.get_children_end("A");
         ++size, ++ch_it)
        assert(ch_it->get_id() == "C");
    assert(size == 1);
    std::vector<Virus::id_type> parents;
    parents.push_back("A");
    parents.push_back("B");
    gen.create("D", parents);
    assert(gen.get_parents("D").size() == parents.size());
    for (size = 0, ch_it = gen.get_children_begin("A");
         ch_it != gen.get_children_end("A");
         ++size, ++ch_it);
    assert(size == 2);
    assert("D" == gen["D"].get_id());
    gen.remove("A");
    assert(!gen.exists("A"));
    assert(gen.exists("B"));
    assert(gen.exists("C"));
    assert(gen.exists("D"));
    gen.remove("B");
    assert(!gen.exists("A"));
    assert(!gen.exists("B"));
    assert(!gen.exists("C"));
    assert(!gen.exists("D"));
    try {
        gen["E"];
    }
    catch (std::exception &e) {
        std::cout << e.what() << std::endl;
    }
    try {
        gen.create("E", "A1H1");
        gen.create("E", "A1H1");
    }
    catch (std::exception &e) {
        std::cout << e.what() << std::endl;
    }
    try {
        gen.remove("A1H1");
    }
    catch (std::exception &e) {
        std::cout << e.what() << std::endl;
    }
}
