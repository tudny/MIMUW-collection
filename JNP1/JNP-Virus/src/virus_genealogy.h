#ifndef VIRUS__VIRUS_GENEALOGY_H
#define VIRUS__VIRUS_GENEALOGY_H

#include <concepts>
#include <vector>
#include <memory>
#include <map>
#include <set>

class VirusAlreadyCreated : public std::exception {
public:
    const char *what() const noexcept override {
        return "VirusAlreadyCreated";
    }
};

class VirusNotFound : public std::exception {
public:
    const char *what() const noexcept override {
        return "VirusNotFound";
    }
};

class TriedToRemoveStemVirus : public std::exception {
public:
    const char *what() const noexcept override {
        return "TriedToRemoveStemVirus";
    }
};

template <typename Virus>
class VirusGenealogy {

    struct virus_order;

    using virus_id_t = typename Virus::id_type;
    using virus_ptr_t = std::shared_ptr<Virus>;
    using virus_set_t = std::set<virus_ptr_t, virus_order>;
    using virus_graph_t = std::map<virus_ptr_t, std::shared_ptr<virus_set_t>, virus_order>;
    using virus_id_to_ptr_t = std::map<virus_id_t, virus_ptr_t>;

    struct virus_order {
        bool operator()(const virus_ptr_t &a, const virus_ptr_t &b) const {
            return a->get_id() < b->get_id();
        }
    };

    /*
     * To jest klasa obsługująca odporność na wyjątki.
     * Używamy jej przy funkcjach connect, remove, create.
     */
    class VirusGuard {
        virus_id_to_ptr_t &mapper;
        virus_graph_t &graph, &transpose;

        size_t to_re_roll_mapper = 0,
                to_re_roll_graph = 0,
                to_re_roll_transpose = 0,
                to_re_roll_graph_edges = 0,
                to_re_roll_transpose_edges = 0;

        // vectory iteratorów elementów, które chcemy usunąć
        // (remove - w przypadku powodzenia, connect i create - w przypadku wyjątku)
        std::vector<typename virus_id_to_ptr_t::iterator> mapper_its;
        std::vector<typename virus_graph_t::iterator> graph_its, transpose_its;
        std::vector<std::pair<typename virus_graph_t::iterator, typename virus_set_t::iterator>> graph_edge_its,
                                                                                                 transpose_edge_its;
        // zmienna informująca czy będzie potrzebne usuwanie
        bool re_roll = true;

        template<typename T>
        T &last_of_vector(std::vector<T> &v) {
            return v[v.size() - 1];
        }

    public:
        VirusGuard() = delete;

        VirusGuard(virus_id_to_ptr_t &m_mapper, virus_graph_t &m_graph, virus_graph_t &m_transpose) noexcept
            : mapper(m_mapper), graph(m_graph), transpose(m_transpose) {}

        // destruktor klasy - jeśli trzeba, usuwa iteratory
        // erase(it) jest no-throw zatem tutaj nie możemy otrzymać wyjątku
        ~VirusGuard() noexcept {
            if (re_roll) {
                for (size_t i = 0; i < to_re_roll_graph_edges; ++i) {
                    const auto &[set, it] = graph_edge_its[i];
                    set->second->erase(it);
                }

                for (size_t i = 0; i < to_re_roll_transpose_edges; ++i) {
                    const auto &[set, it] = transpose_edge_its[i];
                    set->second->erase(it);
                }

                for (size_t i = 0; i < to_re_roll_mapper; ++i) {
                    const auto &it = mapper_its[i];
                    mapper.erase(it);
                }

                for (size_t i = 0; i < to_re_roll_graph; ++i) {
                    const auto &it = graph_its[i];
                    graph.erase(it);
                }

                for (size_t i = 0; i < to_re_roll_transpose; ++i) {
                    const auto &it = transpose_its[i];
                    transpose.erase(it);
                }
            }
        }

        void cancel_re_roll() noexcept {
            re_roll = false;
        }

        void set_erasing(bool rollback) noexcept {
            re_roll = rollback;
        }

        // bezpiecznie tworzymy miejsce na nowy węzeł
        // dzięki temu nie dostaniemy wyjątku po dodaniu wirusa do grafu
        void create_place_for_node() {
            mapper_its.emplace_back();
            graph_its.emplace_back();
            transpose_its.emplace_back();
        }

        // bezpiecznie tworzymy miejsce na nową krawędź
        void create_place_for_edge() {
            create_place_for_graph_edge();
            create_place_for_transpose_edge();
        }

        void create_place_for_transpose_edge() {
            transpose_edge_its.emplace_back();
        }

        void create_place_for_graph_edge() {
            graph_edge_its.emplace_back();
        }

        // funkcje dodające iteratory do wektorów
        // swap jest no-throw
        void add_mapper_it(typename virus_id_to_ptr_t::iterator it) {
            swap(last_of_vector(mapper_its), it);
            ++to_re_roll_mapper;
        }

        void add_graph_it(typename virus_graph_t::iterator it) {
            std::swap(last_of_vector(graph_its), it);
            ++to_re_roll_graph;
        }

        void add_transpose_it(typename virus_graph_t::iterator it) {
            std::swap(last_of_vector(transpose_its), it);
            ++to_re_roll_transpose;
        }

        void add_graph_edge_it(typename virus_graph_t::iterator g_it, typename virus_set_t::iterator s_it) {
            auto pair_to_add = std::make_pair(g_it, s_it);
            std::swap(last_of_vector(graph_edge_its), pair_to_add);
            ++to_re_roll_graph_edges;
        }

        void add_transpose_edge_it(typename virus_graph_t::iterator g_it, typename virus_set_t::iterator s_it) {
            auto pair_to_add = std::make_pair(g_it, s_it);
            std::swap(last_of_vector(transpose_edge_its), pair_to_add);
            ++to_re_roll_transpose_edges;
        }
    };

    typename virus_graph_t::iterator find_or_throw_not_found(virus_id_t const &id, virus_graph_t &graph) {
        auto virus_ptr_it = _id_mapper.find(id);
        if (virus_ptr_it == _id_mapper.end())
            throw VirusNotFound();

        return graph.find(virus_ptr_it->second);
    }

    static typename virus_graph_t::const_iterator find_or_throw_not_found_const(virus_id_t const &id,
                                                                                const virus_graph_t &graph,
                                                                                const virus_id_to_ptr_t &mapper) {
        const auto virus_ptr_it = mapper.find(id);
        if (virus_ptr_it == mapper.end())
            throw VirusNotFound();

        return graph.find(virus_ptr_it->second);
    }

    // funkcja tworząca pusty wierzchołek (bez krawędzi)
    // informacje o iteratorach, które chcemy dodać przechowujemy w guardzie
    // dzięki temu w razie otrzymania wyjątku możemy bezpiecznie je usunąć
    virus_ptr_t insert_empty_node(virus_id_t const virus_id, VirusGuard &guard) {
        guard.create_place_for_node();

        auto virus_ptr = std::make_shared<Virus>(virus_id);
        auto id_to_insert = std::make_pair(virus_id, virus_ptr);
        auto [map_it, was_map_inserted] = _id_mapper.insert(id_to_insert);
        guard.add_mapper_it(map_it);

        auto single_node_for_graph = std::make_pair(virus_ptr, std::make_shared<virus_set_t>());
        auto single_node_for_transpose = std::make_pair(virus_ptr, std::make_shared<virus_set_t>());

        auto [graph_it, was_graph_inserted] = _graph.insert(single_node_for_graph);
        guard.add_graph_it(graph_it);

        auto [transpose_it, was_transpose_inserted] = _transpose.insert(single_node_for_transpose);
        guard.add_transpose_it(transpose_it);

        return virus_ptr;
    }

    // klasa służąca do przechowywania rzeczywistej liczby rodziców wirusa
    // podczas wykonywania funkcji remove
    class ParentCounter {
        const virus_graph_t &transpose;
        const virus_id_to_ptr_t &mapper;
        std::map<virus_id_t, size_t> parent_counter;

        auto get_count_it(const virus_id_t &id) {
            auto count_it = parent_counter.find(id);
            if (count_it == parent_counter.end()) {
                auto transpose_it = find_or_throw_not_found_const(id, transpose, mapper);
                size_t how_many_parents = transpose_it->second->size();
                auto [new_it, added] = parent_counter.insert({ id, how_many_parents });
                count_it = new_it;
            }

            return count_it;
        }

    public:
        explicit ParentCounter(const virus_graph_t &m_transpose, const virus_id_to_ptr_t &m_mapper)
                                                                    : transpose(m_transpose), mapper(m_mapper) {}

        size_t get_parent_count(const virus_id_t &id) {
            return get_count_it(id)->second;
        }

        void remove_parent(const virus_id_t &id) {
            --get_count_it(id)->second;
        }
    };

    // wrzucamy do guarda iteratory, które chcemy usunąć
    void unsafe_remove_parents_connections(auto &transpose_list_it, auto &virus, auto &guard) {
        for (auto &parent : *transpose_list_it->second) {
            auto parent_it = _graph.find(parent);
            auto child_it = parent_it->second->find(virus);
            guard.create_place_for_graph_edge();
            guard.add_graph_edge_it(parent_it, child_it);
        }
    }

    // wrzucamy do guarda iteratory, które chcemy usunąć
    void unsafe_remove_children_connections(auto &graph_list_it, auto &virus, auto &guard, auto &counter) {
        for (auto &child : *graph_list_it->second) {
            auto child_it = _transpose.find(child);
            auto parent_it = child_it->second->find(virus);
            guard.create_place_for_transpose_edge();
            guard.add_transpose_edge_it(child_it, parent_it);
            counter.remove_parent(child_it->first->get_id());
        }
    }

    // w unsafe_remove wszsytko, co chcemy usunąć przechowujemy w obiekcie guard
    // jeśli nie otrzymamy wyjątku, wszytko co trzeba zostanie usunięte w destruktorze obiektu guard
    // w przeciwnym wypadku nic się nie stanie, graf zostanie nienaruszony
    void unsafe_remove(virus_ptr_t virus, VirusGuard &guard, ParentCounter &counter) {
        for (auto &child : *_graph.at(virus)) {
            if (counter.get_parent_count(child->get_id()) == 1) {
                unsafe_remove(child, guard, counter);
            }
        }

        auto mapper_it = _id_mapper.find(virus->get_id());
        auto graph_list_it = _graph.find(virus);
        auto transpose_list_it = _transpose.find(virus);

        unsafe_remove_parents_connections(transpose_list_it, virus, guard);
        unsafe_remove_children_connections(graph_list_it, virus, guard, counter);

        guard.create_place_for_node();
        guard.add_graph_it(graph_list_it);
        guard.add_transpose_it(transpose_list_it);
        guard.add_mapper_it(mapper_it);
    }

    // wrzucamy do guarda iteratory, które chcemy dodać
    // jeśli krawędź już istniała, nie dodajemy iteratora do guarda aby nie usunąć jej w razie wyjątku
    void unsafe_connect(virus_id_t const &child_id, virus_id_t const &parent_id, VirusGuard &guard) {

        guard.create_place_for_edge();

        auto child_it = find_or_throw_not_found(child_id, _transpose);
        auto parent_it = find_or_throw_not_found(parent_id, _graph);

        auto [graph_edge_it, was_graph_added] = parent_it->second->insert(child_it->first);
        if (was_graph_added) guard.add_graph_edge_it(parent_it, graph_edge_it);
        auto [transpose_edge_it, was_transpose_added] = child_it->second->insert(parent_it->first);
        if (was_transpose_added) guard.add_transpose_edge_it(child_it, transpose_edge_it);
    }

public:

    class children_iterator {
        typename virus_set_t::iterator it;

    public:

        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::iter_difference_t<typename virus_set_t::iterator>;
        using value_type = Virus;
        using pointer = Virus *;
        using reference = const Virus &;

        children_iterator() noexcept = default;

        explicit children_iterator(typename virus_set_t::iterator _it) noexcept : it(_it) {}

        children_iterator operator++(int) noexcept {
            children_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        children_iterator operator--(int) noexcept {
            children_iterator tmp = *this;
            --(*this);
            return tmp;
        }

        children_iterator &operator++() noexcept {
            it++;
            return *this;
        }

        children_iterator &operator--() noexcept {
            it--;
            return *this;
        }

        reference operator*() const {
            return **it;
        }

        pointer operator->() {
            return &(**it);
        }

        bool operator==(const children_iterator &other) const {
            return this->it == other.it;
        }

        bool operator!=(const children_iterator &other) const {
            return this->it != other.it;
        }
    };

    // Tworzy nową genealogię.
    // Tworzy także węzeł wirusa macierzystego o identyfikatorze stem_id.
    explicit VirusGenealogy(virus_id_t const &stem_id) {
        VirusGuard guard(_id_mapper, _graph, _transpose);
        _stem = insert_empty_node(stem_id, guard);
        guard.cancel_re_roll();
    }

    /* Według treści:
     * próba użycia konstruktora kopiującego lub operatora przypisania dla
     * obiektów klasy VirusGenealogy powinna zakończyć się błędem kompilacji;
     * */

    VirusGenealogy(const VirusGenealogy&) = delete;

    VirusGenealogy &operator=(const VirusGenealogy &) = delete;
    VirusGenealogy &operator=(VirusGenealogy &&) = delete;

    // Zwraca identyfikator wirusa macierzystego.
    virus_id_t get_stem_id() const {
        return _stem->get_id();
    }

    // Zwraca iterator pozwalający przeglądać listę identyfikatorów
    // bezpośrednich następników wirusa o podanym identyfikatorze.
    // Zgłasza wyjątek VirusNotFound, jeśli dany wirus nie istnieje.
    // Iterator musi spełniać koncept bidirectional_iterator oraz
    // typeid(*v.get_children_begin()) == typeid(const Virus &).
    VirusGenealogy<Virus>::children_iterator get_children_begin(virus_id_t const &id) const {
        auto it = find_or_throw_not_found_const(id, _graph, _id_mapper);
        return children_iterator(it->second->begin());
    }

    // Iterator wskazujący na element za końcem wyżej wspomnianej listy.
    // Zgłasza wyjątek VirusNotFound, jeśli dany wirus nie istnieje.
    VirusGenealogy<Virus>::children_iterator get_children_end(virus_id_t const &id) const {
        auto it = find_or_throw_not_found_const(id, _graph, _id_mapper);
        return children_iterator(it->second->end());
    }

    // Zwraca listę identyfikatorów bezpośrednich poprzedników wirusa
    // o podanym identyfikatorze.
    // Zgłasza wyjątek VirusNotFound, jeśli dany wirus nie istnieje.
    std::vector<virus_id_t> get_parents(virus_id_t const &id) const {
        auto it = find_or_throw_not_found_const(id, _transpose, _id_mapper);
        const auto &parents = *it->second;

        std::vector<virus_id_t> result;

        for (const auto &parent : parents) {
            result.push_back(parent->get_id());
        }

        return result;
    }

    // Sprawdza, czy wirus o podanym identyfikatorze istnieje.
    bool exists(virus_id_t const &id) const {
        return _id_mapper.contains(id);
    }

    // Zwraca referencję do obiektu reprezentującego wirus o podanym
    // identyfikatorze.
    // Zgłasza wyjątek VirusNotFound, jeśli żądany wirus nie istnieje.
    const Virus& operator[](virus_id_t const &id) const {
        auto it = find_or_throw_not_found_const(id, _graph, _id_mapper);
        return *it->first;
    }

    // Tworzy węzeł reprezentujący nowy wirus o identyfikatorze id
    // powstały z wirusów o podanym identyfikatorze parent_id lub
    // podanych identyfikatorach parent_ids.
    // Zgłasza wyjątek VirusAlreadyCreated, jeśli wirus o identyfikatorze
    // id już istnieje.
    // Zgłasza wyjątek VirusNotFound, jeśli któryś z wyspecyfikowanych
    // poprzedników nie istnieje.
    void create(virus_id_t const &id, virus_id_t const &parent_id) {
        std::vector<virus_id_t> singleton = { parent_id };
        create(id, singleton);
    }

    // tworzymy wirus o danym id i liście rodziców
    // wszystkie informacje, które chcemy dodać do grafu przechowujemy w guardzie
    // dzięki temu w razie otrzymania wyjątku zostaną one bezpiecznie usunięte
    void create(virus_id_t const &id, std::vector<virus_id_t> const &parent_ids) {
        if (parent_ids.empty())
            return;

        if (exists(id))
            throw VirusAlreadyCreated();

        for (auto &parent_id : parent_ids)
            find_or_throw_not_found(parent_id, _graph);

        VirusGuard guard(_id_mapper, _graph, _transpose);

        insert_empty_node(id, guard);

        for (auto parent_id : parent_ids)
            unsafe_connect(id, parent_id, guard);

        guard.cancel_re_roll();
    }

    // Dodaje nową krawędź w grafie genealogii.
    // Zgłasza wyjątek VirusNotFound, jeśli któryś z podanych wirusów nie istnieje.

    // łączymy wirusy (dziecko i rodzica) krawędzią
    // wszystkie informacje, które chcemy dodać do grafu przechowujemy w guardzie
    // dzięki temu w razie otrzymania wyjątku zostaną one bezpiecznie usunięte
    void connect(virus_id_t const &child_id, virus_id_t const &parent_id) {
        VirusGuard guard(_id_mapper, _graph, _transpose);

        unsafe_connect(child_id, parent_id, guard);
        guard.cancel_re_roll();
    }

    // Usuwa wirus o podanym identyfikatorze.
    // Zgłasza wyjątek VirusNotFound, jeśli żądany wirus nie istnieje.
    // Zgłasza wyjątek TriedToRemoveStemVirus przy próbie usunięcia
    // wirusa macierzystego.

    // usuwamy wirus o danym identyfikatorze
    // wszystkie informacje, które chcemy usunąć z grafu przechowujemy w guardzie
    // dzięki temu w razie otrzymania wyjątku graf pozostanie nienaruszony
    void remove(virus_id_t const &id) {
        auto virus_it = find_or_throw_not_found(id, _graph);

        if (virus_it->first->get_id() == get_stem_id())
            throw TriedToRemoveStemVirus();

        VirusGuard guard(_id_mapper, _graph, _transpose);
        guard.set_erasing(false);
        ParentCounter counter(_transpose, _id_mapper);

        unsafe_remove(virus_it->first, guard, counter);

        guard.set_erasing(true);
    }

private:
    virus_graph_t _graph, _transpose;
    virus_ptr_t _stem;
    virus_id_to_ptr_t _id_mapper;
};

#endif // VIRUS__VIRUS_GENEALOGY_H