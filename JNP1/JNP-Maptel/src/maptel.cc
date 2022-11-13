#include "maptel.h"
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <cassert>
#include <vector>
#include <sstream>

#define SAFE_ASSERT(pred) do {              \
    bool correct = (pred);                  \
    assert(correct && #pred);               \
    if (!(correct))                         \
        throw std::logic_error(#pred);      \
} while (false)

using phone_number_t = std::string;
using maptel_t = std::unordered_map<phone_number_t, phone_number_t>;
using map_id_t = unsigned long;
using maptels_set_t = std::unordered_map<map_id_t, maptel_t>;
using visited_phones_t = std::unordered_set<phone_number_t>;
using std::to_string;

using std::string;
using message_t = string;
using message_stream_t = std::stringstream;
using argument_list_t = std::vector<string>;

namespace {

#ifdef NDEBUG
    const bool debug = false;
#else
    const bool debug = true;
#endif

    bool is_char_end_of_string(const char &character) {
        return character == '\0';
    }

    bool check_tel_correct(const char *tel) {
        bool has_end_of_string;

        size_t i = 0;
        for (; i <= jnp1::TEL_NUM_MAX_LEN && !(has_end_of_string = is_char_end_of_string(tel[i])); ++i) {
            if (!isdigit(tel[i]))
                return false;
        }

        // Pusty numer nie jest poprawny, więc i musi być różne od 0.
        return has_end_of_string && i != 0;
    }

    template<typename K, typename V, typename I = typename std::unordered_map<K,V>::iterator>
    bool contains_iterator(const std::unordered_map<K, V> &m, const I &key_iterator) {
        return key_iterator != m.end();
    }

    template<typename K, typename V>
    bool contains(const std::unordered_map<K, V> &m, const K &key) {
        return m.find(key) != m.end();
    }

    template<typename T, typename I = typename std::unordered_set<T>::iterator>
    bool contains_iterator(const std::unordered_set<T> &m, const I &element_iterator) {
        return element_iterator != m.end();
    }

    template<typename T>
    bool contains(const std::unordered_set<T> &m, const T &element) {
        return m.find(element) != m.end();
    }

    template<typename T>
    string pointer_to_string(const T *pointer) {
        std::stringstream stream;
        stream << static_cast<const void *>(pointer);
        return stream.str();
    }

    maptels_set_t &get_maptels() {
        static maptels_set_t maptel;
        return maptel;
    }

    maptels_set_t::iterator get_maptels_iterator(const map_id_t &id) {
        auto maptel_iterator = get_maptels().find(id);

        SAFE_ASSERT(contains_iterator(get_maptels(), maptel_iterator));

        return maptel_iterator;
    }

    maptel_t &maptel_of_iterator(const maptels_set_t::iterator &maptels_iterator) {
        return maptels_iterator->second;
    }

    map_id_t next_id() {
        static map_id_t current_id = 0;
        return current_id++;
    }

    void print_debug_info(const message_t &message) {
        if (debug) {
            std::cerr << "maptel: " << message << "\n";
        }
    }

    void print_function_call(const string &function_name, const argument_list_t &arguments) {
        if (debug) {
            message_stream_t message_stream;

            message_stream << function_name << "(";
            for (size_t i = 0; i < arguments.size(); i++) {
                message_stream << arguments[i] << (i != arguments.size() - 1 ? ", " : "");
            }
            message_stream << ")";

            message_t message = message_stream.str();
            print_debug_info(message);
        }
    }

    void print_message_from_function(const string &function_name, const message_t &message) {
        if (debug) {
            message_stream_t new_message_stream;

            new_message_stream << function_name << ": " << message;

            message_t new_message = new_message_stream.str();
            print_debug_info(new_message);
        }
    }

    phone_number_t maptel_transform_unsafe(const maptel_t &phone_book, const phone_number_t &source,
                                           const string &function_name) {

        visited_phones_t visited_phones;
        phone_number_t pointer = source;

        while (contains(phone_book, pointer) && !contains(visited_phones, pointer)) {
            visited_phones.insert(pointer);
            pointer = phone_book.at(pointer);
        }

        if (contains(visited_phones, pointer)) {
            pointer = source;

            print_message_from_function(function_name, "cycle detected");
        }

        return pointer;
    }
}

map_id_t jnp1::maptel_create(void) {
    static const string function_name(__FUNCTION__);
    argument_list_t arguments_array;
    print_function_call(function_name, arguments_array);

    map_id_t id = next_id();

    SAFE_ASSERT(!contains(get_maptels(), id));

    get_maptels()[id];

    print_message_from_function(function_name, "new map id = " + to_string(id));

    return id;
}

void jnp1::maptel_delete(map_id_t id) {
    static const string function_name(__FUNCTION__);
    argument_list_t arguments_array{ to_string(id) };
    print_function_call(function_name, arguments_array);

    auto maptel_iterator = get_maptels_iterator(id);
    SAFE_ASSERT(contains_iterator(get_maptels(), maptel_iterator));

    get_maptels().erase(maptel_iterator);

    print_message_from_function(function_name, "map " + to_string(id) + " deleted");
}

void jnp1::maptel_insert(map_id_t id, const char *tel_src, const char *tel_dst) {
    static const string function_name(__FUNCTION__);
    argument_list_t arguments_array{ to_string(id), string(tel_src), string(tel_dst) };
    print_function_call(function_name, arguments_array);

    auto maptel_iterator = get_maptels_iterator(id);

    SAFE_ASSERT(tel_dst);
    SAFE_ASSERT(contains_iterator(get_maptels(), maptel_iterator));
    SAFE_ASSERT(check_tel_correct(tel_dst));
    SAFE_ASSERT(check_tel_correct(tel_src));

    phone_number_t src(tel_src);
    phone_number_t dest(tel_dst);

    maptel_of_iterator(maptel_iterator)[src] = dest;

    print_message_from_function(function_name, "inserted");
}

void jnp1::maptel_erase(map_id_t id, const char *tel_src) {
    static const string function_name(__FUNCTION__);
    argument_list_t arguments_array{ to_string(id), string(tel_src) };
    print_function_call(function_name, arguments_array);

    auto maptel_iterator = get_maptels_iterator(id);

    SAFE_ASSERT(tel_src);
    SAFE_ASSERT(contains_iterator(get_maptels(), maptel_iterator));
    SAFE_ASSERT(check_tel_correct(tel_src));

    phone_number_t src(tel_src);

    bool erased = maptel_of_iterator(maptel_iterator).erase(src);

    print_message_from_function(function_name, erased ? "erased" : "nothing to erase");
}

void jnp1::maptel_transform(map_id_t id, const char *tel_src, char *tel_dst, size_t len) {
    static const string function_name(__FUNCTION__);
    argument_list_t arguments_array{ to_string(id), string(tel_src), pointer_to_string(tel_dst), to_string(len) };
    print_function_call(function_name, arguments_array);

    const auto maptel_iterator = get_maptels_iterator(id);

    SAFE_ASSERT(tel_src && tel_dst);
    SAFE_ASSERT(contains_iterator(get_maptels(), maptel_iterator));
    SAFE_ASSERT(check_tel_correct(tel_src));

    const maptel_t &phone_book = maptel_of_iterator(maptel_iterator);
    const phone_number_t source(tel_src);
    const phone_number_t destination = maptel_transform_unsafe(phone_book, source, function_name);

    SAFE_ASSERT(destination.size() < len);

    auto number = destination.c_str();
    strncpy(tel_dst, number, len); // Ucinanie po osiągnięciu limitu (włącznie z \0)

    print_message_from_function(function_name, source + " -> " + number);
}


