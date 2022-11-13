#ifndef SIK_BOMBERMAN_SET_H
#define SIK_BOMBERMAN_SET_H

#include "../utils/utils.h"

/* STD */
#include <set>
#include <shared_mutex>

namespace safe {

    template<typename T>
    class set {
        using insert_return_type = std::pair<typename std::set<T>::iterator, bool>;

    public:
        set() = default;

        insert_return_type add(T &&element) {
            std::lock_guard write_lock(_mutex);
            auto to_ret = _set.insert(std::move(element));
            on_change();
            return to_ret;
        }

        insert_return_type add(const T &element) {
            std::lock_guard write_lock(_mutex);
            auto to_ret = _set.insert(element);
            on_change();
            return to_ret;
        }

        bool remove(const T &element) {
            std::lock_guard write_lock(_mutex);
            auto rem = _set.erase(element);
            on_change();
            return rem;
        }

        [[nodiscard]] size_t size() {
            std::shared_lock read_lock(_mutex);
            return _set.size();
        }

        bool contains(const T &element) {
            std::shared_lock read_lock(_mutex);
            return _set.contains(element);
        }

        void for_each(std::function<void(T&)> applier) {
            std::lock_guard write_lock(_mutex);
            for (T element : _set) {
                applier(element);
            }
            on_change();
        }

        void for_each_unsafe(std::function<void(T&)> applier) {
            for (T element : _set) {
                applier(element);
            }
        }

        void perform(std::function<void(std::set<T>&)> applier) {
            applier(_set);
            on_change();
        }

        void set_listener(const callback_t &_listener) {
            std::lock_guard write_lock(_mutex);
            listener = _listener;
            has_listener = true;
        }

        void on_change() {
            if (has_listener) {
                listener();
            }
        }

        void lock() {
            _mutex.lock();
        }

        void unlock() {
            _mutex.unlock();
        }

    private:

        std::set<T> _set;
        std::shared_mutex _mutex;
        bool has_listener{false};
        callback_t listener;
    };
}


#endif //SIK_BOMBERMAN_SET_H
