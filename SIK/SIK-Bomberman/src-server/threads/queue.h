#ifndef SIK_BOMBERMAN_QUEUE_H
#define SIK_BOMBERMAN_QUEUE_H

/* STD */
#include <shared_mutex>
#include <queue>
#include <functional>

namespace safe {

    template<typename T>
    class queue {
    public:
        queue() = default;

        void push(const T &element) {
            std::lock_guard write_lock(_mutex);
            _queue.push(element);
            inform_listener();
        }

        void push(T &&element) {
            std::lock_guard write_lock(_mutex);
            _queue.push(element);
            inform_listener();
        }

        T front() {
            std::lock_guard write_lock(_mutex);
            return _queue.front();
        }

        T take() {
            std::lock_guard write_lock(_mutex);
            if (_queue.empty()) {
                return T{};
            }
            T to_return = _queue.front();
            _queue.pop();
            inform_listener();
            return to_return;
        }

        void pop() {
            std::lock_guard write_lock(_mutex);
            _queue.pop();
            inform_listener();
        }

        [[nodiscard]] bool empty() {
            std::shared_lock read_lock(_mutex);
            return _queue.empty();
        }

        [[nodiscard]] size_t size() {
            std::shared_lock read_lock(_mutex);
            return _queue.size();
        }

        void clear() {
            std::lock_guard write_lock(_mutex);
            while (!_queue.empty()) _queue.pop();
            inform_listener();
        }

        void perform(std::function<void(std::queue<T> &)> func) {
            std::lock_guard write_lock(_mutex);
            func(_queue);
            inform_listener();
        }

        void lock() {
            _mutex.lock();
        }

        void unlock() {
            _mutex.unlock();
        }

        void set_listener_on_change(const callback_t &callback) {
            std::lock_guard write_lock(_mutex);
            listener = callback;
            has_listener = true;
        }

    private:
        void inform_listener() {
            if (has_listener)
                listener();
        }

        std::queue<T> _queue;
        std::shared_mutex _mutex;
        bool has_listener{false};
        callback_t listener;
    };
}


#endif //SIK_BOMBERMAN_QUEUE_H
