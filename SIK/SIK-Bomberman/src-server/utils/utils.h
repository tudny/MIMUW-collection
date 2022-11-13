#ifndef SIK_BOMBERMAN_UTILS_H
#define SIK_BOMBERMAN_UTILS_H

/* STD */
#include <functional>

/* BOOST */
#include <boost/log/trivial.hpp>
#include <boost/asio.hpp>

#ifndef NDEBUG
#define MAX_LOG_LEVEL trace
#else
#define MAX_LOG_LEVEL info
#endif

#define LOG(lvl) \
    if constexpr (boost::log::trivial::lvl >= boost::log::trivial::MAX_LOG_LEVEL) \
        BOOST_LOG_TRIVIAL(lvl)

#define try_forever(X) try { while (true) { X } }

#define ignore_exception(X) try { X } catch (...) {}

using callback_t = std::function<void(void)>;
using new_connection_callback_t = std::function<void(std::shared_ptr<boost::asio::ip::tcp::socket>)>;

#endif //SIK_BOMBERMAN_UTILS_H
