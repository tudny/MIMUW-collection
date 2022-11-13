#ifndef SIK_BOMBERMAN_RECEIVER_H
#define SIK_BOMBERMAN_RECEIVER_H


#include "Connection.h"
#include "PlayerInstance.h"

namespace Client {

    class Receiver {
    public:
        explicit Receiver(const std::shared_ptr<Client::Connection> &connection,
                          const std::shared_ptr<PlayerInstance> &instance);

        void operator()();

        std::shared_ptr<Client::Connection> connection;
    private:
        std::shared_ptr<PlayerInstance> instance;
    };
}


#endif //SIK_BOMBERMAN_RECEIVER_H
