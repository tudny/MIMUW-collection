#ifndef SIK_BOMBERMAN_SENDER_H
#define SIK_BOMBERMAN_SENDER_H


#include "Connection.h"
#include "PlayerInstance.h"

namespace Client {

    class Sender {
    public:
        explicit Sender(std::shared_ptr<Client::Connection> connection, std::shared_ptr<PlayerInstance> instance);

        void operator()();

    private:
        std::shared_ptr<Client::Connection> connection;
        std::shared_ptr<PlayerInstance> instance;
    };
}


#endif //SIK_BOMBERMAN_SENDER_H
