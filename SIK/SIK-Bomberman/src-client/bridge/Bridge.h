#ifndef SIK_BOMBERMAN_BRIDGE_H
#define SIK_BOMBERMAN_BRIDGE_H

/* LOCAL */
#include "../gui/Connection.h"
#include "../server/Connection.h"
#include "../utils/utils.h"
#include "../game/GameState.h"

/**
 * This class represents the synchronization between
 * Server and GUI. The object require GUI and Server
 * connection objects to perform operations. It runs
 * in an infinite loop, so it should be run as a thread
 * job. As the job runs in an infinite loop when
 * an exception occurs (e.g. socket is closed) callback
 * function should be called.
 */
class Bridge {
public:

    /**
     * Constructor of the bridge synchronizer.
     * @param gui_connection connection to GUI
     * @param server_connection connection to Server
     * @param callback callback function to call on exception
     * @param game_state state of the game shared between all threads
     */
    Bridge(
        GUI::Connection &gui_connection,
        Server::Connection &server_connection,
        const callback_t &callback,
        GameState &game_state);

    /**
     * This is a job to be done by thread.
     */
    virtual void operator()() = 0;

    /**
     * Clone all connections.
     */
    void close();

protected:
    GUI::Connection &gui_connection_; ///< Connection with GUI.
    Server::Connection &server_connection_; ///< Connection with Server.
    const callback_t &callback_; ///< Callback after finishing.
    GameState &game_state; ///< State of the game.
};


#endif //SIK_BOMBERMAN_BRIDGE_H
