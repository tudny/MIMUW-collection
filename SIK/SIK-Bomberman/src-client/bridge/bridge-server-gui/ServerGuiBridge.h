#ifndef SIK_BOMBERMAN_SERVERGUIBRIDGE_H
#define SIK_BOMBERMAN_SERVERGUIBRIDGE_H

/* LOCAL */
#include "../Bridge.h"

/**
 * This class represents an implementation
 * of synchronized connection between the Server
 * and the GUI.
 */
class ServerGuiBridge : public Bridge {
public:

    /**
     * Constructor of the bridge synchronizer.
     * @param gui_connection connection to GUI
     * @param server_connection connection to Server
     * @param callback callback function to call on exception
     * @param game_state state of the game shared between all threads
     */
    ServerGuiBridge(GUI::Connection &gui_connection,
                    Server::Connection &server_connection,
                    const callback_t &callback,
                    GameState &game_state);

    /**
     * Job for thread to run. It is a infinite loop, but
     * when an exception occurs a callback function is
     * called.
     */
    void operator()() override;
};


#endif //SIK_BOMBERMAN_SERVERGUIBRIDGE_H
