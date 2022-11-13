#ifndef SIK_BOMBERMAN_GUISERVERBRIDGE_H
#define SIK_BOMBERMAN_GUISERVERBRIDGE_H

/* LOCAL */
#include "../Bridge.h"

/**
 * This class represents an implementation
 * of synchronized connection between the GUI
 * and the Server.
 */
class GuiServerBridge : public Bridge {
public:
    /**
     * Constructor of the bridge synchronizer.
     * @param gui_connection connection to GUI
     * @param server_connection connection to Server
     * @param callback callback function to call on exception
     * @param game_state state of the game shared between all threads
     * @param player_name name of the player
     */
    GuiServerBridge(GUI::Connection &gui_connection,
                    Server::Connection &server_connection,
                    const callback_t &callback,
                    GameState &game_state,
                    const std::string &player_name);

    /**
     * Job for thread to run. It is a infinite loop, but
     * when an exception occurs a callback function is
     * called.
     */
    void operator()() override;

private:
    /**
     * Players name is being send after each game in
     * a JOIN message.
     */
    const std::string &player_name;
};


#endif //SIK_BOMBERMAN_GUISERVERBRIDGE_H
