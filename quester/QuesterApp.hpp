/**
 *
 *
 *
 */

#ifndef QUESTERAPP_HPP
#define QUESTERAPP_HPP

namespace quester {

class QuesterApp {
public:
    QuesterApp (int argc, char** argv);
    void run();
private:
    void init ();
};

}

#endif



