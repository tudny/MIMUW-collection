#include "../fuzzy.h"
#include "../logger.h"

int main() {

    Logger::log("Checking if classes are properly included.");

    TriFuzzyNum number(1.0, 2.0, 3.0);
    TriFuzzyNumSet set;

    Logger::log("They are.");


    return 0;
}
