#include "../fuzzy.h"
#include "../set_utils.h"
#include "../utils.h"
#include "../logger.h"
#include <cassert>

int main() {

    {
        TriFuzzyNumSet s1 = SetUtils::random_set(10);

        TriFuzzyNum a(10, 20, 30);
        TriFuzzyNum b(40, 50, 60);
        TriFuzzyNum c(70, 80, 90);

        s1.insert(a);
        s1.insert(std::move(b));
        s1.insert(TriFuzzyNum(100, 110, 120));

        s1.remove(c);
        s1.remove(a);
        s1.remove(b);
    }

    {
        TriFuzzyNumSet s1;

        try {
            TriFuzzyNum mean = s1.arithmetic_mean();
            std::cout << mean << "\n";
            assert("Never been here");
        } catch(const std::length_error &error) {
            Logger::log(error.what());
        }
    }

    return 0;
}
