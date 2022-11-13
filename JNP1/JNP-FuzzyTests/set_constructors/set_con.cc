#include "../fuzzy.h"
#include "../utils.h"
#include "../logger.h"

int main() {

    {
        TriFuzzyNumSet set1;
        TriFuzzyNumSet set2;
        TriFuzzyNumSet set3 ({});
        TriFuzzyNumSet set4 ({TriFuzzyNum(10, 20, 30)});
        TriFuzzyNumSet set5 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60)});
        TriFuzzyNumSet set6 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90)});
        TriFuzzyNumSet set7 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120)});
        TriFuzzyNumSet set8 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150)});
        TriFuzzyNumSet set9 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180)});
        TriFuzzyNumSet set10 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210)});
        TriFuzzyNumSet set11 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240)});
        TriFuzzyNumSet set12 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270)});
        TriFuzzyNumSet set13 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300)});
        TriFuzzyNumSet set14 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330)});
        TriFuzzyNumSet set15 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360)});
        TriFuzzyNumSet set16 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390)});
        TriFuzzyNumSet set17 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390), TriFuzzyNum(400, 410, 420)});
        TriFuzzyNumSet set18 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390), TriFuzzyNum(400, 410, 420), TriFuzzyNum(430, 440, 450)});
        TriFuzzyNumSet set19 ({TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390), TriFuzzyNum(400, 410, 420), TriFuzzyNum(430, 440, 450), TriFuzzyNum(460, 470, 480)});
    }

    {
        TriFuzzyNumSet set1;
        TriFuzzyNumSet set2;
        TriFuzzyNumSet set3 = {};
        TriFuzzyNumSet set4 = {TriFuzzyNum(10, 20, 30)};
        TriFuzzyNumSet set5 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60)};
        TriFuzzyNumSet set6 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90)};
        TriFuzzyNumSet set7 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120)};
        TriFuzzyNumSet set8 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150)};
        TriFuzzyNumSet set9 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180)};
        TriFuzzyNumSet set10 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210)};
        TriFuzzyNumSet set11 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240)};
        TriFuzzyNumSet set12 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270)};
        TriFuzzyNumSet set13 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300)};
        TriFuzzyNumSet set14 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330)};
        TriFuzzyNumSet set15 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360)};
        TriFuzzyNumSet set16 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390)};
        TriFuzzyNumSet set17 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390), TriFuzzyNum(400, 410, 420)};
        TriFuzzyNumSet set18 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390), TriFuzzyNum(400, 410, 420), TriFuzzyNum(430, 440, 450)};
        TriFuzzyNumSet set19 = {TriFuzzyNum(10, 20, 30), TriFuzzyNum(40, 50, 60), TriFuzzyNum(70, 80, 90), TriFuzzyNum(100, 110, 120), TriFuzzyNum(130, 140, 150), TriFuzzyNum(160, 170, 180), TriFuzzyNum(190, 200, 210), TriFuzzyNum(220, 230, 240), TriFuzzyNum(250, 260, 270), TriFuzzyNum(280, 290, 300), TriFuzzyNum(310, 320, 330), TriFuzzyNum(340, 350, 360), TriFuzzyNum(370, 380, 390), TriFuzzyNum(400, 410, 420), TriFuzzyNum(430, 440, 450), TriFuzzyNum(460, 470, 480)};
    
        TriFuzzyNumSet set20(set1);
        TriFuzzyNumSet set21(set2);
        TriFuzzyNumSet set22(set3);
        TriFuzzyNumSet set23(set4);
        TriFuzzyNumSet set24(set5);

        TriFuzzyNumSet set25(std::move(set1));
        TriFuzzyNumSet set26(std::move(set2));
        TriFuzzyNumSet set27(std::move(set3));
        TriFuzzyNumSet set28(std::move(set4));
        TriFuzzyNumSet set29(std::move(set5));

        set20 = set25;
        set21 = set26;
        set22 = set27;
        set23 = set28;
        set24 = set29;

        set20 = std::move(set25);
        set21 = std::move(set26);
        set22 = std::move(set27);
        set23 = std::move(set28);
        set24 = std::move(set29);
    }

    return 0;
}
