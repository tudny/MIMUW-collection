#include <unistd.h>
#include <stdio.h>

int main(){
    printf("Hello from synced test! time!\n");

    pid_t pidd = getlcapid(196, 248);

    printf("Result=%d\n", pidd);

    printf("Hello from the end!\n");

    return 0;
}
