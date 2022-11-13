#include <stdio.h>
#include <unistd.h>

int main ()
{
    printf("Hello from %d\n", getpid());

    if (fork() == -1) printf("Error in fork\n");  /* powstaje nowy proces */
    printf("Inside %d\n", getpid());

    sleep(60);

    printf("Goodbye from %d\n", getpid());

    return 0;
}

