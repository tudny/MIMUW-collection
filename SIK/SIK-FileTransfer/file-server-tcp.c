#include <sys/types.h>
#include <netinet/in.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "common.h"

#define LINE_SIZE 100
#define QUEUE_LENGTH 5

#define pthread_printf(...) printf("[thread %lu, pid %d] ", (unsigned long) pthread_self(), getpid()); \
        printf(__VA_ARGS__)

char *bools[] = { "FALSE", "TRUE" };

pthread_mutex_t mutex;
uint64_t total_size = 0;

void init_mutex() {
    CHECK(pthread_mutex_init(&mutex, 0));
}

void rem_mutex() {
    CHECK(pthread_mutex_destroy(&mutex));
}

void add_to_global(uint64_t value) {
    CHECK(pthread_mutex_lock(&mutex));
    total_size += value;
    pthread_printf("total size of uploaded files %lu\n", total_size);
    CHECK(pthread_mutex_unlock(&mutex));
}

void *handle_connection(void *client_fd_ptr) {
    int client_fd = *(int *) client_fd_ptr;
    free(client_fd_ptr);

    char *ip = get_ip_from_socket(client_fd);
    int port = get_port_from_socket(client_fd);

    struct file_info_t file_info;
    size_t read = receive_message(client_fd, &file_info, sizeof(file_info), MSG_WAITALL);
    ENSURE(read == sizeof(file_info));

    file_info.file_name_size = ntohs(file_info.file_name_size);
    file_info.file_size = be64toh(file_info.file_size);

    char file_name[file_info.file_name_size + 1];
    memset(file_name, 0, sizeof(file_name));
    read = receive_message(client_fd, &file_name, sizeof(file_name) - 1, MSG_WAITALL);
    ENSURE(read + 1 == sizeof(file_name));

    pthread_printf("new client %s:%d size=%lu file=%s\n", ip, port, file_info.file_size, file_name);
    sleep(1);

    FILE *file = fopen(file_name, "wb");
    uint64_t total_size_read = 0;

    char line[LINE_SIZE];
    for (;;) {
        memset(line, 0, sizeof(line));
        read = receive_message(client_fd, line, sizeof(line) - 1, MSG_WAITALL);
        if (read == 0)
            break;

        total_size_read += read;
        fwrite(line, 1, read, file);
    }

    CHECK_ERRNO(fclose(file));

    pthread_printf("client %s:%d has sent its file of size=%lu\n", ip, port, total_size_read);
    pthread_printf("sent file is of expected size=%s\n", bools[total_size_read == file_info.file_size]);

    add_to_global(total_size_read);

    pthread_printf("connection closed\n");
    CHECK_ERRNO(close(client_fd));
    return 0;
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
        fatal("Usage: %s <port>", argv[0]);
    }

    init_mutex();

    uint16_t port = read_port(argv[1]);

    int socket_fd = open_socket();

    bind_socket(socket_fd, port);
    printf("Listening on port %d\n", port);

    start_listening(socket_fd, QUEUE_LENGTH);

    for (;;) {
        struct sockaddr_in client_addr;

        int client_fd = accept_connection(socket_fd, &client_addr);

        // Arguments for the thread must be passed by pointer
        int *client_fd_pointer = malloc(sizeof(int));
        if (client_fd_pointer == NULL) {
            fatal("malloc");
        }
        *client_fd_pointer = client_fd;

        pthread_t thread;
        CHECK_ERRNO(pthread_create(&thread, 0, handle_connection, client_fd_pointer));
        CHECK_ERRNO(pthread_detach(thread));
    }

    rem_mutex();
}

