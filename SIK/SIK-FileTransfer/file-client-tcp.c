#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <arpa/inet.h>
#define _FILE_OFFSET_BITS 64
#include <sys/stat.h>
#include <libgen.h>

#include "err.h"
#include "common.h"

#define BUFFER_SIZE 65536

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fatal("Usage: %s host port file_name \n", argv[0]);
    }

    char *host = argv[1];
    uint16_t port = read_port(argv[2]);
    char *file_path = argv[3];

    FILE *file = fopen(file_path, "rb");
    if (!file) {
        fatal("No such file! %s", file_path);
    }

    char *file_name = basename(file_path);

    struct stat st;
    if (stat(file_path, &st) != 0) {
        perror("Error calling stat.");
    }

    uint64_t file_size = st.st_size;

    struct sockaddr_in server_address = get_address(host, port);

    int socket_fd = open_socket();

    connect_socket(socket_fd, &server_address);

    char *server_ip = get_ip(&server_address);
    uint16_t server_port = get_port(&server_address);

    printf("Connected to %s:%d\n", server_ip, server_port);

    char buffer[BUFFER_SIZE];
    memset(buffer, 0, BUFFER_SIZE);

    uint16_t file_name_size = (uint16_t) strnlen(file_name, UINT16_MAX);

    struct file_info_t info;
    info.file_name_size = htons(file_name_size);
    info.file_size = htobe64(file_size);

    printf("Sending file '%.*s' of size=%lu\n", file_name_size, file_name, file_size);

    send_message(socket_fd, &info, sizeof(struct file_info_t), NO_FLAGS);
    send_message(socket_fd, file_name, file_name_size, NO_FLAGS);

    size_t read_len;
    while ((read_len = fread(buffer, 1, sizeof(buffer), file)) > 0) {
        send_message(socket_fd, buffer, read_len, NO_FLAGS);
        printf("Sent %zu bytes to %s:%d.\n", read_len, server_ip, server_port);
    }

    fclose(file);
}
