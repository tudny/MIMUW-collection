#define _GNU_SOURCE
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <stdlib.h>
#include <inttypes.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <arpa/inet.h>
#include <endian.h>
#include <time.h>

#define print_err(...) fprintf(stderr, __VA_ARGS__)

#ifdef NDEBUG
const int debug = false;
#define print_debug(...)
#else
const int debug = true;
#define print_debug(...) fprintf(stdout, __VA_ARGS__)
#endif

#define SIZE(X) (sizeof(X) / sizeof((X)[0]))

// Evaluate `x`: if false, print error message and proper usage and exit with an error.
#define assure_usage(X)                                                   \
    do {                                                                  \
        if (!(X)) {                                                       \
            print_err("'"#X"' not assured.\n");                           \
            print_wrong_usage(argv[0]);                                   \
        }                                                                 \
    } while (0)

// Evaluate `x`: if false, print an error message and exit with an error.
#define ENSURE(x)                                                         \
    do {                                                                  \
        bool result = (x);                                                \
        if (!result) {                                                    \
            fprintf(stderr, "Error: %s was false in %s at %s:%d\n",       \
                #x, __func__, __FILE__, __LINE__);                        \
            exit(EXIT_FAILURE);                                           \
        }                                                                 \
    } while (0)

// Check if errno is non-zero, and if so, print an error message and exit with an error.
#define PRINT_ERRNO()                                                  \
    do {                                                               \
        if (errno != 0) {                                              \
            fprintf(stderr, "Error: errno %d in %s at %s:%d\n%s\n",    \
              errno, __func__, __FILE__, __LINE__, strerror(errno));   \
            exit(EXIT_FAILURE);                                        \
        }                                                              \
    } while (0)


// Set `errno` to 0 and evaluate `x`. If `errno` changed, describe it and exit.
#define CHECK_ERRNO(x)                                                             \
    do {                                                                           \
        errno = 0;                                                                 \
        (void) (x);                                                                \
        PRINT_ERRNO();                                                             \
    } while (0)

// Note: the while loop above wraps the statements so that the macro can be used with a semicolon
// for example: if (a) CHECK(x); else CHECK(y);

// ===================== CONSTANTS =====================

#define MESSAGE_ID_GET_EVENTS       1
#define MESSAGE_ID_EVENTS           2
#define MESSAGE_ID_GET_RESERVATION  3
#define MESSAGE_ID_RESERVATION      4
#define MESSAGE_ID_GET_TICKETS      5
#define MESSAGE_ID_TICKETS          6
#define MESSAGE_ID_BAD_REQUEST      255

#define COOKIE_SIZE 48

// ===================== PARSING PART =====================

const uint64_t PORT_LOWER = 0;
const uint64_t PORT_UPPER = 65535;
const uint64_t PORT_DEFAULT = 2022;

const uint64_t TIMEOUT_LOWER = 1;
const uint64_t TIMEOUT_UPPER = 86400;
const uint64_t TIMEOUT_DEFAULT = 5;

void print_wrong_usage(char *run_arg) {
    (void) run_arg;
    print_err("Wrong usage of program.\n");
    print_err("Usage: %s -f <file_name> [-p <port>] [-t <timeout>]\n", run_arg);
    exit(1);
}

bool parse_file_name(char *file_name, char **file_name_param) {
    *file_name_param = file_name;
    return (access(file_name, F_OK | R_OK) == 0);
}

bool is_digit(char c) {
    return '0' <= c && c <= '9';
}

bool has_only_digits(char *str) {
    for (; *str; str++)
        if (!is_digit(*str))
            return false;

    return true;
}

bool parse_uint64_t(char *str, uint64_t *value, uint64_t lower, uint64_t upper) {
    char *end_ptr;
    errno = 0;
    *value = strtoull(str, &end_ptr, 10);
    return (*end_ptr == '\0') && (*value >= lower) && (*value <= upper) && (errno == 0) && has_only_digits(str);
}

bool parse_port(char *port_str, uint64_t *port) {
    return parse_uint64_t(port_str, port, PORT_LOWER, PORT_UPPER);
}

bool parse_timeout(char *timeout_str, uint64_t *timeout) {
    return parse_uint64_t(timeout_str, timeout, TIMEOUT_LOWER, TIMEOUT_UPPER);
}

bool check_single_argv(char *arg) {
    for (; *arg; arg++)
        if (*arg != '-')
            return true;

    return false;
}

void check_strange_parameters(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        assure_usage(check_single_argv(argv[i]));
    }
}

typedef struct {
    char *file_name;
    uint64_t port;
    uint64_t timeout;
} parameters_t;

parameters_t load_params(int argc, char *argv[]) {
    bool file_param_is_present = false;

    check_strange_parameters(argc, argv);

    parameters_t parameters = {
            .port = PORT_DEFAULT,
            .timeout = TIMEOUT_DEFAULT
    };

    int opt;
    while ((opt = getopt(argc, argv, "f:p:t:")) != -1) {
        switch (opt) {
            case 'f':
                assure_usage(parse_file_name(optarg, &parameters.file_name));
                file_param_is_present = true;
                break;
            case 'p':
                assure_usage(parse_port(optarg, &parameters.port));
                break;
            case 't':
                assure_usage(parse_timeout(optarg, &parameters.timeout));
                break;
            default:
                print_wrong_usage(argv[0]);
                break;
        }
    }

    assure_usage(file_param_is_present);

    return parameters;
}

// ===================== PARSING PART : END =====================

// ===================== SAFE MEMORY =====================

void *require_not_null(void *ptr) {
    if (ptr == NULL) {
        print_err("Memory allocation error.");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void *safe_malloc(size_t memory_block_size) {
    return require_not_null(malloc(memory_block_size));
}

void *safe_calloc(size_t number_of_elements, size_t size_of_element) {
    return require_not_null(calloc(number_of_elements, size_of_element));
}

void *safe_realloc(void *ptr, size_t new_size) {
    return require_not_null(realloc(ptr, new_size));
}

void safe_free(void **ptr) {
    free(*ptr);
    *ptr = NULL;
}

// ===================== SAFE MEMORY : END =====================

// ===================== CUSTOM VECTOR PART =====================

const size_t VECTOR_INIT_SIZE = 4;

typedef struct {
    size_t size;
    size_t malloced_size;
    size_t nulled;
    void **data;
} vector_t;

vector_t *vector_new() {
    vector_t *vector = safe_malloc(sizeof(vector_t));
    memset(vector, 0, sizeof(vector_t));

    vector->malloced_size = VECTOR_INIT_SIZE;
    vector->size = vector->nulled = 0;
    vector->data = safe_malloc(VECTOR_INIT_SIZE * sizeof(void *));

    return vector;
}

void vector_free(vector_t *vector) {
    safe_free((void **) &vector->data);
    safe_free((void **) &vector);
}

void vector_double_size(vector_t *vector) {
    vector->malloced_size <<= 1;
    vector->data = safe_realloc(vector->data, vector->malloced_size * sizeof(void *));
}

void vector_add(vector_t *vector, void *element) {
    if (vector->size == vector->malloced_size)
        vector_double_size(vector);

    vector->data[vector->size++] = element;
}

// Only for reservations
void vector_clear_nulled(vector_t *vector);

// ===================== CUSTOM VECTOR PART : END =====================

// ===================== FILE PART =====================

const size_t DATAGRAM_MAX_SIZE = 65507;

typedef struct {
    uint32_t id;
    uint32_t ticket_count;
    uint16_t description_length;
    uint64_t prefix_size;
    char *description;
} event_t;

typedef struct __attribute__((__packed__)) {
    uint32_t event_id;
    uint16_t ticket_count;
    uint8_t description_length;
} event_info_t;

#define EV(C, X) ((event_t *) (C)->data[X])

void free_event(event_t *event) {
    safe_free((void **) &event->description);
    safe_free((void **) &event);
}

const size_t EVENT_DATA_SIZE = 4 + 2 + 1;
size_t events_datagram_size = 1;

vector_t *load_events(char *file_name) {
    FILE *events_file = fopen(file_name, "r");
    vector_t *events = vector_new();

    event_t *event = NULL;
    size_t line_number = 0;
    size_t len = 0;
    ssize_t read;
    char *line = NULL;
    char *end_ptr;

    while ((read = getline(&line, &len, events_file)) != EOF) {
        if (line_number & 1) { // Tickets
            errno = 0;
            event->ticket_count = strtoull(line, &end_ptr, 10);
            ENSURE(errno == 0);
            event->id = line_number >> 1;

            size_t additional_size = event->description_length + EVENT_DATA_SIZE;
            if (events_datagram_size + additional_size > DATAGRAM_MAX_SIZE) {
                free_event(event);
                break;
            }

            events_datagram_size += additional_size;

            vector_add(events, event);
        }
        else { // Description
            if (line[read - 1] == '\n')
                line[read - 1] = '\0';

            event = safe_malloc(sizeof(event_t));
            event->description = safe_malloc(sizeof(char) * (read + 1));
            memset(event->description, 0, sizeof(char) * (read + 1));
            strcpy(event->description, line);
            event->description_length = strlen(event->description);
            event->prefix_size = ((line_number >> 1) > 0) ? EV(events, (line_number >> 1) - 1)->prefix_size
                            + sizeof(event_info_t) + EV(events, (line_number >> 1) - 1)->description_length : 0;
        }

        ++line_number;
    }

    safe_free((void **) &line);
    fclose(events_file);

    return events;
}

typedef struct {
    size_t size;
    char *data;
} datagram_t;

// We can clear all descriptions which are too long to be sent.
datagram_t create_events_datagram(vector_t *events, size_t datagram_size) {
    datagram_t datagram = {
            .size = datagram_size,
            .data = safe_calloc(datagram_size, sizeof(char))
    };

    *datagram.data = MESSAGE_ID_EVENTS;

    char *datagram_ptr = datagram.data + 1;
    for (size_t i = 0; i < events->size; ++i) {
        event_info_t event_info = {
                .event_id = htonl(EV(events, i)->id),
                .ticket_count = htons(EV(events, i)->ticket_count),
                .description_length = EV(events, i)->description_length
        };

        *((event_info_t *) datagram_ptr) = event_info;
        datagram_ptr += sizeof(event_info_t);

        for (size_t d = 0; d < EV(events, i)->description_length; ++d) {
            *datagram_ptr = EV(events, i)->description[d];
            datagram_ptr++;
        }
    }

    return datagram;
}

datagram_t init_common_datagram() {
    datagram_t datagram = {
        .size = 0,
        .data = safe_calloc(DATAGRAM_MAX_SIZE, sizeof(char))
    };

    return datagram;
}

datagram_t clear_common_datagram(datagram_t datagram) {
    memset(datagram.data, 0, DATAGRAM_MAX_SIZE);
    datagram.size = 0;
    return datagram;
}

// ===================== FILE PART : END =====================

// ===================== SERVER =====================

int bind_socket(uint16_t port) {
    int socket_fd = socket(AF_INET, SOCK_DGRAM, 0); // creating IPv4 UDP socket
    ENSURE(socket_fd > 0);
    // after socket() call; we should close(sock) on any execution path;

    struct sockaddr_in server_address;
    server_address.sin_family = AF_INET; // IPv4
    server_address.sin_addr.s_addr = htonl(INADDR_ANY); // listening on all interfaces
    server_address.sin_port = htons(port);

    // bind the socket to a concrete address
    CHECK_ERRNO(bind(socket_fd, (struct sockaddr *) &server_address,
                     (socklen_t) sizeof(server_address)));

    return socket_fd;
}

size_t read_message(int socket_fd, struct sockaddr_in *client_address, char *buffer, size_t max_length) {
    socklen_t address_length = (socklen_t) sizeof(*client_address);
    int flags = 0; // we do not request anything special
    errno = 0;
    ssize_t len = recvfrom(socket_fd, buffer, max_length, flags,
                           (struct sockaddr *) client_address, &address_length);
    if (len < 0) {
        PRINT_ERRNO();
    }
    return (size_t) len;
}

void send_message(int socket_fd, const struct sockaddr_in *client_address, const char *message, size_t length) {
    socklen_t address_length = (socklen_t) sizeof(*client_address);
    int flags = 0;
    ssize_t sent_length = sendto(socket_fd, message, length, flags,
                                 (struct sockaddr *) client_address, address_length);
    ENSURE(sent_length == (ssize_t) length);
}

datagram_t events_datagram, common_datagram;
vector_t *events;
parameters_t parameters;

datagram_t handle_message(char *message, size_t message_length);

void run_server() {
    int socket_fd = bind_socket(parameters.port);

    char shared_buffer[DATAGRAM_MAX_SIZE + 1];

    print_debug("Started listening.\n");

    struct sockaddr_in client_address;
    size_t read_length;
    while (true) {
        read_length = read_message(socket_fd, &client_address, shared_buffer, sizeof(shared_buffer));
        print_debug("received %zd bytes from client %s:%u.\n",
                    read_length,
                    inet_ntoa(client_address.sin_addr),
                    ntohs(client_address.sin_port));

        datagram_t response = handle_message(shared_buffer, read_length);

        if (response.data != NULL) {
            send_message(socket_fd, &client_address, response.data, response.size);
        }
    }

    print_debug("Finished exchange.\n");

    CHECK_ERRNO(close(socket_fd));
}

// ===================== SERVER : END =====================

// ===================== TICKETS =====================

#define TICKET_SIZE 7

// 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ  <- available characters
const char * const alphabet = "4SBFOZHTQVXMJ20NPIW85L7CEKR6A1U9DY3G";

char common_ticket[TICKET_SIZE];
size_t first_free_ticket = 0;

char *generate_ticket(size_t number) {
    const size_t alphabet_size = strlen(alphabet);
    memset(common_ticket, 0, TICKET_SIZE);

    for (size_t i = 0; i < TICKET_SIZE; ++i) {
        common_ticket[i] = alphabet[((number % alphabet_size) + i) % alphabet_size];
        number /= alphabet_size;
    }

    return common_ticket;
}

// ===================== TICKETS : END =====================

// ===================== HANDLERS =====================

#define H_TO_N_64(X) (X) = htobe64(X)
#define N_TO_H_64(X) (X) = be64toh(X)
#define H_TO_N_32(X) (X) = htonl(X)
#define N_TO_H_32(X) (X) = ntohl(X)
#define H_TO_N_16(X) (X) = htons(X)
#define N_TO_H_16(X) (X) = ntohs(X)

datagram_t handle_bad_request(datagram_t datagram);

typedef struct __attribute__((__packed__)) {
    uint8_t message_id;
} get_events_t;

datagram_t handle_get_events(datagram_t datagram) {
    get_events_t *get_events = (get_events_t *) datagram.data;
    ENSURE(get_events->message_id == MESSAGE_ID_GET_EVENTS);
    ENSURE(datagram.size == sizeof(get_events_t));

    return events_datagram;
}

const size_t TICKETS_DATA_SIZE = 1 + 4 + 2;
const size_t SINGLE_TICKET_SIZE = 7;
const size_t RESERVATION_ID_LENGTH = 32;

vector_t *reservations;
uint32_t next_reservation_id = 1000000;

typedef struct __attribute__((__packed__)) {
    uint8_t message_id;
    uint32_t reservation_id;
    uint32_t event_id;
    uint16_t ticket_count;
    char cookie[COOKIE_SIZE];
    uint64_t expiration_time;
} reservation_t;

typedef enum {
    EXPIRED,
    AWAITING,
    COLLECTED
} reservation_status_t;

// Reservation wrapper
typedef struct {
    reservation_t *reservation;
    reservation_status_t status;
    size_t first_ticket;
} reservation_info_t;

#define RES(C, I) ((reservation_info_t *) (C)->data[I])

typedef struct __attribute__((__packed__)) {
    uint8_t message_id;
    uint32_t event_id;
    uint16_t ticket_count;
} get_reservation_t;

#define COOKIE_LOWER 33
#define COOKIE_UPPER 126

char random_cookie_char() {
    return (char) (rand() % (COOKIE_UPPER - COOKIE_LOWER + 1)) + COOKIE_LOWER;
}

void generate_cookie(reservation_t *reservation) {
    char *cookie = reservation->cookie;

    for (size_t i = 0; i < RESERVATION_ID_LENGTH; ++i) {
        cookie[i] = (reservation->reservation_id & (1l << i)) ? '1' : '0';
    }

    for (size_t i = RESERVATION_ID_LENGTH; i < COOKIE_SIZE; ++i) {
        cookie[i] = random_cookie_char();
    }
}

void update_event_ticket_count(uint32_t event_id) {
    char *event_datagram_ptr = events_datagram.data + 1 + EV(events, event_id)->prefix_size;
    ((event_info_t *) event_datagram_ptr)->ticket_count = htons(EV(events, event_id)->ticket_count);
}

void free_reservation(reservation_info_t *reservation_info) {
    safe_free((void **) &reservation_info->reservation);
    safe_free((void **) &reservation_info);
}

size_t reservation_to_check = 0;

// Only for reservations
void vector_clear_nulled(vector_t *vector) {

    if (vector->nulled * 2 <= vector->malloced_size)
        return;

    size_t new_size = vector->size - vector->nulled;
    size_t new_malloc = VECTOR_INIT_SIZE;
    while (new_malloc < new_size) new_malloc <<= 1;

    size_t new_iterator = 0;
    for (size_t i = 0; i < vector->size; ++i) {
        if (RES(reservations, i)->status != EXPIRED) {
            vector->data[new_iterator++] = vector->data[i];
        }
        else {
            free_reservation(RES(reservations, i));
        }
    }

    for (size_t i = new_iterator; i < vector->malloced_size; ++i) {
        vector->data[i] = NULL;
    }

    vector->malloced_size = new_malloc;
    vector->size = new_size;
    vector->nulled = 0;
    reservation_to_check = vector->size;

    vector->data = safe_realloc(vector->data, vector->malloced_size * sizeof(void *));
}

void prune_reservations() {
    uint64_t current_time = time(NULL);

    for (size_t i = reservation_to_check; i < reservations->size; ++i, ++reservation_to_check) {
        reservation_info_t *reservation_info = RES(reservations, i);
        if (reservation_info->status != AWAITING)
            continue;

        if (reservation_info->reservation->expiration_time <= current_time) {
            reservation_info->status = EXPIRED;
            ++reservations->nulled;

            EV(events, reservation_info->reservation->event_id)->ticket_count
                    += reservation_info->reservation->ticket_count;
            update_event_ticket_count(reservation_info->reservation->event_id);
        }
        else {
            break;
        }
    }

    vector_clear_nulled(reservations);
}

reservation_info_t *find_reservation(uint32_t reservation_id) {
    size_t begin = 0, midd,
            end = reservations->size + 1;

    while (end - begin > 1) {
        midd = (begin + end) / 2;
        if (RES(reservations, midd - 1)->reservation->reservation_id <= reservation_id)
            begin = midd;
        else
            end = midd;
    }

    if (begin == 0 || RES(reservations, begin - 1)->reservation->reservation_id != reservation_id)
        return NULL;

    return RES(reservations, begin - 1);
}

reservation_t *new_reservation(get_reservation_t *get_reservation) {
    reservation_t *reservation = safe_calloc(1, sizeof(reservation_t));
    reservation->message_id = MESSAGE_ID_RESERVATION;
    reservation->reservation_id = next_reservation_id++;
    reservation->event_id = get_reservation->event_id;
    reservation->ticket_count = get_reservation->ticket_count;
    generate_cookie(reservation);
    reservation->expiration_time = time(NULL) + parameters.timeout;

    reservation_info_t *reservation_info = safe_calloc(1, sizeof(reservation_info_t));
    reservation_info->reservation = reservation;
    reservation_info->status = AWAITING;
    vector_add(reservations, reservation_info);

    EV(events, reservation->event_id)->ticket_count -= reservation->ticket_count;
    update_event_ticket_count(reservation->event_id);

    return reservation;
}

datagram_t reserve_tickets(get_reservation_t *get_reservation) {
    reservation_t *reservation = new_reservation(get_reservation);

    H_TO_N_32(reservation->reservation_id);
    H_TO_N_32(reservation->event_id);
    H_TO_N_16(reservation->ticket_count);
    H_TO_N_64(reservation->expiration_time);

    common_datagram = clear_common_datagram(common_datagram);
    *((reservation_t *) common_datagram.data) = *reservation;
    common_datagram.size = sizeof(reservation_t);

    H_TO_N_32(reservation->reservation_id);
    H_TO_N_32(reservation->event_id);
    H_TO_N_16(reservation->ticket_count);
    H_TO_N_64(reservation->expiration_time);

    return common_datagram;
}

datagram_t handle_get_reservation(datagram_t datagram) {
    get_reservation_t *get_reservation = (get_reservation_t *) datagram.data;
    ENSURE(get_reservation->message_id == MESSAGE_ID_GET_RESERVATION);
    ENSURE(datagram.size == sizeof(get_reservation_t));

    N_TO_H_32(get_reservation->event_id);
    N_TO_H_16(get_reservation->ticket_count);

    print_debug("Asked for reservation: event_id=%d, ticket_count=%d\n", get_reservation->event_id, get_reservation->ticket_count);

    event_t *event = get_reservation->event_id >= events->size ? NULL : EV(events, get_reservation->event_id);

    // Tickets must fit into single datagram.
    size_t size_of_requested_tickets = TICKETS_DATA_SIZE + get_reservation->ticket_count * SINGLE_TICKET_SIZE;

    // No event of such id was found.
    if (!event || !get_reservation->ticket_count || get_reservation->ticket_count > event->ticket_count ||
            size_of_requested_tickets > DATAGRAM_MAX_SIZE) {

        print_debug("Request constraints not kept.\n");
        return handle_bad_request(datagram);
    }

    return reserve_tickets(get_reservation);
}

typedef struct __attribute__((__packed__)) {
    uint8_t message_id;
    uint32_t reservation_id;
    char cookie[COOKIE_SIZE];
} get_tickets_t;

const size_t TICKETS_DATAGRAM_SIZE = 1 + 4 + 2;

typedef struct __attribute__((__packed__)) {
    uint8_t message_id;
    uint32_t reservation_id;
    uint16_t ticket_count;
} tickets_t;

datagram_t handle_tickets(reservation_info_t *reservation_info) {
    common_datagram = clear_common_datagram(common_datagram);

    if (reservation_info->status != COLLECTED) {
        reservation_info->status = COLLECTED;
        reservation_info->first_ticket = first_free_ticket;
        first_free_ticket += reservation_info->reservation->ticket_count;
    }

    common_datagram.size = TICKETS_DATAGRAM_SIZE + reservation_info->reservation->ticket_count * TICKET_SIZE;
    tickets_t *tickets = (tickets_t *) common_datagram.data;
    tickets->message_id = MESSAGE_ID_TICKETS;
    tickets->ticket_count = reservation_info->reservation->ticket_count;
    tickets->reservation_id = reservation_info->reservation->reservation_id;

    H_TO_N_32(tickets->reservation_id);
    H_TO_N_16(tickets->ticket_count);

    char *tickets_ptr = common_datagram.data + sizeof(tickets_t);
    for (size_t i = 0; i < reservation_info->reservation->ticket_count; ++i) {
        char *ticket = generate_ticket(reservation_info->first_ticket + i);
        for (size_t d = 0; d < TICKET_SIZE; ++d) {
            *tickets_ptr = ticket[d];
            ++tickets_ptr;
        }
    }

    return common_datagram;
}

bool cookie_matches(const char *cookie_a, const char *cookie_b) {
    for (size_t i = 0; i < COOKIE_SIZE; ++i)
        if (cookie_a[i] != cookie_b[i])
            return false;

    return true;
}

datagram_t handle_get_tickets(datagram_t datagram) {
    get_tickets_t *get_tickets = (get_tickets_t *) datagram.data;
    ENSURE(get_tickets->message_id == MESSAGE_ID_GET_TICKETS);
    ENSURE(datagram.size == sizeof(get_tickets_t));

    N_TO_H_32(get_tickets->reservation_id);

    print_debug("Asked for tickets: reservation_id=%d\n", get_tickets->reservation_id);

    reservation_info_t *reservation_info = find_reservation(get_tickets->reservation_id);
    if (reservation_info == NULL ||
          !cookie_matches(get_tickets->cookie, reservation_info->reservation->cookie) ||
          reservation_info->status == EXPIRED) {

        return handle_bad_request(datagram);
    }

    return handle_tickets(reservation_info);
}

typedef struct __attribute__((__packed__)) {
    uint8_t message_id;
    uint32_t id;
} bad_request_t;

datagram_t handle_bad_request(datagram_t datagram) {
    common_datagram = clear_common_datagram(common_datagram);

    uint32_t id = *((uint32_t *) (datagram.data + 1));

    bad_request_t *bad_request = (bad_request_t *) common_datagram.data;
    bad_request->message_id = MESSAGE_ID_BAD_REQUEST;
    bad_request->id = H_TO_N_32(id);
    common_datagram.size = sizeof(bad_request_t);

    return common_datagram;
}

typedef struct {
    uint8_t id;
    datagram_t (*handler)(datagram_t);
    size_t expected_datagram_size;
} handler_t;

handler_t handlers[] = {
        {MESSAGE_ID_GET_EVENTS,         handle_get_events, sizeof(get_events_t)},
        {MESSAGE_ID_GET_RESERVATION,    handle_get_reservation, sizeof(get_reservation_t)},
        {MESSAGE_ID_GET_TICKETS,        handle_get_tickets, sizeof(get_tickets_t)},
};

#define SIZE(X) (sizeof(X) / sizeof((X)[0]))

datagram_t handle_message(char *message, size_t message_length) {
    datagram_t datagram = {
            .data = message,
            .size = message_length
    };

    prune_reservations();

    for (size_t i = 0; i < SIZE(handlers); ++i) {
        if (handlers[i].id == *datagram.data) {
            if (datagram.size == handlers[i].expected_datagram_size) {
                return handlers[i].handler(datagram);
            }
            else {
                // We found the appropriate handler, but the size of the datagram is wrong.
                break;
            }
        }
    }

    datagram_t result = {
            .size = 0,
            .data = NULL
    };

    return result;
}

// ===================== HANDLERS : END =====================

int main(int argc, char *argv[]) {

    parameters = load_params(argc, argv);
    print_debug("Loaded with parameters: filename='%s', port=%lu, timeout=%lu\n",
                parameters.file_name, parameters.port, parameters.timeout);

    events = load_events(parameters.file_name);
    events_datagram = create_events_datagram(events, events_datagram_size);

    reservations = vector_new();
    common_datagram = init_common_datagram();

    run_server();

    for (size_t i = 0; i < events->size; ++i) {
        free_event(events->data[i]);
    }

    vector_free(events);
    safe_free((void **) &events_datagram.data); // Free events datagram
    safe_free((void **) &common_datagram.data);

    for (size_t i = 0; i < reservations->size; ++i) {
        free_reservation(RES(reservations, i));
    }

    vector_free(reservations);

    return 0;
}
