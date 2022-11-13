#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include "read_write_lock.h"

#define P(X) pthread_mutex_lock(X)
#define V(X) pthread_mutex_unlock(X)
#define wait(X, Y) pthread_cond_wait(X, Y)
#define signal_all(X) pthread_cond_broadcast(X)

struct rw_lock_t {
    pthread_mutex_t mutex;
    pthread_cond_t readers;
    pthread_cond_t writers;
    int readers_count;
    int readers_waiting;
    int writers_count;
    int writers_waiting;
    int change;
};

rw_lock_t *rw_lock_new() {
    rw_lock_t *lock = malloc(sizeof(rw_lock_t));
    if (lock == NULL)
        return NULL; // Memory error

    if (pthread_mutex_init(&lock->mutex, 0) != 0) {
        free(lock);
        return NULL;
    }

    if (pthread_cond_init(&lock->readers, 0) != 0) {
        pthread_mutex_destroy(&lock->mutex);
        free(lock);
        return NULL;
    }

    if (pthread_cond_init(&lock->writers, 0) != 0) {
        pthread_mutex_destroy(&lock->mutex);
        pthread_cond_destroy(&lock->readers);
        free(lock);
        return NULL;
    }

    lock->readers_count = lock->writers_count = lock->readers_waiting = lock->writers_waiting = lock->change = 0;

    return lock;
}

int rw_lock_free(rw_lock_t *lock) {
    int err;
    if ((err = pthread_mutex_destroy(&lock->mutex)) != 0) {
        return err;
    }
    if ((err = pthread_cond_destroy(&lock->readers)) != 0) {
        return err;
    }
    if ((err = pthread_cond_destroy(&lock->writers)) != 0) {
        return err;
    }

    free(lock);
    return 0;
}

int read_lock(rw_lock_t *lock) {
    int err;
    if ((err = P(&lock->mutex)) != 0) {
        return err;
    }

    ++lock->readers_waiting;

    if (lock->writers_count > 0 || lock->writers_waiting != 0) {
        do {
            if ((err = wait(&lock->readers, &lock->mutex)) != 0) {
                return err;
            }
        } while (lock->change == 0);
        --lock->change;
    }

    assert(lock->readers_waiting > 0); // TODO safety assertion
    --lock->readers_waiting;
    assert(lock->writers_count == 0);
    ++lock->readers_count;

    if ((err = V(&lock->mutex)) != 0) {
        return err;
    }

    return 0;
}

int read_unlock(rw_lock_t *lock) {
    int err;
    if ((err = P(&lock->mutex)) != 0) {
        return err;
    }

    assert(lock->readers_count > 0);
    --lock->readers_count;
    if (lock->readers_count == 0) {
        if ((err = signal_all(&lock->writers)) != 0) {
            return err;
        }
    }

    if ((err = V(&lock->mutex)) != 0) {
        return err;
    }

    return 0;
}

int write_lock(rw_lock_t *lock) {
    int err;
    if ((err = P(&lock->mutex)) != 0) {
        return err;
    }

    ++lock->writers_waiting;

    while (lock->readers_count + lock->writers_count > 0 || lock->change > 0) {
        if ((err = wait(&lock->writers, &lock->mutex)) != 0) {
            return err;
        }
    }

    assert(lock->writers_waiting > 0); // TODO safety assertion
    --lock->writers_waiting;
    assert(lock->readers_count == 0);
    assert(lock->writers_count == 0);
    ++lock->writers_count;

    if ((err = V(&lock->mutex)) != 0) {
        return err;
    }

    return 0;
}

int writer_unlock(rw_lock_t *lock) {
    int err;
    if ((err = P(&lock->mutex)) != 0) {
        return err;
    }

    assert(lock->writers_count > 0); // TODO safety assertion

    --lock->writers_count;
    if (lock->readers_waiting > 0) {
        lock->change = lock->readers_waiting;
        if ((err = signal_all(&lock->readers)) != 0) {
            return err;
        }
    }
    else {
        if ((err = signal_all(&lock->writers)) != 0) {
            return err;
        }
    }

    if ((err = V(&lock->mutex)) != 0) {
        return err;
    }

    return 0;
}
