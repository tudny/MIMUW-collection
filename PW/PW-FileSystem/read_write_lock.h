#pragma once

typedef struct rw_lock_t rw_lock_t;

rw_lock_t *rw_lock_new();

int rw_lock_free(rw_lock_t *lock);

int read_lock(rw_lock_t *lock);

int read_unlock(rw_lock_t *lock);

int write_lock(rw_lock_t *lock);

int writer_unlock(rw_lock_t *lock);
