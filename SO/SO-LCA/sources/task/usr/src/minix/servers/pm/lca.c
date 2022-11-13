#include "pm.h"
#include <minix/callnr.h>
#include <sys/types.h>
#include <stdio.h>
#include "mproc.h"

int is_in_use(size_t index) {
    return mproc[index].mp_flags & IN_USE;
}

void fill_path(size_t index, pid_t *path, size_t *size) {
    *size = 0;
    while (is_in_use(index)) {
        path[(*size)++] = mproc[index].mp_pid;

        size_t parent_id = mproc[index].mp_parent;
        if (mproc[index].mp_pid == mproc[parent_id].mp_pid) break;
        index = parent_id;
    }
}

int calc_lcs(size_t proc1_counter, pid_t *proc1_path, size_t proc2_counter, pid_t *proc2_path) {
    if (proc1_counter < proc2_counter) {
        return calc_lcs(proc2_counter, proc2_path, proc1_counter, proc1_path);
    }

    if (proc1_counter == 1 || proc2_counter == 1) {
        /* one of them must be the root, so it has no ancestors */
        return ESRCH;
    }

    for (size_t i = 0; i < proc1_counter; ++i) {
        if (proc1_path[i] == proc2_path[0]) {
            /* if one is the ancestor of the other, then next is an ancestor of both */
            return proc1_path[i + 1];
        }
    }

    while (proc1_counter != 0 && proc2_counter != 0) {
        if (proc1_path[proc1_counter - 1] == proc2_path[proc2_counter - 1]) {
            --proc1_counter;
            --proc2_counter;
        }
        else {
            return proc1_path[proc1_counter];
        }
    }

    return ESRCH;
}

int do_getlcapid(void) {

    pid_t pid1 = m_in.m_u32.data[0];
    pid_t pid2 = m_in.m_u32.data[1];

    size_t proc1, proc2;
    int was_proc1_set = 0, was_proc2_set = 0;

    for (size_t i = 0; i < NR_PROCS; ++i) {
        if (mproc[i].mp_pid == pid1) {
            proc1 = i;
            was_proc1_set = 1;
        }
        if (mproc[i].mp_pid == pid2) {
            proc2 = i;
            was_proc2_set = 1;
        }
    }

    if (!was_proc1_set || !was_proc2_set || !is_in_use(proc1) || !is_in_use(proc2)) {
        return EINVAL;
    }

    size_t proc1_counter = 0;
    size_t proc2_counter = 0;
    pid_t proc1_path[NR_PROCS];
    pid_t proc2_path[NR_PROCS];

    fill_path(proc1, proc1_path, &proc1_counter);
    fill_path(proc2, proc2_path, &proc2_counter);

    return calc_lcs(proc1_counter, proc1_path, proc2_counter, proc2_path);
}
