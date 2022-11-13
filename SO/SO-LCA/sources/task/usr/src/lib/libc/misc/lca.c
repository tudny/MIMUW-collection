#include <lib.h>
#include <minix/rs.h>

int get_mp_endpt(endpoint_t *pm_ep) {
    return minix_rs_lookup("pm", pm_ep);
}

pid_t getlcapid(pid_t pid_1, pid_t pid_2) {
    endpoint_t pm_ep;
    message m;
    m.m_u32.data[0] = pid_1;
    m.m_u32.data[1] = pid_2;

    if (get_mp_endpt(&pm_ep) != 0) {
        errno = ENOSYS;
        return -1;
    }

    return (_syscall(pm_ep, PM_GETLCAPID, &m));
}
