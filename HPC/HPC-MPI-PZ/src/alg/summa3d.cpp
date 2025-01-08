#include "summa3d.hpp"
#include "summa2d.hpp"
#include "../utils/tags.hpp"
#include "../driver.hpp"
#include "../worker.hpp"

std::vector<uint32_t> calc_bounds(
        uint32_t l,
        uint32_t n,
        uint32_t m
) {
    Topology local_topology{l * l, std::nullopt};
    std::vector<uint32_t> bounds;
    for (uint32_t i = 0; i < l; i++) {
        auto [start, end] = local_topology.tell_my_range(0, i, n, m, 1);
        bounds.push_back(start);
    }
    bounds.push_back(m);
    return bounds;
}

std::vector<CompressedSparseRowsMatrixPtr> all_to_all(
        const Topology &topology,
        const ProcessIdentifier &id,
        std::vector<CompressedSparseRowsMatrixPtr> &Df_k
) {
    auto &[id_l, id_i, id_j] = id;

    std::vector<CompressedSparseRowsMatrixPtr> Cf_k(topology.l);
    Cf_k[id_l] = Df_k[id_l];

    for (uint32_t offset_l = 1; offset_l < topology.l; offset_l++) {
        auto send_to_layer = (id_l + offset_l) % topology.l;
        auto recv_from_layer = (id_l - offset_l + topology.l) % topology.l;

        ProcessIdentifier send_to{send_to_layer, id_i, id_j};
        ProcessIdentifier recv_from{recv_from_layer, id_i, id_j};

        PackedMatrixMetadata metadata = Df_k[send_to_layer]->metadata();
        RequestsPtr requests = send_sub_matrix(topology, *Df_k[send_to_layer], metadata, send_to);

        auto Cf_k_recv = recv_sub_matrix_common(topology, id, recv_from);
        Cf_k[recv_from_layer] = Cf_k_recv;

        MPI_Waitall(static_cast<int>(requests->size()), requests->data(), MPI_STATUSES_IGNORE);
    }

    return Cf_k;
}

CompressedSparseRowsMatrixPtr run_summa3d(
        const MatmulParams::ParamsPtr &params,
        const Topology &topology,
        const CompressedSparseRowsMatrixPtr &A,
        const CompressedSparseRowsMatrixPtr &B,
        int rank, int size
) {
    ProcessIdentifier id = topology.get_process_identifier(static_cast<uint32_t>(rank));
    const auto &[id_l, id_i, id_j] = id;

    auto D_k = run_summa2d_without_converting(params, topology, A, B, rank, size, id_l);
    auto &[_, n, m] = *D_k;
    auto Df_k = convert_to_n_csr_s(*D_k, topology.l, calc_bounds(topology.l, n, m));
    auto Cf_k = all_to_all(topology, id, Df_k);
    auto C_k = add_sub_csr_s_hashing(Cf_k);
    return C_k;
}
