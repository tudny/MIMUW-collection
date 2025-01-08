from pathlib import Path
import os
from time import sleep

SBATCH_PATH = 'sbatch'
TESTING_DATA_PATH = '/home/tudny/HPC-MPI/testing-data'

CASES = [
    ("small-0", 12),
    ("medium-0", 120),
    ("large-0", 125),
]

GRIDS = [
    (4, 1, 1),
    (9, 1, 1),
    (16, 1, 1),
    (12, 3, 1),
    (24, 6, 1),
    (24, 6, 4),
    (24, 6, 9),
]

ALGS = [
    "2D",
    "3D",
]


def create_sbatch_script(
        case_name, case_g,
        grid_tasks, grid_nodes, grid_layers,
        alg
):
    run_name = 'matmul_{}_t{}_n{}_l{}_a{}'.format(case_name, grid_tasks, grid_nodes, grid_layers, alg)
    a_matrix_file = '{}/{}_A.csr'.format(TESTING_DATA_PATH, case_name)
    b_matrix_file = '{}/{}_B.csr'.format(TESTING_DATA_PATH, case_name)
    sbatch_script = '''#!/bin/bash
#SBATCH --job-name {}
#SBATCH --output "output/{}-%j.out"
#SBATCH --error "output/{}-%j.err"
#SBATCH --account "g96-1905"
#SBATCH --nodes {}
#SBATCH --tasks-per-node {}
#SBATCH --time 00:01:30

srun build/matmul -a {} -b {} -g {} -l {} -t {}
    '''.format(run_name, run_name, run_name, grid_nodes, grid_tasks, a_matrix_file, b_matrix_file, case_g, grid_layers,
               alg)

    sbatch_file = ('{}/{}.sbatch'.format(SBATCH_PATH, run_name))
    with open(sbatch_file, 'w') as f:
        f.write(sbatch_script)
    return sbatch_file


def is_perfect_square(n):
    return n == int(n ** 0.5) ** 2


def main():
    Path(SBATCH_PATH).mkdir(parents=True, exist_ok=True)
    for grid_tasks, grid_nodes, grid_layers in GRIDS:
        np_2d = grid_tasks * grid_nodes
        assert is_perfect_square(np_2d)
        assert np_2d % grid_layers == 0
        assert is_perfect_square(np_2d // grid_layers)

    for case_name, case_g in CASES:
        for grid_tasks, grid_nodes, grid_layers in GRIDS:
            for alg in ALGS:
                sbatch_script = create_sbatch_script(
                    case_name, case_g,
                    grid_tasks, grid_nodes, grid_layers,
                    alg
                )
                print('[NEW SBATCH SCRIPT]')
                print('  CARD: {}_{}'.format(case_name, case_g))
                print('  GRID: {} tasks, {} nodes, {} layers'.format(grid_tasks, grid_nodes, grid_layers))
                print('  ALGORITHM:', alg)
                print('  ', sbatch_script)
                print('[END SBATCH SCRIPT]')

                os.system('sbatch {}'.format(sbatch_script))
                sleep(10)


if __name__ == '__main__':
    main()
