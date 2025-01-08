#!/usr/bin/env python3
# Keep file in Python 3.6 compatible

import logging
import os
import tempfile
import numpy as np
from scipy.sparse import random
import subprocess
import shutil
import scipy.sparse as sp
from pathlib import Path

SEED = 69  # 42
VALUE_RANGE = 100
LAYERS = 5
NUM_PROC = 4
EXEC_COMMAND = ['mpiexec', '-n', str(NUM_PROC * LAYERS), '--hostfile', 'hostfile']
EXEC_PATH = './build/matmul'
DUMP_PATH = './test-data-dump'
ALGORITHM = '3D'
EPS = 1e-6
CHECK_ONLY_G = False

np.random.seed(SEED)

logFormatter = logging.Formatter(fmt=' %(name)s :: %(levelname)-8s :: %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)


def random_floats_in_range(a: float, b: float):
    def _loc(n: int):
        return a + (b - a) * np.random.random(n)

    return _loc


def random_dense_matrix(n: int, m: int, dens: float, values_range: float):
    matrix = random(n, m, density=dens, format='csr', data_rvs=random_floats_in_range(-values_range, values_range))
    return matrix


def dump_dense_matrix_to_file(matrix, output: str, filename: str):
    r, c = matrix.shape
    nnz = matrix.nnz
    with open(os.path.join(output, filename), 'w') as f:
        f.write('{} {} {} {}\n'.format(r, c, nnz, matrix.getnnz(axis=1).max()))

        def dump_array(array):
            for v in array:
                f.write('{} '.format(v))
            f.write('\n')

        dump_array(matrix.data)
        dump_array(matrix.indices)
        dump_array(matrix.indptr)


def generate_test_data(test_dir, test_name, size, dens, g=None):
    A = random_dense_matrix(size, size, dens, VALUE_RANGE)
    B = random_dense_matrix(size, size, dens, VALUE_RANGE)
    C = A @ B

    dump_dense_matrix_to_file(A, test_dir, '{}_A.csr'.format(test_name))
    dump_dense_matrix_to_file(B, test_dir, '{}_B.csr'.format(test_name))

    if CHECK_ONLY_G:
        assert g is not None
        g_count = np.sum(C.data > g)
        with open(os.path.join(test_dir, '{}_g_count.txt'.format(test_name)), 'w') as f:
            f.write('{}\n'.format(g_count))
    else:
        dump_dense_matrix_to_file(C, test_dir, '{}_C.csr'.format(test_name))


def load_matrix(file_path):
    with open(file_path, 'r') as f:
        num_rows, num_cols, nnz, _ = map(int, f.readline().split())
        values = np.array(list(map(float, f.readline().split())), dtype=float)
        cols = np.array(list(map(int, f.readline().split())), dtype=int)
        rows = np.array(list(map(int, f.readline().split())), dtype=int)

        assert len(values) == nnz
        assert len(cols) == nnz
        assert len(rows) == num_rows + 1

        assert rows[0] == 0
        assert rows[-1] == nnz

        assert all(rows[i] <= rows[i + 1] for i in range(num_rows))
        assert all(cols[i] < num_cols for i in range(nnz))

        matrix = sp.csr_matrix((values, cols, rows), shape=(num_rows, num_cols))
        return matrix


def parse_dense_matrix(stdout):
    lines = stdout.split('\n')
    n, m = map(int, lines[0].split())
    matrix = np.zeros((n, m))
    for i in range(n):
        matrix[i] = list(map(float, lines[i + 1].split()))
    return matrix


def matrix_error(expected, obtained):
    return np.max(np.abs(expected - obtained))


def run_code(test_dir, test_name, g=None):
    proc = subprocess.run(
        EXEC_COMMAND + [
            EXEC_PATH, '-a', '{}/{}_A.csr'.format(test_dir, test_name), '-b',
            '{}/{}_B.csr'.format(test_dir, test_name),
            '-t', ALGORITHM, '-l', str(LAYERS)
        ] + (['-g', str(g)] if CHECK_ONLY_G else ['-v']), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    obtained_output = proc.stdout
    with open('{}/out.txt'.format(DUMP_PATH), 'w') as f:
        f.write(obtained_output)
    with open('{}/err.txt'.format(DUMP_PATH), 'w') as f:
        f.write(proc.stderr)
    if CHECK_ONLY_G:
        expected_g_count = int(open('{}/{}_g_count.txt'.format(test_dir, test_name)).read())
        obtained_g_count = int(obtained_output)
        logger.info('Expected g count: {}, obtained g count: {}'.format(expected_g_count, obtained_g_count))
        return expected_g_count == obtained_g_count
    else:
        expected_matrix = load_matrix('{}/{}_C.csr'.format(test_dir, test_name))
        obtained_matrix = parse_dense_matrix(obtained_output)
        error = matrix_error(expected_matrix.toarray(), obtained_matrix)
        logger.info('Error: {}'.format(error))
        return error < EPS


def run_tests(test_dir):
    cases = [
        (20, 'small-{}', 10, 0.3, 11.1),
        (10, 'medium-{}', 100, 0.3, 12.3),
        (2, 'large-{}', 1000, 0.3, 120.4),
        (2, 'large-large-{}', 10000, 0.0001, 140.8),
    ]
    Path(DUMP_PATH).mkdir(parents=True, exist_ok=True)
    for iterations, name, size, dens, g in cases:
        for i in range(iterations):
            test_name = name.format(i)
            logger.info('Running test: {}'.format(test_name))
            generate_test_data(test_dir, test_name, size, dens, g)
            result = run_code(test_dir, test_name, g)
            if not result:
                Path(DUMP_PATH).mkdir(parents=True, exist_ok=True)
                logger.error('Test failed [ERROR]: {}'.format(test_name))
                shutil.copy('{}/{}_A.csr'.format(test_dir, test_name), '{}/{}_A.csr'.format(DUMP_PATH, test_name))
                shutil.copy('{}/{}_B.csr'.format(test_dir, test_name), '{}/{}_B.csr'.format(DUMP_PATH, test_name))
                if CHECK_ONLY_G:
                    shutil.copy('{}/{}_g_count.txt'.format(test_dir, test_name),
                                '{}/{}_g_count.txt'.format(DUMP_PATH, test_name))
                else:
                    shutil.copy('{}/{}_C.csr'.format(test_dir, test_name), '{}/{}_C.csr'.format(DUMP_PATH, test_name))
                return
            else:
                logger.info('Test passed [OK]: {}'.format(test_name))


def cleanup(test_dir):
    shutil.rmtree(test_dir, ignore_errors=True)
    logger.info('Removed temporary directory: {}'.format(test_dir))


def setup_testing_environment():
    test_dir = tempfile.mkdtemp()
    logger.info('Created temporary directory: {}'.format(test_dir))
    return test_dir


def main():
    test_dir = setup_testing_environment()
    run_tests(test_dir)
    cleanup(test_dir)


if __name__ == '__main__':
    main()
