#!/usr/bin/env python
import click
import os
import numpy as np
from scipy.sparse import random


def random_floats_in_range(a: float, b: float):
    def _loc(n: int):
        return a + (b - a) * np.random.random(n)

    return _loc


def random_dense_matrix(n: int, m: int, dens: float, values_range: float, verbose: bool):
    matrix = random(n, m, density=dens, format='csr', data_rvs=random_floats_in_range(-values_range, values_range))
    if verbose:
        print(matrix.todense())
    nnz = matrix.nnz

    print('Generating matrix with shape {}x{} and {} non-zero elements'.format(n, m, nnz))

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


@click.command()
@click.option('--n', default=10, help='First dimension of the matrix')
@click.option('--m', default=10, help='Second dimension of the matrix')
@click.option('--k', default=10, help='Third dimension of the matrix')
@click.option('--dens', default=0.3, help='Density of the matrix')
@click.option('--values-range', default=100, help='Max value of the matrix')
@click.option('--seed', default=None, help='Random seed', type=int)
@click.option('--output', default='output', help='Output directory')
@click.option('--verbose', default=False, help='Verbose mode')
@click.option('--name', default='{}.csr', help='Name of the matrix')
def main(n: int, m: int, k: int, dens: int, values_range: float, seed: int | None, output: str, verbose: bool, name: str):
    if seed is not None:
        np.random.seed(seed)

    if '{}' not in name:
        raise ValueError('Name must contain {}')

    A = random_dense_matrix(n, m, dens, values_range, verbose)
    B = random_dense_matrix(m, k, dens, values_range, verbose)
    C = A @ B

    if not os.path.exists(output):
        os.makedirs(output)

    dump_dense_matrix_to_file(A, output, name.format('A'))
    dump_dense_matrix_to_file(B, output, name.format('B'))
    dump_dense_matrix_to_file(C, output, name.format('C'))


if __name__ == '__main__':
    main()
