#!/usr/bin/env python
import click
import os
import numpy as np
from scipy.sparse import random


def random_floats_in_range(a: float, b: float):
    def _loc(n: int):
        return a + (b - a) * np.random.random(n)

    return _loc


# The first row contains 4 integers: the number of rows, the number of columns, the total number of non-zero elements,
# and the maximum number of non-zero elements in each row. The following 3 rows specify values
# (array V in wikipedia's description); column indices (array COL_INDEX); and row offsets (array ROW_INDEX).
# Values may be integers or doubles in format 12.345.
# e.g.
# 6 6 12 2
# 5 3 7 4 8 1 2 9 1 6 2 6
# 1 4 0 5 3 0 2 5 1 5 0 4
# 0 2 4 5 8 10 12
@click.command()
@click.option('--rows', default=10, help='Max number of rows')
@click.option('--cols', default=10, help='Max number of columns')
@click.option('--dens', default=0.3, help='Density of the matrix')
@click.option('--num', default=10, help='Number of matrices')
@click.option('--output', default='output', help='Output directory')
@click.option('--name', default='matrix-{}.csr', help='Name of the matrix')
@click.option('--value', default=100, help='Max value of the matrix')
@click.option('--seed', default=None, help='Random seed')
@click.option('--random-size', default=False, help='Random size of the matrix')
def main(rows: int, cols: int, dens: int, num: int, output: str, name: str, value: float, seed: int | None, random_size: bool):
    if num != 1 and '{}' not in name:
        raise ValueError('Name must contain {}')

    if seed is not None:
        np.random.seed(seed)

    if not os.path.exists(output):
        os.makedirs(output)

    for i in range(num):
        if random_size:
            r = np.random.randint(1, rows + 1)
            c = np.random.randint(1, cols + 1)
        else:
            r = rows
            c = cols

        matrix = random(r, c, density=dens, format='csr', data_rvs=random_floats_in_range(-value, value))
        print(matrix.todense())
        nnz = matrix.nnz

        print('Generating matrix {} with shape {}x{} and {} non-zero elements'.format(i, r, c, nnz))

        with open(os.path.join(output, name.format(i)), 'w') as f:
            f.write('{} {} {} {}\n'.format(r, c, nnz, matrix.getnnz(axis=1).max()))

            def dump_array(array):
                for v in array:
                    f.write('{} '.format(v))
                f.write('\n')

            dump_array(matrix.data)
            dump_array(matrix.indices)
            dump_array(matrix.indptr)


if __name__ == '__main__':
    main()
