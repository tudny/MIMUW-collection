#!/usr/bin/env python
import click
import scipy.sparse as sp
import numpy as np
import scipy


@click.command()
@click.argument('matrix', type=click.Path(exists=True))
def main(matrix: str):
    with open(matrix, 'r') as f:
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

        print('=' * 200)
        print('Matrix with shape {}x{} and {} non-zero elements'.format(num_rows, num_cols, nnz))
        print('-- Values:', values)
        print('-- Columns:', cols)
        print('-- Rows:', rows)

        matrix = sp.csr_matrix((values, cols, rows), shape=(num_rows, num_cols))
        print('=' * 200)
        print('-- Matrix:')
        np.set_printoptions(linewidth=200, precision=3)
        print(np.array(matrix.todense()))


if __name__ == '__main__':
    main()
