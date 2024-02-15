import unittest

import numpy as np
from scipy.linalg import block_diag


def sliding_window(collection, window_size):
    if len(collection) < window_size:
        return
    for i in range(len(collection) - window_size + 1):
        yield collection[i : i + window_size]


def is_square(A: np.ndarray) -> bool:
    """
    Checks if object A is a square matrix (NxN)
    @param A: Object to be checked
    @return: True if A is a square matrix, False otherwise
    """
    return A.ndim == 2 and A.shape[0] == A.shape[1]


def inverse_diagonal(D: np.ndarray, tol: float = 1e-6) -> np.ndarray:
    """
    Computes the inverse of a diagonal matrix
    @param D: Diagonal matrix to be inverted
    @return: Inverse of D
    """
    assert is_square(D), "Matrix D must be square"
    return np.diag([1 / d if abs(d) > tol else 0 for d in np.diagonal(D)])


def divide_and_conquer_svd(
    A: np.ndarray, tol: float = 1e-6
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Divide and conquer SVD algorithm
    Decomposes a matrix A into U, S, V^T such that A = U S V^T
    Matrix U is orthogonal, matrix S is diagonal, and matrix V is orthogonal
    @param A: Matrix to be decomposed (m x n)
    @param tol: Tolerance for singular values
    @return: U, S, V^T such that A = U S V^T (m x m, m x n, n x n)
    """
    assert A.dtype == float or A.dtype == np.float
    u, s, vt = _divide_and_conquer_svd(A, tol=tol)
    # assert np.allclose(A, u @ s @ vt), "Decomposition is incorrect"
    return u, s, vt


def _divide_and_conquer_svd(
    A: np.ndarray, tol: float = 1e-6
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Divide and conquer SVD algorithm
    @param A: Matrix to be decomposed (m x n)
    @return: U, S, V^T such that A = U S V^T (m x m, m x n, n x n)
    """
    m, n = A.shape
    if n > m:
        u, s, v = _divide_and_conquer_svd(A.T)
        return v.T, s.T, u.T

    # Step 1: Reduce to bidiagonal form
    u, b, v = _reduce_to_bidiagonal_form(A)
    assert np.allclose(A, u @ b @ v), "Bidiagonalization is incorrect"
    b = b[:n, :n]  # Remove extra rows and columns with zeros

    # Step 2: Compute SVD of bidiagonal matrix
    U, S, V = _divide_and_conquer_svd_bidiagonal(b, tol=tol)
    U = block_diag(U, np.eye(m - n))
    S = np.block([[S], [np.zeros((m - n, n))]])

    return u @ U, S, V @ v


def _row_moving_matrix(k: int, n: int) -> np.ndarray:
    """
    Shifts the k-th row to the top of the matrix
    Moves the first row to the second position, the second row to the third position, and so on
    @param k: Row to be moved
    @param n: Number of rows in the matrix
    @return: Row moving matrix (n x n)
    """
    row_switching_matrix = np.zeros((n, n))
    for i in range(k):
        row_switching_matrix[i, i + 1] = 1
    row_switching_matrix[k, 0] = 1
    for i in range(k + 1, n):
        row_switching_matrix[i, i] = 1

    return row_switching_matrix


def _full_eigen_problem_d_zzt(
    C_dash: np.ndarray, tol: float = 1e-6
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes the eigenvalues and eigenvectors of the matrix D + Z Z^T
    @param C_dash: Matrix D + Z Z^T (n x n)
    @return: Eigenvalues and eigenvectors of D + Z Z^T (S, V.T)
    """
    D = C_dash.copy()
    z = D[0, :].copy()
    D[0, :] = 0
    d = D.diagonal().copy() ** 2

    eigenvalues, eigenvectors = _deflation(d, z, tol=tol)
    assert np.all(eigenvectors is not None)
    YT = np.stack(eigenvectors)
    S = np.diag(eigenvalues)

    return YT, S


def _divide_and_conquer_svd_bidiagonal(
    B: np.ndarray,
    tol: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the SVD of a bidiagonal matrix
    @param B: Bidiagonal matrix (m x n)
    @return: U, S, V^T such that A = U S V^T (m x m, m x n, n x n)
    """
    m, n = B.shape
    if n > m:
        u, s, v = _divide_and_conquer_svd_bidiagonal(B.T, tol=tol)

        return v.T, s, u.T

    if n == 0:
        return np.zeros((m, m)), np.zeros((m, n)), np.zeros((n, n))

    if n == 1:
        return np.eye(m), B, np.eye(n)

    # Decompose B
    k = n // 2
    B_1 = B[:k, : k + 1]
    B_2 = B[k + 1 :, k + 1 :]
    q_k = B[k, k]
    r_k = B[k, k + 1] if k + 1 < n else 0

    U_1, D_1, V_1T = _divide_and_conquer_svd(B_1, tol=tol)
    D_1 = D_1[:, :-1]  # Remove last column of D_1
    U_2, D_2, V_2T = _divide_and_conquer_svd(B_2, tol=tol)

    assert is_square(D_1), "D_1 is not square"
    assert is_square(D_2), "D_2 is not square"

    U_dash = block_diag(U_1, np.eye(1), U_2)
    V_dashT = block_diag(V_1T, V_2T)

    P_k = _row_moving_matrix(k, n)

    lambda_1 = V_1T.T[-1, -1].copy()
    l_1T = V_1T.T[-1, :-1].copy()
    f_2T = V_2T.T[0, :].copy() if k + 1 < n else np.zeros((1, 0))

    C_dash = block_diag(np.zeros((1, 1)), D_1, D_2)

    C_dash[0, :] = np.concatenate(
        [
            np.array([lambda_1 * q_k]),
            q_k * l_1T,
            (r_k * f_2T if k + 1 < n else np.array([])),
        ]
    )

    YT, S = _full_eigen_problem_d_zzt(C_dash, tol=tol)
    S = np.sqrt(S)

    S_inv = inverse_diagonal(S)
    X = C_dash @ YT.T @ S_inv

    return U_dash @ P_k @ X, S, YT @ P_k.T @ V_dashT


def _reduce_to_bidiagonal_form(
    A: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reduces a matrix A to bidiagonal form using Householder reflections
    @param A: Matrix to be reduced to bidiagonal form (m x n)
    @return: U, B, V^T such that A = U B V^T (m x m, m x n, n x n)
    """
    m, n = A.shape
    assert m >= n, "Matrix must have more rows than columns"

    u = np.eye(m)
    v = np.eye(n)
    b = A.copy()

    for k in range(0, n):
        H_k = _householder_from_k(b[:, k], k + 1)
        u = u @ H_k.T
        b = H_k @ b
        H__k = _householder_from_k(b[k, :], k + 2)
        v = H__k @ v
        b = b @ H__k.T

    return u, b, v


def _householder_from_k(x: np.ndarray, k: int) -> np.ndarray:
    """
    Constructs a Householder matrix such that multiplied by the vector `x` it will
    produce a vector with k non-zero elements (at the top of the vector).
    @param x: Vector to construct Householder matrix from
    @param k: Number of non-zero elements in the resulting vector
    @return: Householder matrix H such that H @ x = [*, ..., *, 0, ..., 0]^T where there are k '*'
    """
    upper_left_identity = np.eye(k - 1)
    lower_right_householder = _householder_from(x[k - 1 :])
    return block_diag(upper_left_identity, lower_right_householder)


def _householder_from(x: np.ndarray, k: int = 0) -> np.ndarray:
    """
    Constructs a Householder matrix such that multiplied by the vector `v` it will
    produce a vector with at most on non-zero element (at the top of the vector).
    @param x: Vector to construct Householder matrix from
    @return: Householder matrix H such that H @ x = [*, 0, 0, ..., 0]^T
    """
    assert len(x.shape) == 1, "x must be a vector"
    n = len(x)
    if n == 0:
        return np.eye(n)
    v = x.copy()
    v_norm = np.linalg.norm(v)
    v[k] += np.sign(v[k]) * v_norm
    v_norm = np.linalg.norm(v)
    if n == 0 or not np.nonzero(v_norm):
        return np.eye(n)
    if (v_norm := np.linalg.norm(v)) != 0:
        v /= v_norm
    return np.eye(n) - 2 * np.outer(v, v)


# ==================================================================================================
# ======================================= Finding zeros ============================================
# ==================================================================================================


class LambdaFunction:
    def __init__(self, z: np.ndarray, d: np.ndarray):
        """
        Represents a function of the form $f(x) = 1 + \sum_i^N z_i^2 / (d_i - x)$
        @param z: Parameter vector (N)
        @param d: Parameter vector (N)
        """
        d_z = zip(d, z)
        d_z = sorted(d_z, key=lambda x: x[0])
        self.z = np.array([x[1] for x in d_z])
        self.d = np.array([x[0] for x in d_z])
        # Assert z and d and one-dimensional vectors
        assert len(z.shape) == 1, "z must be a vector"
        assert len(d.shape) == 1, "d must be a vector"
        # Assert z and d have the same length
        assert len(z) == len(d), "z and d must have the same length"
        assert len(z) > 0
        assert len(d) > 0

    def __call__(self, *args, **kwargs):
        """
        Evaluates the function at x
        @param x: Value to evaluate function at
        @return: f(x)
        """
        x = args[0]
        if np.any(np.abs(self.d - x) < 1e-10):
            raise RuntimeError("Zero in divide")
        return 1 + np.sum(self.z**2 / (self.d - x))


def _find_zero_in(
    f: LambdaFunction, a: float, b: float, tol: float = 1e-6, is_last=False
) -> float:
    """
    Finds a zero of the function f in the interval [a, b]
    @param f: Function to find zero of
    @param a: Left endpoint of interval
    @param b: Right endpoint of interval
    @param tol: Tolerance for finding zero
    @return: Zero of f in [a, b]
    """
    assert a < b, "a must be less than b"

    f_a = -1  # We only care about the sign of f(a)
    c = (a + b) / 2
    max_iteration = 1000
    while abs(f_c := f(c)) > tol:
        if (max_iteration := max_iteration - 1) == 0:
            break
        if f_a * f_c < 0:
            b = c
        else:
            a = c
        c = (a + b) / 2
    return c


def _find_zero_after(f: LambdaFunction, a: float, tol: float = 1e-6) -> float:
    """
    Finds a zero of the function f after the point a.
    @param f: Function to find zero of.
    @param a: Point after which to find zero.
    @param tol: Tolerance for finding zero.
    @return: Zero of f after a
    """
    maybe_b = _predict_last_zero(f, tol=tol)
    step = 1
    while f(maybe_b) <= 0:
        maybe_b += step
        step *= 2
    assert a < maybe_b
    return _find_zero_in(f, a, maybe_b, tol, is_last=True)


def _predict_last_zero(f: LambdaFunction, tol: float = 1e-6, model="linear") -> float:
    """
    Predicts the last zero of the function f
    @param f: Function to predict last zero of
    @param tol: Tolerance for finding zero
    @return: Predicted last zero of f
    """

    def constant():
        return f.d[-1] + 10

    def linear():
        return f.d[-1] + np.abs(f.z[-1])

    def square():
        try:
            d_N = f.d[-1]
            d_Nm1 = f.d[-2]
            z_N = f.z[-1]
            z_Nm1 = f.z[-2]
            a = 1
            b = -(d_N + d_Nm1 + z_N + z_Nm1)
            c = d_N * d_Nm1 + d_N * z_Nm1 + d_Nm1 * z_N
            if (discriminant := b**2 - 4 * a * c) >= 0:
                return (-b + np.sqrt(discriminant)) / (2 * a)
        except:
            pass
        return linear()

    models = {
        "const": constant,
        "linear": linear,
        "square": square,
    }

    return models[model]()


def _find_all_zeros(f: LambdaFunction, tol: float = 1e-6) -> np.ndarray:
    """
    Finds all zeros of the function f. There are N zeros.
    N-1 are between d_i and d_{i+1} and the last one is after d_N.
    @param f: Function to find zeros of
    @param tol: Tolerance for finding zeros
    @return: Vector of all zeros of f
    """
    zeros = []
    for d_i, d_i1 in sliding_window(f.d, 2):
        zero_i = _find_zero_in(f, d_i, d_i1, tol)
        zeros.append(zero_i)
    zeros.append(_find_zero_after(f, f.d[-1], tol))
    return np.array(zeros)


def _find_eigenvalue_of_d_z_matrix(
    d: np.ndarray, z: np.ndarray, tol: float = 1e-6
) -> np.ndarray:
    """
    Finds the eigenvalues of the matrix D - Z^T Z
    @param d: Parameter vector (N)
    @param z: Parameter vector (N)
    @param tol: Tolerance for finding zeros
    @return: Eigenvalues of D - Z^T Z
    """
    f = LambdaFunction(z, d)
    zeros = _find_all_zeros(f, tol)
    return zeros


def _find_kth_eigenvector_of_d_z_matrix(
    d: np.ndarray, z: np.ndarray, kth_eigenvalue: float, tol: float = 1e-6
) -> np.ndarray:
    """
    Finds the kth eigenvector of the matrix D + z^T z
    @param d: parameter vector (N)
    @param z: parameter vector (N)
    @param kth_eigenvalue: kth eigenvalue of D + z^T z
    @param tol: Tolerance for finding eigenvectors
    @return: kth eigenvector of D + z^T z
    """
    m_rev = np.array([1 / (d_i - kth_eigenvalue) for d_i in d])
    return m_rev * z


def _find_eigenvectors_of_d_z_matrix(
    d: np.ndarray, z: np.ndarray, eigenvalues: np.ndarray, tol: float = 1e-6
) -> list[np.ndarray]:
    """
    Finds the eigenvectors of the matrix D + z^T z
    @param d: Parameter vector (N)
    @param z: Parameter vector (N)
    @param eigenvalues: Eigenvalues of D + z^T z
    @param tol: Tolerance for finding eigenvectors
    @return: Eigenvectors of D + z^T z
    """
    eigenvectors = []
    for eigenvalue in eigenvalues:
        eigenvector = _find_kth_eigenvector_of_d_z_matrix(d, z, eigenvalue, tol)
        if (v_e_norm := np.linalg.norm(eigenvector)) > tol:
            eigenvector /= v_e_norm
        eigenvectors.append(eigenvector)
    return eigenvectors


def find_eignepairs_of_d_z_matrix(
    d: np.ndarray, z: np.ndarray, tol: float = 1e-6
) -> tuple[np.ndarray, list[np.ndarray]]:
    """
    Finds the eigenvalues and eigenvectors of the matrix D + z^T z
    @param d: Parameter vector (N)
    @param z: Parameter vector (N)
    @param tol: Tolerance for finding eigenvalues and eigenvectors
    @return: Eigenvalues and eigenvectors of D + z^T z
    """
    # d, z = _deflation(d, z, tol)
    eigenvalues = _find_eigenvalue_of_d_z_matrix(d, z, tol)
    eigenvectors = _find_eigenvectors_of_d_z_matrix(d, z, eigenvalues, tol)
    return eigenvalues, eigenvectors


# ==================================================================================================
# ========================================= Deflation ==============================================
# ==================================================================================================


def _has_zeros(z: np.ndarray) -> list[int]:
    """
    Finds the indices of the zeros in z
    @param z: Parameter vector (N)
    @return: Indices of the zeros in z
    """
    return np.count_nonzero(np.array([z])) > 0


def _is_zero(a: float, tol: float = 1e-6) -> bool:
    return np.abs(a) < 1e-3


def _generate_permutation_matrix_from_zeros(
    z: np.ndarray, tol: float = 1e-6
) -> np.ndarray:
    """
    Constructs a permutation matrix that brings zeros in z to the top.
    @param z: Vector with zeros to bring to the top.
    @return: Described permutation matrix.
    """
    assert len(z.shape) == 1
    (n,) = z.shape

    first_free_index = 0
    p = [-1 for _ in range(len(z))]
    for idx, z_i in enumerate(z):
        if _is_zero(z_i, tol=tol):
            p[idx] = first_free_index
            first_free_index += 1
    for idx, z_i in enumerate(z):
        if not _is_zero(z_i, tol=tol):
            p[idx] = first_free_index
            first_free_index += 1
    P = np.zeros((n, n))
    for i, j in enumerate(p):
        P[j, i] = 1
    return P


def e_i(i: int, n: int) -> np.ndarray:
    """
    Constructs the i-th unit vector of dimension n.
    @param i: Index of the unit vector.
    @param n: Dimension of the unit vector.
    @return: i-th unit vector of dimension n.
    """
    assert 0 <= i < n
    e = np.zeros(n)
    e[i] = 1
    return e


def extend_top_vector(v: np.ndarray, n: int) -> np.ndarray:
    """
    Extends the vector v to dimension n by adding zeros to the top.
    @param v: Vector to be extended.
    @param n: Dimension of the extended vector.
    @return: Extended vector.
    """
    assert len(v.shape) == 1
    assert len(v) <= n
    return np.concatenate((np.zeros(n - len(v)), v))


def _deflation(
    d: np.ndarray, z: np.ndarray, tol: float = 1e-6
) -> tuple[np.ndarray, list[np.ndarray]]:
    P = _generate_permutation_matrix_from_zeros(z, tol=tol)

    part_of_eigenvalues = []
    part_of_eigenvactors = []

    dp = P @ np.diag(d) @ P.T
    zp = P @ z

    shorter_d = []
    shorter_z = []

    for idx, z_i in enumerate(zp):
        if _is_zero(z_i):
            part_of_eigenvalues.append(dp[idx, idx])
            part_of_eigenvactors.append(P.T @ e_i(idx, len(zp)))
        else:
            shorter_d.append(dp[idx, idx])
            shorter_z.append(z_i)

    if len(shorter_z) == 0:
        return np.array(part_of_eigenvalues), part_of_eigenvactors

    eigenvalues, eigenvectors = _step_2(
        np.array(shorter_d), np.array(shorter_z), tol=tol
    )

    eigenvectors = [extend_top_vector(v, len(z)) for v in eigenvectors]
    eigenvectors = [P.T @ v for v in eigenvectors]
    return (
        np.array(part_of_eigenvalues + list(eigenvalues)),
        part_of_eigenvactors + eigenvectors,
    )


def _almost_equal(a, b, tol=1e-6) -> bool:
    if a is None or b is None:
        return False
    return np.abs(a - b) < tol


def _step_2(
    d: np.ndarray, z: np.ndarray, tol=1e-6
) -> tuple[np.ndarray, list[np.ndarray]]:
    d_with_indexes = list(enumerate(d))
    d_with_indexes.sort(key=lambda d_p: d_p[1])

    (n,) = d.shape
    P = np.zeros((n, n))
    for i, (j, _) in enumerate(d_with_indexes):
        P[i, j] = 1

    sorted_d = [d[1] for d in d_with_indexes] + [None]
    k = 0
    idx = 0
    res = []
    while sorted_d[idx] is not None:
        if _almost_equal(sorted_d[idx], sorted_d[idx + k], tol=tol):
            k += 1
        else:
            res.append((idx, k - 1))
            idx += k
            k = 0

    k_res = -1  # d_i
    idx_res = -1  # d_{i+k}
    for idx, k in res:
        if k > 0:
            idx_res, k_res = idx, k
            break

    sorted_d = sorted_d[:-1]

    if k_res == -1:
        new_d = np.array(sorted_d)
        new_z = P @ z
        eigenvalues, eigenvactors = find_eignepairs_of_d_z_matrix(new_d, new_z, tol=tol)
        new_eigenvactors = [P.T @ v for v in eigenvactors]
        return eigenvalues, new_eigenvactors

    # Run step one again after modification
    z_perm = P @ z
    d_perm = P @ d
    H_i = _householder_from(z_perm[idx_res : idx_res + k_res + 1], k_res)
    H_m = block_diag(np.eye(idx_res), H_i, np.eye(n - idx_res - k_res - 1))

    M_P = np.zeros((n, n))
    iterator = 0
    a = idx_res
    b = idx_res + k_res
    for i in list(range(a, b + 1)) + list(range(0, a)) + list(range(b + 1, n)):
        M_P[i, iterator] = 1
        iterator += 1

    H = M_P @ H_m

    new_d = np.diagonal(H @ np.diag(d_perm) @ H.T)
    new_z = H @ z_perm
    eigenvalues, eigenvactors = _deflation(new_d, new_z, tol=tol)
    new_eigenvactors = [P.T @ H.T @ v for v in eigenvactors]
    return eigenvalues, new_eigenvactors


# ==================================================================================================
# ========================================= Unit Tests =============================================
# ==================================================================================================


class HouseholderFromTest(unittest.TestCase):
    def _is_single_element(self, H, x):
        check = H @ x
        self.assertTrue(np.allclose(H @ x, check) or np.allclose(H @ x, -check))

    def test_householder_from_1(self):
        x = np.array([1, 2, 3, 4], dtype=float)
        H = _householder_from(x)
        self._is_single_element(H, x)

    def test_householder_from_2(self):
        x = np.array([-1, 2, 3, 4], dtype=float)
        H = _householder_from(x)
        self._is_single_element(H, x)


class HouseholderFromKTest(unittest.TestCase):
    def _is_k_elements(self, H, x, k):
        check = H @ x
        self.assertTrue(np.allclose(H @ x, check) or np.allclose(H @ x, -check))
        n = len(x)
        for i in range(k, n):
            self.assertTrue(np.isclose(check[i], 0))

    def test_householder_from_k_1(self):
        x = np.array([1, 2, 3, 4], dtype=float)
        k = 2
        H = _householder_from_k(x, k)
        self._is_k_elements(H, x, k)

    def test_householder_from_k_2(self):
        x = np.array([-1, 2, 3, 4], dtype=float)
        k = 3
        H = _householder_from_k(x, k)
        self._is_k_elements(H, x, k)


class ReduceToBidiagonalFormTest(unittest.TestCase):
    def test_reduce_to_bidiagonal_form_1(self):
        A = np.array(
            [
                [1, 2, 3],
                [4, 5, 6],
                [7, 8, 9],
                [10, 11, 12],
            ],
            dtype=float,
        )
        u, b, v = _reduce_to_bidiagonal_form(A)
        self.assertTrue(np.allclose(A, u @ b @ v))
        self.assertTrue(np.allclose(u.T @ u, np.eye(4)))
        self.assertTrue(np.allclose(v.T @ v, np.eye(3)))

    def test_reduce_to_bidiagonal_form_2(self):
        A = np.array(
            [
                [1, 2, 3, 4, 5],
                [4, 5, 6, 7, 8],
                [7, 8, 9, 10, 11],
                [10, 11, 12, 13, 14],
                [13, 14, 15, 16, 17],
                [16, 17, 18, 19, 20],
            ],
            dtype=float,
        )
        u, b, v = _reduce_to_bidiagonal_form(A)
        self.assertTrue(np.allclose(A, u @ b @ v))
        self.assertTrue(np.allclose(u.T @ u, np.eye(6)))
        self.assertTrue(np.allclose(v.T @ v, np.eye(5)))

    def test_reduce_to_bidiagonal_form_3(self):
        A = np.fromiter(range(1, 101), dtype=float).reshape((10, 10))
        u, b, v = _reduce_to_bidiagonal_form(A)
        self.assertTrue(np.allclose(A, u @ b @ v))
        self.assertTrue(np.allclose(u.T @ u, np.eye(10)))
        self.assertTrue(np.allclose(v.T @ v, np.eye(10)))
