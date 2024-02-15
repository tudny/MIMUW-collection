
import numpy as np

m, n = 3, 5

# 1) Zrobić macierz losową (m x n)
A = np.random.rand(m, n)

print(f'{A=}')

# 2) Zrobić losowy wektor b (m)
b = np.random.rand(m)

print(f'{b=}')

# 3) Rozwiązać układ ||b - Ax||_2 -> min!
x = np.linalg.lstsq(A, b)

print(f'{x=}')


# Sprawdzenie
print(f'{np.linalg.norm(b - A @ x[0])=}')
