import numpy as np


def print_c_array(a):
  print('{')
  for r in range(a.shape[0]):
    for c in range(a.shape[1]):
      print(f'{a[r, c]:+.17e},', end='')
      if c < a.shape[1] - 1: print(end=' ')
    print()
  print('};\n')


rng = np.random.default_rng(42)


A = rng.uniform(size=(10, 10))
B = rng.uniform(size=(10, 20))
C = A @ B
D = np.linalg.inv(A)

print('A:')
print_c_array(A)
print('B:')
print_c_array(B)
print('A * B:')
print_c_array(C)
print('A ** -1:')
print_c_array(D)
