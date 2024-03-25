# Generates quadrature rules for each polynomial order (p) assuming we want to
# integrate polynomials of order 2 * p + 1 exactly.


import scipy as sp


n_max = 16


def print_array(a):
  for i in range(len(a)):
    print(f'{a[i]:+.17e}', end='')
    if i < len(a) - 1:
      print(',', end=' ')


for n in range(1, n_max + 1):

  x, w = sp.special.roots_legendre(n)

  print(f'const usize  gauss_legendre_n{n}_n                      = {n};')
  print(f'const double gauss_legendre_n{n}_x[gauss_legendre_n{n}_n] = {{', end='')
  print_array(x)
  print('};')
  print(f'const double gauss_legendre_n{n}_w[gauss_legendre_n{n}_n] = {{', end='')
  print_array(w)
  print('};')

  print()
