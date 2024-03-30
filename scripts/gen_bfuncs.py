import numpy as np


def lagrange1d(i, p, x):

  nodes = np.cos(((2. * np.arange(1, p + 2) - 1.) / (2. * (p + 1.))) * np.pi)

  xi   = nodes[i]
  eval = 1.
  for j in range(p + 1):
    if i != j:
      xj    = nodes[j]
      eval *= (x - xj) / (xi - xj)
  return eval


if __name__ == '__main__':

  import sympy as sy
  import sympy.codegen.rewriting as rw

  expand_pow = rw.create_expand_pow_optimization(128)

  p_max = 15

  x = sy.Symbol('x')

  # generate set of symbolic basis functions

  sym_evals = []

  for p in range(0, p_max + 1):
    sym_evals.append([])
    for bf in range(p + 1):
      sym_evals[-1].append(sy.simplify(lagrange1d(bf, p, x)))

  # generate C code for basis fucntions

  print('BASIS:\n')
  for p in range(0, p_max + 1):
    evals = sym_evals[p]
    print(f'case {p}:')
    for i in range(len(evals)):
      print(f'b[{i}] =', end=' ')
      print(sy.ccode(expand_pow(evals[i])), end=';\n')
    print('break;\n')

  # generate derivatives

  for d in range(1, 4):
    print(f'DERIVATIVE {d}\n')
    for p in range(0, p_max + 1):
      evals = sym_evals[p]
      print(f'case {p}:')
      for i in range(len(evals)):
        print(f'b[{i}] =', end=' ')
        print(sy.ccode(expand_pow(sy.simplify(sy.diff(evals[i], x, d)))), end=';\n')
      print('break;\n')
