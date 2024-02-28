import sympy as sy


# Constructs a polynomial up to order 15 based on fixed coefficeints.

def polynomial(order, x):

  coeff = [0.49330289, 0.90463143, 0.18816676, 0.30646078, 0.18025477, 0.77162931, 0.03942296, 0.33686216, 0.26149018, 0.49742156, 0.47202936, 0.32424461, 0.25503326, 0.26378048, 0.50724305, 0.94544638]

  result = 0.
  for i in range(order + 1):
    result += coeff[i] * (x ** float(i))

  return result


x = sy.Symbol('x')

for o in range(15 + 1):
  integral = sy.integrate(polynomial(o, x), (x, -1., +1.))
  print(f'{integral:+.17e}')
