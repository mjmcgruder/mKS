import math
import matplotlib.pyplot as plt
import numpy as np


def rates(h, e):
  return (np.log(e[1:]) - np.log(e[:-1])) / (np.log(h[1:]) - np.log(h[:-1]))


l = 2. * math.pi
h = l / np.array([8, 16, 32, 64, 128])


# p = 2

e2 = np.array([
  2.89803071e-01,
  3.73783835e-02,
  1.00276612e-02,
  2.55815799e-03,
  6.42903961e-04
])

# p = 3

e3 = np.array([
  1.07293305e-02,
  8.93218560e-04,
  4.65481393e-05,
  2.73963422e-06,
  1.68519391e-07
])

# p = 4

e4 = np.array([
  1.03814351e-03,
  4.16751989e-05,
  1.52718051e-06,
  5.02681749e-08,
  1.59423091e-09
])


rates2 = rates(h, e2)
rates3 = rates(h, e3)
rates4 = rates(h, e4)

plt.loglog(h, e2, '+-k', label=f'p = 2, rate = {rates2[-1]:.1f}')
plt.loglog(h, e3, 'x-k', label=f'p = 3, rate = {rates3[-1]:.1f}')
plt.loglog(h, e4, '*-k', label=f'p = 4, rate = {rates4[-1]:.1f}')
plt.xticks(h, [r'$1 / 2^3$', r'$1 / 2^4$', r'$1 / 2^5$', r'$1 / 2^6$', r'$1 / 2^7$'])
plt.minorticks_off()
plt.xlabel(r'$h / 2 \pi$')
plt.ylabel(r'$|| u - u_{exact} ||_{2}$')
plt.legend()

plt.savefig('convergence_uxxxx.pdf', bbox_inches='tight')

plt.show()
