import math
import matplotlib.pyplot as plt
import numpy as np


def rates(h, e):
  return (np.log(e[1:]) - np.log(e[:-1])) / (np.log(h[1:]) - np.log(h[:-1]))


l = 2. * math.pi
h = l / np.array([8, 16, 32, 64, 128, 256])


# p = 1

e1 = np.array([
  9.47535021e-02,
  4.11502035e-02,
  1.11351840e-02,
  2.84704554e-03,
  7.15915222e-04,
  1.79241766e-04
])

# p = 2

e2 = np.array([
  3.26611394e-02,
  2.51445409e-03,
  2.80913336e-04,
  3.39369109e-05,
  4.20445024e-06,
  5.24372849e-07
])

# p = 3

e3 = np.array([
  1.02025762e-03,
  2.60323670e-04,
  1.88150839e-05,
  1.19630362e-06,
  7.45475321e-08,
  4.65169139e-09
])

# p = 4

e4 = np.array([
  1.20879781e-02,
  1.65423599e-03,
  8.65550483e-06,
  7.45133234e-08,
  2.14717465e-09,
  6.58298022e-11
])


rates1 = rates(h, e1)
rates2 = rates(h, e2)
rates3 = rates(h, e3)
rates4 = rates(h, e4)

plt.loglog(h, e1, 'D-k', label=f'p = 1, rate = {rates1[-1]:.1f}')
plt.loglog(h, e2, '+-k', label=f'p = 2, rate = {rates2[-1]:.1f}')
plt.loglog(h, e3, 'x-k', label=f'p = 3, rate = {rates3[-1]:.1f}')
plt.loglog(h, e4, '*-k', label=f'p = 4, rate = {rates4[-1]:.1f}')
plt.xticks(h, [r'$1 / 2^3$', r'$1 / 2^4$', r'$1 / 2^5$', r'$1 / 2^6$', r'$1 / 2^7$', r'$1 / 2^8$'])
plt.minorticks_off()
plt.xlabel(r'$h / 2 \pi$')
plt.ylabel(r'$|| u - u_{exact} ||_{2}$')
plt.legend()

plt.savefig('convergence_uxx.pdf', bbox_inches='tight')

plt.show()
