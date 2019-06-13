import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp, log

def V(r, f=1, s=1):
    try:
        return np.array([V(x, f=f, s=s) for x in r])
    except TypeError:
        if r <= s:
            return ((5./18)*pow(f, 3./2.)*(
                -log(r/s)+
                pow((1+sqrt(f)/2), -1)))
        else:
            return ((5./18)*
                    pow(f, 3./2)*
                    pow((1+sqrt(f)/2), -1)*
                    (s/r)*
                    exp(-sqrt(f)*(r-s)/(2*s)))
        if False:
            sig = s
            dist = r
            prefactor = 5./18 * pow(f, 3./2)
            _1_1psf2 = 1./(1.+sqrt(f)/2)
            if dist <= sig:
                return prefactor * (-log(dist/sig) + _1_1psf2)
            else:
                return prefactor * ((sig/dist)*_1_1psf2 *
                        exp(-sqrt(f)/(2*sig)*(dist-sig)))

x = np.linspace(0.00001, 3, 1000)

plt.plot(x, V(x, f=18), '--', color='black')
plt.plot(x, V(x, f=32), '--', color='black')
plt.plot(x, V(x, f=64), '--', color='black')
plt.plot(x, V(x, f=128), '--', color='black')
plt.plot(x, V(x, f=256), '--', color='black')
plt.xlim([0, 3])
plt.ylim([0, 140])
plt.xlabel(r'$r/\sigma$')
plt.ylabel(r'$\beta{}V(r)$')

plt.title('''The pair potential for f = 18, 32,
64, 128, and 256 (from left to right)
as a function of the centerto-center separation r.''')

plt.show()
