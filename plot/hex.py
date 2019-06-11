import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def do(fn, label, **kw):
    offsets = np.loadtxt(fn)
    m = offsets.mean(axis=0) / 3
    print(offsets.shape)
    m = np.roll(m, 5)
    m[10:] += 1
    def sg(x, s, m, d):
        return s * np.arctan(np.exp(m*x+d))
    cells = np.arange(19)
    (S, M, D), _cov = curve_fit(sg, cells, m, [0.5, 1, -10])
    # print('COV')
    # print(_cov)
    # print('scale', S)
    # print('m gamma', M)
    # print('delta', D)
    x = np.linspace(0, 19, 1000)
    plt.plot(cells, m, 'o', label=label, **kw)
    plt.plot(x, sg(x, S, M, D), '-', **kw)

do('saved/hex_offsets0', 'beta=0.000', color='#000000')
do('saved/hex_offsets1', 'beta=0.001', color='#222222')
do('saved/hex_offsets2', 'beta=0.002', color='#444444')
do('saved/hex_offsets3', 'beta=0.003', color='#666666')
do('saved/hex_offsets4', 'beta=0.004', color='#888888')
plt.legend()
plt.title('Hexagonal vacancy - projection on p3 with sine-gordon fit')
plt.xlabel('site')
plt.ylabel('<u>/a')
plt.show()
