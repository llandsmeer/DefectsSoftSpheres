import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Axis:
    def __init__(self, spec):
        parts = spec.split()
        self.label = parts[0]
        self.state = []
        self.offset = []
        for i in range(1, len(parts), 2):
            self.state.append(float(parts[i]))
            self.offset.append(float(parts[i+1]))
        self.state = np.array(self.state)
        self.offset = np.array(self.offset)
        assert len(self.state) == len(self.offset)

    def __repr__(self):
        return f'<Axis {self.label} {len(self.state)}>'


m = 0
m2 = 0
m3 = 0
N = 0

for line in open('./sim/test/sc_offsets'):
    specs = line.split(' : ')
    axes = [Axis(spec) for spec in specs]
    best = axes[-1]
    m = m + best.offset
    m2 += axes[-2].offset
    m3 += axes[-3].offset
    N += 1

m = m / N
m2 /= N
m3 /= N

a = 3
print('a =', a)

m[len(m)//2:] += a
m2[len(m)//2:] += a
m3[len(m)//2:] += a

cells = np.arange(len(m))
cells[len(m)//2:] += 1
plt.plot(cells, -m3/a, 'o', color='#999999')
plt.plot(cells, -m2/a, 'o', color='#666666')
plt.plot(cells, -m/a, 'o', color='black')

def sg(x, s, m, d):
    return s * np.arctan(np.exp(m*x+d))

(S, M, D), _cov = curve_fit(sg, cells, m, [0.5, 1, -10])
x = np.linspace(0, len(m)+1, 1000)
plt.plot(x, -sg(x, S, M, D)/a, '--', color='black')

plt.ylabel('$<u>/a$')
plt.xlabel('Site')
plt.xticks(cells)
plt.title('SC vacancy, mean displacement and sine-Gordon fit \n$\\rho{} \\sigma^3=5.2$, $k_BT/\\epsilon$=0.001')

plt.show()

