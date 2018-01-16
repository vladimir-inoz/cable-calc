#проверка экспериментальных данных от статей ИПМТ
import static_cable
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

inp = dict(vx = 0.5, \
                vy = 0, \
                vz = 0, \
                Fx = 0, \
                Fy = 5.6, \
                Fz = 0, \
                cable_len = 6.8, \
                segment_count = 10, \
                Cxa = 0.2, \
                Cza = 1.05, \
                Sa = 0.0185**0.6, \
                Dk = 0.007, \
                cable_mass = 0.0414
           )

inp2 = dict(vx = 0.5, \
                vy = 0, \
                vz = 0, \
                Fx = 0, \
                Fy = 20, \
                Fz = 0, \
                cable_len = 6.8, \
                segment_count = 100, \
                Cxa = 0.2, \
                Cza = 0.2, \
                Sa = 0.0699, \
                Dk = 0.007, \
                cable_mass = 0.3
           )


res = static_cable.cable(inp2)
inp2['vx'] = 1.0
res1 = static_cable.cable(inp2)
inp2['vx'] = 1.5
res2 = static_cable.cable(inp2)
inp2['vx'] = 2.0
res3 = static_cable.cable(inp2)

static_cable.print_cable(res)
static_cable.print_cable(res1)
static_cable.print_cable(res2)
static_cable.print_cable(res3)

plt.show()
