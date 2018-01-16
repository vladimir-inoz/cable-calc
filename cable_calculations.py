import static_cable
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

#максимальный упор движков 150 Н
Pmax = 150
#максимальные упоры X,Y
Fxmax = 259.8
Fymax = 300

#аппарат на глубине 300 м
#под носителем
#считаем конфигурации кабеля на разных скоростях движения
inp = dict(vx = 6, \
                vy = 0, \
                vz = 0, \
                Fx = -83.72, \
                Fy = -22104.8, \
                Fz = 0, \
                cable_len = 550/1000, \
                segment_count = 1000, \
                Cxa = 1.25, \
                Cza = 1.05, \
                Sa = 0.181, \
                Dk = 0.016, \
                cable_mass = 0
           )

res = static_cable.cable(inp)
static_cable.print_cable(res)

plt.show()
