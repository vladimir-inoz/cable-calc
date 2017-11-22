import static_cable
import matplotlib.pyplot as plt

#аппарат на глубине 300 м
#под носителем
#считаем конфигурации кабеля на разных скоростях движения
inp = dict(vx = 0.52, \
                vy = 0, \
                vz = 0, \
                Fx = 223.1, \
                Fy = -20.5, \
                Fz = 0, \
                segment_len = 5.0, \
                segment_count = 40, \
           )
res = static_cable.cable(inp)
static_cable.print_cable(res['joints'])

plt.show()
