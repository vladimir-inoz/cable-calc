from static_cable import cable, cplot, mplot
inp = {'vx': 0.5, \
       'vy': 0, \
       'vz': 0, \
       'Dk': 0.005, \
       'Cxa': 0.7, \
       'Cza': 0, \
       'Sa': 0.5574, \
       'Fx': 120, \
       'Fy': -120, \
       'Fz': 0, \
       'cable_mass': 0.0255, \
       'x0': 0, \
       'y0': -300, \
       'z0': 0, \
       'Fxmax': 250, \
       'Fymax': 350, \
       'Fzmax': 0, \
       'cable_len' : 350, \
       'segment_count' : 100}
cplot(cable(inp))
