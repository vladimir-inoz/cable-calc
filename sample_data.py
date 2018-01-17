from static_cable import cable, cplot, mplot, calculate_forces
from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
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
       'segment_count' : 50}

#рассчитываем конфигурации кабеля
#для нескольких координат ходового конца
#а затем рисуем их
def series(settings, coords):
    for coord in coords:
        settings['x0']=coord[0]
        settings['y0']=coord[1]
        res = calculate_forces(settings)
        if not (res is None):
            cplot(res,False)
    plt.show()

#строим графики зависимости Tx,Ty на
#корневом конце в зависимости от координаты x
def TxTySeries(settings, xr, y):
    Tx = []
    Ty = []
    xl = []
    for x in xr:
        settings['x0']=x
        settings['y0']=y
        res = calculate_forces(inp)
        if not (res is None):
            Tx.append(res['F_winch'].item(0,0))
            Ty.append(res['F_winch'].item(1,0))
            xl.append(x)
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(xl,Tx)
    axarr[0].set_title('Tx')
    axarr[0].grid(visible=True)
    axarr[1].plot(xl,Ty)
    axarr[1].set_title('Ty')
    axarr[1].grid(visible=True)
    plt.grid(True)
    plt.show()
    return [xl,Tx,Ty]
    
    
c = cable(inp)
cplot(c)
