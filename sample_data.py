from static_cable import cable, cplot, mplot, calculate_forces
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import math

inp = {'vx': 0.5, \
       'vy': 0, \
       'vz': 0, \
       'Dk': 0.005, \
       'Cxa': 0.7, \
       'Cza': 0, \
       'NPA_buoyancy': -20.0, \
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
       'cable_len': 350, \
       'segment_count': 50}


# рассчитываем конфигурации кабеля
# для нескольких координат ходового конца
# а затем рисуем их
def series(settings, coords):
    for coord in coords:
        settings['x0'] = coord[0]
        settings['y0'] = coord[1]
        res = calculate_forces(settings)
        if not (res is None):
            cplot(res, False)
    plt.show()


# математическая модель потребления аппаратом
# мощности
# аппарат с двумя вертикальными и двумя горизонтальными
# движителями
def NPAPower(settings):
    Kn = 0.9
    # КПД винтомотора
    efficiency = 0.4
    Nvert = 0.707 * Kn * (math.fabs(settings['Fy']) ** 1.5)
    Nhor = math.fabs(settings['Fx']) * math.fabs(settings['vx']) / efficiency
    Nsum = Nvert + Nhor
    return Nsum


# строим графики зависимости Fx,Fy на корневом
# конце в зависимости от координаты x
#
def fx_fy_series(settings, x, y):
    Fx = np.zeros(x.size)
    Fy = np.zeros(x.size)
    Fsum = np.zeros(x.size)
    for i in range(0, x.size):
        settings['x0'] = x[i]
        settings['y0'] = y
        res = calculate_forces(inp)
        if not (res is None):
            Fx_res = res['input']['Fx']
            Fy_res = res['input']['Fy']
            Fsum_res = math.fabs(Fx_res) + math.fabs(Fy_res)
            Fx[i] = Fx_res
            Fy[i] = Fy_res
            Fsum[i] = Fsum_res
    f, axarr = plt.subplots(3, sharex=True)
    plt.title('Упор движителей')
    axarr[0].plot(x, Fx)
    axarr[0].grid(visible=True)
    axarr[0].set_xlabel('x, м')
    axarr[0].set_ylabel('Fx, Н')
    axarr[1].plot(x, Fy)
    axarr[1].grid(visible=True)
    axarr[1].set_xlabel('x, м')
    axarr[1].set_ylabel('Fy, Н')
    axarr[2].plot(x, Fsum)
    axarr[2].grid(visible=True)
    axarr[2].set_xlabel('x, м')
    axarr[2].set_ylabel('|Fx|+|Fy|, Н')

    plt.show()
    return [Fx, Fy, Fsum]


# строим графики зависимости Tx,Ty на
# ходовом конце в зависимости от координаты x
# и для разных скоростей движения
def tx_ty_npa_series(settings, x, y, v):
    Tx = np.ndarray(shape=(v.size, x.size), dtype=float)
    Ty = np.ndarray(shape=(v.size, x.size), dtype=float)
    Tsum = np.ndarray(shape=(v.size, x.size), dtype=float)
    for i in range(0, v.size):
        for j in range(0, x.size):
            settings['x0'] = x[j]
            settings['y0'] = y
            settings['vx'] = v[i]
            res = calculate_forces(inp)
            if not (res is None):
                Tx_res = res['F_NPA'].item(0, 0)
                Ty_res = res['F_NPA'].item(1, 0)
                Tsum_res = math.sqrt(Tx_res * Tx_res + Ty_res * Ty_res)
                Tx[i][j] = Tx_res
                Ty[i][j] = Ty_res
                Tsum[i][j] = Tsum_res
    f, axarr = plt.subplots(3, sharex=True)
    plt.title('Натяжение кабеля на ходовом конце')
    for i in range(0, v.size):
        axarr[0].plot(x, Tx[i])
        axarr[0].text(x[-2], Tx[i][-2] + 1, 'v=%2.1f' % v[i])
    axarr[0].grid(visible=True)
    axarr[0].set_xlabel('x, м')
    axarr[0].set_ylabel('Tx, Н')
    for i in range(0, v.size):
        axarr[1].plot(x, Ty[i])
        axarr[1].text(x[-2], Ty[i][-2] + 1, 'v=%2.1f' % v[i])
    axarr[1].grid(visible=True)
    axarr[1].set_xlabel('x, м')
    axarr[1].set_ylabel('Ty, Н')
    for i in range(0, v.size):
        axarr[2].plot(x, Tsum[i])
        axarr[2].text(x[-2], Tsum[i][-2] + 1, 'v=%2.1f' % v[i])
    axarr[2].grid(visible=True)
    axarr[2].set_xlabel('x, м')
    axarr[2].set_ylabel('T, Н')

    plt.show()
    return [Tx, Ty, Tsum]

# строим графики зависимости Tx,Ty на
# ходовом конце в зависимости от координаты x
def tx_ty_anpa_series(settings, x, y):
    Tx = np.zeros(x.size)
    Ty = np.zeros(x.size)
    Tsum = np.zeros(x.size)
    for i in range(0, x.size):
        settings['x0'] = x[i]
        settings['y0'] = y
        res = calculate_forces(inp)
        if not (res is None):
            Tx_res = res['F_winch'].item(0, 0)
            Ty_res = res['F_winch'].item(1, 0)
            Tsum_res = math.sqrt(Tx_res * Tx_res + Ty_res * Ty_res)
            Tx[i] = Tx_res
            Ty[i] = Ty_res
            Tsum[i] = Tsum_res
    f, axarr = plt.subplots(3, sharex=True)
    plt.title('Натяжение кабеля на корневом конце')
    axarr[0].plot(x, Tx)
    axarr[0].grid(visible=True)
    axarr[0].set_xlabel('x, м')
    axarr[0].set_ylabel('Tx, Н')
    axarr[1].plot(x, Ty)
    axarr[1].grid(visible=True)
    axarr[1].set_xlabel('x, м')
    axarr[1].set_ylabel('Ty, Н')
    axarr[2].plot(x, Tsum)
    axarr[2].grid(visible=True)
    axarr[2].set_xlabel('x, м')
    axarr[2].set_ylabel('T, Н')

    plt.show()
    return [Tx, Ty, Tsum]


# 2D график распределения натяжения кабеля на
# ходовом конце
def plot_npa_2d(settings):
    nx = 200
    ny = 200
    u, v = np.meshgrid(np.linspace(-settings['Fxmax'], settings['Fxmax'], nx),
                       np.linspace(-settings['Fymax'], settings['Fymax'], ny),
                       sparse=False, indexing='ij')
    x_t = []
    y_t = []
    z_t = []
    for i in range(nx):
        for j in range(ny):
            fx = u[i, j]
            fy = v[i, j]
            settings['Fx'] = fx
            settings['Fy'] = fy
            res = cable(settings)
            x_npa = res['joints_nos'][0].item(0, 0)
            y_npa = res['joints_nos'][0].item(1, 0)
            x_t.append(x_npa)
            y_t.append(y_npa)
            Tx_res = res['F_NPA'].item(0, 0)
            Ty_res = res['F_NPA'].item(1, 0)
            Tsum_res = math.sqrt(Tx_res * Tx_res + Ty_res * Ty_res)
            z_t.append(Tsum_res)
    x = np.array(x_t)
    y = np.array(y_t)
    z = np.array(z_t)
    # define grid.
    xmin = -150
    xmax = 150
    ymin = -400
    ymax = 0
    xi = np.linspace(xmin, xmax, 100)
    yi = np.linspace(ymin, ymax, 100)
    # grid the data.
    zi = griddata(x, y, z, xi, yi, interp='linear')
    # contour the gridded data, plotting dots at the nonuniform data points.
    CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    CS = plt.contourf(xi, yi, zi, 15,
                      vmax=abs(zi).max(), vmin=-abs(zi).max())
    plt.colorbar()  # draw colorbar
    # plot data points.
    plt.scatter(x, y, marker='o', s=5, zorder=10)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.title('griddata test (%d points)' % 400)
    plt.show()

Tx, Ty, Tsum = tx_ty_npa_series(inp,np.linspace(-100,90,10),-300, np.array((0.2,0.3,0.4,0.5)))
np.savetxt('data.txt',[Tx,Ty,Tsum])
