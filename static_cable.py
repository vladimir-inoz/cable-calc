import numpy as np
import math
from math import pi
import pylab
import matplotlib.pyplot as plt
from math import pi
from scipy import optimize
import matplotlib.cm as cm

#расчет направляющих косинусов звена по силе
def calcOrientCos(F):
    #считаем натяжение в шарнире
    T = np.linalg.norm(F)
    C = F / T
    return C

#расчет параметров кабеля

#нормальный и касательный гидродинамические коэффициенты кабеля
Cnorm = 1.2
Ctau = 0.02
#Cnorm = 0.2
#Ctau = 0.2

def cable(settings):
    #результат
    result = dict()

    #длина сегмента
    len = settings['segment_len']
    #число сегментов
    nseg = settings['segment_count']
    #диаметр кабеля
    Dk = settings['Dk']
    
    #курсовой угол корабля (радианы)
    psi = 0

    #коэффициенты гидродинамические для ТПА
    #принимаем ТПА как кубик
    Cxa = settings['Cxa']
    Cza = settings['Cza']

    #плотность воды (кг\м3)
    ro = 1025.0

    #характерная площадь НПА
    #примем, что аппарат 1x1x1м
    #характерная площадь - 2/3 от объема
    Sa = settings['Sa']

    #Скорость движения ТПС отностиельно воды
    vx = settings['vx']
    vy = settings['vy']
    vz = settings['vz']
    V = np.matrix([[vx],[vy],[vz]])

    #матрица В
    B = np.matrix([[math.cos(psi), -math.sin(psi)],
    [math.sin(psi), math.cos(psi)]])

    #
    Va = -B.transpose() * np.matrix([[vx],[vz]])
    Vxa = Va.item((0,0))
    Vza = Va.item((1,0))

    #силы, действующие на аппарат
    Fxa = ro * Cxa * Sa * Vxa * abs(Vxa) / 2
    Fza = ro * Cza * Sa * Vza * abs(Vza) / 2

    tR1 = -B * np.matrix([[Fxa],[Fza]])
    R1 = -np.matrix([[tR1.item(0,0)],[0],[tR1.item(1,0)]])

    #силы, развиваемые ДРК ТПА
    result['input'] = settings
    P1 = np.matrix([[settings['Fx']], \
                    [settings['Fy']], \
                    [settings['Fz']]])

    #остаточная плавучесть ТПА
    G1 = np.matrix([[0],[0],[0]])

    #массив сил в шарнирах
    F=[]

    #массив направляющих косинусов звеньев
    C=[]

    #Сила, приложенная к 1-му шарниру
    F.append(R1 + P1 + G1)
    C.append(calcOrientCos(F[0]))

    #массив гидродинамических сил на звеньях
    R=[]

    G=[]

    #плавучесть сегмент кабеля
    cable_mass = settings['cable_mass']

    Ft= np.matrix([[0],[0],[0]])

    #пересчитываем силы для всех сегментов
    for i in range(1, nseg):
        Cphix = C[i-1].item(0,0)
        Cphiy = C[i-1].item(1,0)
        Cphiz = C[i-1].item(2,0)

        Fprev = F[i-1]
        Gk = 9.81 * (cable_mass - ro * pi * Dk * Dk / 4)
        Gprev = np.matrix([[0],[-Gk * len],[0]])

        #гидродинамические силы на звене кабеля
        Rnx = Cnorm * len * Dk * ro * \
              (vx * math.sqrt(1-Cphix ** 2)) ** 2 / 2
        Rnz = Cnorm * len * Dk * ro * \
              (vz * math.sqrt(1-Cphiz ** 2)) ** 2 / 2
        Rtaux = Ctau * Cnorm * len * Dk * math.pi * \
                np.sign(Cphix)*(vx*Cphix)**2/2
        Rtauz = Ctau * Cnorm * len * Dk * math.pi * \
                np.sign(Cphiz)*(vz*Cphiz)**2/2

        Rprev = np.matrix([[Rnx],[Rtaux],[Rnz],[Rtauz]])

        Cyz = -math.sqrt(Cphiy ** 2 + Cphiz ** 2)
        Cyx = math.sqrt(Cphix ** 2 + Cphiy ** 2)
        
        if (math.fabs(Cyz) > 1e-3):
            Cxy = -Cphix * Cphiy / Cyz
            Cxz = -Cphix * Cphiz / Cyz
        else:
            Cxy = 0
            Cxz = 0
        if (math.fabs(Cyx) > 1e-3):
            Czx = Cphiz * Cphix / Cyx
            Czy = Cphiz * Cphiy / Cyx
        else:
            Czx = 0
            Cxz = 0

        Aprev = np.matrix([
            [Cyz, -Cphix, Czx, -Cphix],
            [Cxy, -Cphiy, Czy, -Cphiy],
            [Cxz, -Cphiz, -Cyx, -Cphiz]
            ])

        #гидродинамические силы на i-1 звено
        Fgydro = Aprev * Rprev
        
        #результирующая сила на i-м шарнире
        Fi = Fprev + Fgydro + Gprev + Ft
        #ориентация звена
        Ci = calcOrientCos(Fi)
        #добавляем значения в массив
        F.append(Fi)
        C.append(Ci)

    result['forces'] = F
    result['orient'] = C

    #считаем силу на лебедке корабля
    Fsum=[0,0,0]
    for force in F:
        Fsum[0]+=force.item(0,0) ** 2
        Fsum[1]+=force.item(1,0) ** 2
        Fsum[2]+=force.item(2,0) ** 2
    Fsums = np.matrix([[math.sqrt(Fsum[0])],
                       [math.sqrt(Fsum[1])],
                       [math.sqrt(Fsum[2])]])
    result['f_winch'] = Fsums

    #рассчитываем координаты шарниров
    #в системе координат НПА
    #первый шарнир в нулях
    r = [np.matrix([[0],[0],[0]])]
    for i in range(1, nseg):
        rnew = r[i-1] - len * C[i]
        r.append(rnew)

    result['joints'] = r

    #рассчитываем координаты шарниров
    #в системе координат носителя
    rnos=[]
    for point in r:
        rnosnew = point - r[-1]
        rnos.append(rnosnew)
    
    result['joints_nos'] = rnos

    return result

#рисование кабеля на графике
def print_cable(res):
    xs = []
    ys = []
    zs = []
    for i in res['joints_nos']:
        xs.append(i.item(0,0))
        ys.append(i.item(1,0))
        zs.append(i.item(2,0))
        
    #plt.plot(xs,ys,'bo')
    plt.plot(xs,ys,'k')
    

    #plt.subplot(212)

    #axes = plt.gca()
    #axes.set_xlim([-100,50])
    #axes.set_ylim([-100,5])
    
    #plt.plot(xs,zs,'bo')
    #plt.plot(xs,zs,'k')
    
#нахождение упоров движителей НПА
#для определенного положения (X,Y,Z)
#делается путем минимизации функционала
def calculate_forces(vx,Fxmax,Fymax,Fzmax,x0,y0,z0,segment_len = 5.0,\
                     segment_count = 20):
    #минимизируемый функционал
    #входные параметры - упоры движителей аппарата
    #выходной - квадрат дистанции от заданной точки
    def force_func(f):
        inp = dict(vx = vx, \
                   vy = 0, \
                   vz = 0, \
                   Fx = f[0], \
                   Fy = f[1], \
                   Fz = f[2], \
                   segment_len = segment_len, \
                   segment_count = segment_count)
        joints = cable(inp)['joints']
        return (joints[-1].item(0,0)-x0)**2 + \
               (joints[-1].item(1,0)-y0)**2 + \
               (joints[-1].item(2,0)-z0)**2
    #осуществляем минимизацию
    #методом Бройдена-Флетчера-Гольдфарба-Шанно #http://www.scipy-lectures.org/advanced/mathematical_optimization/
    #у нас нет якобиана и градиента
    max_force = math.sqrt(Fxmax**2 + Fymax**2 + Fzmax**2)
    res = optimize.minimize( \
        force_func, [0,0,0], method = "BFGS")
    if (res.success == True):
        #рисуем получившийся кабель
        inp = dict(vx = vx, \
                   vy = 0, \
                   vz = 0, \
                   Fx = res.x[0], \
                   Fy = res.x[1], \
                   Fz = res.x[2], \
                   segment_len = 5.0, \
                   segment_count = 20)
        #возвращаем рассчитанную конфигурацию кабеля
        return cable(inp)
    else:
        return None

#нахождение рабочей зоны аппарата
def calculate_workzone(vx,Fxmax,Fymax,Fzmax,segment_len,segment_count):
    xs = []
    ys = []
    
    def getxy(fx,fy):
        inp = dict(vx = vx, \
                vy = 0, \
                vz = 0, \
                Fx = fx, \
                Fy = fy, \
                Fz = 0, \
                segment_len = segment_len, \
                segment_count = segment_count)
        npa_x = cable(inp)['joints'][-1].item(0,0)
        npa_y = cable(inp)['joints'][-1].item(1,0)
        xs.append(npa_x)
        ys.append(npa_y)
    
    n = 100
    #длина кабеля
    cable_len = segment_len * segment_count;
    #делаем все возможные упоры движителей и считаем координаты ТПА
    Fmax = math.sqrt(Fxmax**2 + Fymax**2 + Fzmax**2)
    for x in np.linspace(-Fxmax,Fxmax,100):
        for y in np.linspace(-Fymax,Fymax,100):
            getxy(x,y)
        
    plt.plot(xs,ys,'go')
