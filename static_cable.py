import numpy as np
import math
import pylab
import matplotlib.pyplot as plt
from scipy import optimize

#расчет направляющих косинусов звена по силе
def calcOrientCos(F):
    #считаем натяжение в шарнире
    T = np.linalg.norm(F)
    C = F / T
    return C

#расчет параметров кабеля

#диаметр кабеля (метры)
Dk = 0.05

#нормальный и касательный гидродинамические коэффициенты кабеля
Cnorm = 1.0
Ctau = 1.0

#число сегментов
nSegments = 10

def cable(settings):
    #результат
    result = dict()

    #длина сегмента
    len = settings['segment_len']
    
    #курсовой угол корабля (радианы)
    psi = 0

    #коэффициенты гидродинамические для ТПА
    #принимаем ТПА как кубик
    Cxa = 1.05
    Cza = 1.05

    #плотность воды (кг\м3)
    ro = 1000.0

    #характерная площадь НПА
    #примем, что аппарат 1x1x1м
    #характерная площадь - 2/3 от объема
    Sa = 0.666

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

    Ft= np.matrix([[0],[0],[0]])

    #пересчитываем силы для всех сегментов
    for i in range(1, nSegments):
        Cphix = C[i-1].item(0,0)
        Cphiy = C[i-1].item(1,0)
        Cphiz = C[i-1].item(2,0)
        Fprev = F[i-1]
        Gprev = np.matrix([[0],[0],[0]])
        
        Rnx = Cnorm * len * Dk * ro * \
              (vx * math.sqrt(1-Cphix ** 2)) ** 2 / 2
        Rnz = Cnorm * len * Dk * ro * \
              (vx * math.sqrt(1-Cphiz ** 2)) ** 2 / 2
        Rtaux = Ctau * Cnorm * len * Dk * math.pi * \
                np.sign(Cphix)*(vx*Cphix)**2/2
        Rtauz = Ctau * Cnorm * len * Dk * math.pi * \
                np.sign(Cphiz)*(vz*Cphiz)**2/2

        Rprev = np.matrix([[Rnx],[Rtaux],[Rnz],[Rtauz]])

        Cyz = -math.sqrt(Cphix ** 2 + Cphiy ** 2 + Cphiz ** 2)
        Cyx = math.sqrt(Cphix ** 2 + Cphiy ** 2 + Cphiz ** 2)
        Cxy = -Cphix * Cphiy / Cyz
        Cxz = -Cphix * Cphiz / Cyz
        Czx = Cphiz * Cphix / Cyx
        Czy = Cphiz * Cphiy / Cyx

        Aprev = np.matrix([
            [Cyz, -Cphix, Czx, -Cphix],
            [Cxy, -Cphiy, Czy, -Cphiy],
            [Cxz, -Cphiz, -Cyx, -Cphiz]
            ])
        
        #результирующая сила на i-м шарнире
        Fi = Fprev + Aprev * Rprev + Gprev + Ft
        #ориентация звена
        Ci = calcOrientCos(Fi)
        #добавляем значения в массив
        F.append(Fi)
        C.append(Ci)

    #считаем силу на лебедке корабля
    Fsum=[0,0,0]
    for force in F:
        Fsum[0]+=force.item(0,0) ** 2
        Fsum[1]+=force.item(1,0) ** 2
        Fsum[2]+=force.item(2,0) ** 2
    Fsums = np.matrix([[math.sqrt(Fsum[0])],
                       [math.sqrt(Fsum[1])],
                       [math.sqrt(Fsum[2])]])
    result['fsum'] = Fsums

    #рассчитываем координаты шарниров
    #в системе координат НПА
    #первый шарнир в нулях
    r = [np.matrix([[0],[0],[0]])]
    for i in range(1, nSegments):
        rnew = r[i-1] + len * C[i]
        r.append(rnew)

    result['joints'] = r

    #рассчитываем координаты шарниров
    #в системе координат носителя
    rnos=[]
    for point in r:
        rnos.append(r[-1] - point)
    result['joints-nos'] = rnos

    return result

#рисование кабеля на графике
def print_cable(joints):
    xs = []
    ys = []
    for i in joints:
        xs.append(i.item(0,0))
        ys.append(i.item(1,0))
    plt.plot(xs,ys,'bo')
    plt.plot(xs,ys,'k')

#нахождение упоров движителей НПА
#для определенного положения (X,Y,Z)
#делается путем минимизации функционала
def calculate_forces(Vx,Fxmax,Fymax,Fzmax,x0,y0,z0):
    #минимизируемый функционал
    #входные параметры - упоры движителей аппарата
    #выходной - квадрат дистанции от заданной точки
    def force_func(f):
        joints = cable(Vx, 0, f[0], f[1], f[2])
        print(joints)
        return (joints[-1].item(0,0)-x0)**2 + \
               (joints[-1].item(1,0)-y0)**2 + \
               (joints[-1].item(2,0)-z0)**2
    #осуществляем минимизацию
    #методом BFGS
    #http://www.scipy-lectures.org/advanced/mathematical_optimization/
    #у нас нет якобиана и градиента
    max_force = math.sqrt(Fxmax**2 + Fymax**2 + Fzmax**2)
    res = optimize.minimize( \
        force_func, [0,0,0], method = "BFGS")
    print(res)
    #рисуем получившийся кабель

axes = plt.gca()
axes.set_xlim([-100,50])
axes.set_ylim([-100,5])

inp = dict(vx=3.0,vy=0,vz=0,Fx=200,Fy=-200,Fz=0,segment_len=5.0)
for x in range(0,12):
    inp['vx'] = x*0.25
    print_cable(cable(inp)['joints'])

#vx = 1.0
#joints = cable(vx, 0, 100, 0, 0)
#print_cable(joints)

plt.show()
