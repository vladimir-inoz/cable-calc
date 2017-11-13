import numpy as np
import math
import pylab
import matplotlib.pyplot as plt

#расчет направляющих косинусов звена по силе
def calcOrientCos(F):
    #считаем натяжение в шарнире
    T = np.linalg.norm(F)
    C = F / T
    return C

#расчет параметров кабеля

#диаметр кабеля (метры)
Dk = 0.01

#нормальный и касательный гидродинамические коэффициенты кабеля
Cnorm = 1.0
Ctau = 1.0

#число сегментов
nSegments = 50

#длина сегмента кабеля
l = 1

def cable(Vx, Vz, Px, Py, Pz):
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
    V = np.matrix([[Vx],[0],[Vz]])

    #матрица В
    B = np.matrix([[math.cos(psi), -math.sin(psi)],
    [math.sin(psi), math.cos(psi)]])

    #
    Va = -B.transpose() * np.matrix([[Vx],[Vz]])
    Vxa = Va.item((0,0))
    Vza = Va.item((1,0))

    #силы, действующие на аппарат
    Fxa = ro * Cxa * Sa * Vxa * abs(Vxa) / 2
    Fza = ro * Cza * Sa * Vza * abs(Vza) / 2
    print(Fxa)

    tR1 = -B * np.matrix([[Fxa],[Fza]])
    R1 = -np.matrix([[tR1.item(0,0)],[0],[tR1.item(1,0)]])

    #силы, развиваемые ДРК ТПА
    Px1 = Px
    Py1 = Py
    Pz1 = Pz
    P1 = np.matrix([[Px1],[Py1],[Pz1]])

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
        
        Rnx = Cnorm * l * Dk * ro * \
              (Vx * math.sqrt(1-Cphix ** 2)) ** 2 / 2
        Rnz = Cnorm * l * Dk * ro * \
              (Vx * math.sqrt(1-Cphiz ** 2)) ** 2 / 2
        Rtaux = Ctau * Cnorm * l * Dk * math.pi * \
                np.sign(Cphix)*(Vx*Cphix)**2/2
        Rtauz = Ctau * Cnorm * l * Dk * math.pi * \
                np.sign(Cphiz)*(Vz*Cphiz)**2/2

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
        Fsum[0]+=Force.item(0,0) ** 2
        Fsum[1]+=Force.item(1,0) ** 2
        Fsum[2]+=Force.item(2,0) ** 2
    Fsums = np.matrix([[math.sqrt(Fsum[0])],
                       [math.sqrt(Fsum[1])],
                       [math.sqrt(Fsum[2])]])

    #рассчитываем координаты шарниров
    #в системе координат НПА
    #первый шарнир в нулях
    r = [np.matrix([[0],[0],[0]])]
    for i in range(1, nSegments):
        rnew = r[i-1] - l * C[i]
        r.append(rnew)

    #рассчитываем координаты шарниров
    #в системе координат носителя
    rnos=[]
    for i in range(0, nSegments):
        rnos_new = r[-1] - r[i]
        rnos.append(rnos_new)
    
    return rnos

#рисование кабеля на графике
def print_cable(joints):
    xs = []
    ys = []
    for i in joints:
        xs.append(i.item(0,0))
        ys.append(i.item(2,0))
    plt.plot(xs,ys)

#нахождение упоров движителей НПА
#для определенного положения (X,Y,Z)
#делается путем минимизации функционала
#def calculate_forces(Vx,Fxmax,Fymax,Fzmax,x0,y0,z0):
    #минимизируемый функционал
#    def force_func(x):
#        return (x[0]-x0)**2 + (x[1]-y0)**2 + (x[2]-z0)**2
        
    

def functional(Vx,Fxmin,Fxmax,Fzmin,Fzmax,x0,y0,z0):
    for Fx in range(Fxmin,Fxmax):
        for Fz in range(Fzmin,Fzmax):
            [xs,ys] = cable(Vx,0,Fx,Fz)
            x=xs[-1]
            y=ys[-1]
            z=0
            J = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
            print(J)

axes = plt.gca()
axes.set_xlim([0,100])
axes.set_ylim([0,100])

for i in range(0,5):
    joints = cable(i/10, 0, 100, 100, 0)
    print_cable(joints)

plt.show()
