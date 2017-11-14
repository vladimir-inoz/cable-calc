import static_cable


#аппарат на глубине 300 м
#под носителем
#считаем конфигурации кабеля на разных скоростях движения
res = static_cable.calculate_forces(0.1,5000,5000,5000,-100,-100,0,10.0,50)
static_cable.print_cable(res['joints'])
