from math import *
from matplotlib import pyplot as plt
import numpy as np

m = 3920
S = 19.5
Vlof = 46
Czmax_lisse = 1.52
Czmax_TO = 1.93
Czmax_Land = 2.2
k = 0.05
Cx0 = 0.025
g = 9.81


def rho(h):
    T = 288.15 - 0.0065 * h
    P = 101325 * (1 - 22.557 * 10**-6 * h) ** 5.226
    return P / (287.05 * T)


def ft_to_m_or_m_to_ft(l, choice='fttom' or 'mttoft'):
    if choice == 'mtoft':
        l *= 3.281
    elif choice == 'fttom':
        l /= 3.281
    return l


def Cz_Cx_f(ro, V):
    Cz = 2 * m * g / (ro * S * V**2)
    Cx = 0.025 + 0.05 * Cz ** 2
    f = Cz / Cx
    return Cz, Cx, f


def poussée_max(delta_m, ro):
    return delta_m * (ro/1.225)**0.6 * 12980


def poussée_requise(ro, V, Vz):
    Cx = Cz_Cx_f(ro=0.7712, V=V)[1]
    Fn = 0.5 * ro * S * V**2 * Cx + m * g * Vz/V
    return Fn


def poussée_nécessaire(ro, V):
    Cx = Cz_Cx_f(ro=ro, V=V)[1]
    return 0.5 * ro * S * V**2 * Cx


def V_decrochage(ro, Czmax=Czmax_lisse):
    return sqrt(2 * m * g / (ro * S * Czmax))


def Vmax(ro):
    Fmax = poussée_max(1, ro)
    CzVmax = 1/(2 * k) * (Fmax / (m * g) - sqrt((Fmax / (m * g))**2 - 4 * k * Cx0))
    return V_decrochage(ro, Czmax=CzVmax)


def alt_plafond_sustentation():
    pf = 2 * m * g / (1.4 * S * 0.34**2 * Czmax_lisse)
    return -((pf/101325)**(1/5.226) - 1) / (22.557 * 10**-6)


def distance_decollage(r, ro, Czmax):
    Fm = poussée_max(delta_m=1, ro=ro)
    Lroulage = Vlof**2 / (Fm/m - r*g)
    Vs_T0 = V_decrochage(ro, Czmax=Czmax)
    V2 = 1.2 * Vs_T0
    f = Cz_Cx_f(ro, V2)[2]
    Lenvol = (10.5 + (V2**2 - Vlof**2) / (2 * g)) / (Fm / (m*g) - 1 / f)
    return Lroulage + Lenvol


def pentewVz(ro, V, delta_m):
    Cx = Cz_Cx_f(ro=ro, V=V)[1]
    pente = asin((poussée_max(delta_m=delta_m, ro=ro) - 0.5 * ro * S * V**2 * Cx) / (m * g))
    Vz = sin(pente) * V
    return pente * (180/pi), Vz


def enveloppe_vol():
    x = []
    Vm = []
    Vd = []
    for i in range(0, int(alt_plafond_sustentation())):
        x.append(i)
        Vm.append(Vmax(ro=rho(i)))
        Vd.append(V_decrochage(ro=rho(i)))
    plt.plot(Vm, x, label="Vmax", color='purple')
    plt.plot(Vd, x, label="Vdécrochage", color='blue')
    plt.axhline(y=alt_plafond_sustentation(), label="Plafond de sustentation="+str(int(alt_plafond_sustentation()))+"m"
                , color='red')
    plt.title("Enveloppe de vol avec plafond de sustentation")
    plt.xlabel("Vitesse en m/s")
    plt.ylabel("Altitude en m")
    plt.legend(loc='lower right')
    plt.grid()
    plt.show()


def pousséereq_altitude(V=165):
    h = []
    Tu = []
    for i in range(0, int(alt_plafond_sustentation()) + 1):
        h.append(i)
        Tu.append(poussée_nécessaire(ro=rho(i), V=V))
    plt.plot(h, Tu, color='purple')
    plt.title("Poussée requise en fonction de l'altitude pour une vitesse de 165m/s")
    plt.xlabel("Altitude en m")
    plt.ylabel("Poussée requise en N")
    plt.grid()
    plt.show()


def endurance_en_h(ro, Cs, V):
    Cz, Cx = Cz_Cx_f(ro, V)[0:-1]
    return (Cz / (g * Cs/3600 * Cx) * log(m/(m - 1110))) / 3600


def endurance_vitesse():
    v = []
    e = []
    for i in range(int(V_decrochage(ro=rho(10000))), int(Vmax(ro=rho(10000))) + 1):
        v.append(i)
        e.append(endurance_en_h(ro=rho(10000), Cs=0.093, V=i))
    plt.plot(v, e, color='blue')
    plt.plot(Vmax(ro=rho(10000)), endurance_en_h(ro=rho(10000), Cs=0.093, V=Vmax(ro=rho(10000))),
             'ko', label='Vmax à 10000m = ' + str(int(Vmax(ro=rho(10000))))+'m/s')
    plt.plot(V_decrochage(ro=rho(10000)), endurance_en_h(ro=rho(10000), Cs=0.093, V=V_decrochage(ro=rho(10000))),
             'ro', label='VDécrochage à 10000m = ' + str(int(V_decrochage(ro=rho(10000))))+'m/s')
    plt.title("Endurance à une altitude de 10000m" + '\n' + "selon la vitesse avec Cs = 0.093kg/N.h")
    plt.xlabel("Vitesse en m/s")
    plt.ylabel("Endurance en h")
    plt.legend()
    plt.grid()
    plt.show()


def rayon_en_km(ro, Cs, V):
    Cz, Cx = Cz_Cx_f(ro, V)[0:-1]
    return 2 / (g * Cs/3600 * Cx) * sqrt(2 * Cz / (ro * S)) * (sqrt(m * g) - sqrt((m - 1110) * g)) * 10**-3


def rayon_vitesse():
    v = []
    r = []
    for i in range(int(V_decrochage(ro=rho(10000))), int(Vmax(ro=rho(10000))) + 1):
        v.append(i)
        r.append(rayon_en_km(ro=rho(10000), Cs=0.093, V=i))
    plt.plot(v, r, color='blue')
    plt.plot(Vmax(ro=rho(10000)), rayon_en_km(ro=rho(10000), Cs=0.093, V=Vmax(ro=rho(10000))),
             'ko', label='Vmax à 10000m = ' + str(int(Vmax(ro=rho(10000)))) + 'm/s')
    plt.plot(V_decrochage(ro=rho(10000)), rayon_en_km(ro=rho(10000), Cs=0.093, V=V_decrochage(ro=rho(10000))),
             'ro', label='VDécrochage à 10000m = ' + str(int(V_decrochage(ro=rho(10000)))) + 'm/s')
    plt.title("Rayon d'action en km à une altitude de 10000m" + '\n' + "selon la vitesse avec Cs = 0.093kg/N.h")
    plt.xlabel("Vitesse en m/s")
    plt.ylabel("Rayon d'action en km")
    plt.legend()
    plt.grid()
    plt.show()


def puissance_virage(ro, V, phi):
    n = 1/cos(phi)
    veq = sqrt(n)*V
    Cx = Cz_Cx_f(ro, V)[1]
    return 0.5 * ro * S * veq**2 * Cx


def poussée_inclinaison():
    phi = []
    F = []
    h = ft_to_m_or_m_to_ft(4500, choice='fttom')
    for i in range(0, 85):
        phi.append(i)
        F.append(puissance_virage(ro=rho(h), V=81, phi=radians(i))*10**-3)
    plt.plot(phi, F, color='blue')
    plt.title("Poussée requise en fonction de l'angle de virage")
    plt.xlabel("Inclinaison en degrés")
    plt.ylabel("Poussée requise en kN")
    plt.grid()
    plt.show()


def Vdecrovirage(ro, Czmax, phi):
    vs = V_decrochage(ro, Czmax)
    n = 1/cos(phi)
    return sqrt(n) * vs


def vdecrovirage_phi():
    phi = []
    vTO = []
    vLand = []
    vlisse = []
    h = ft_to_m_or_m_to_ft(4500, choice='fttom')
    for i in range(0, 85):
        phi.append(i)
        vTO.append(Vdecrovirage(ro=rho(h), Czmax=Czmax_TO, phi=radians(i)))
        vLand.append(Vdecrovirage(ro=rho(h), Czmax=Czmax_Land, phi=radians(i)))
        vlisse.append(Vdecrovirage(ro=rho(h), Czmax=Czmax_lisse, phi=radians(i)))
    plt.plot(phi, vTO, color='blue', label='Czmax_TO =' + str(Czmax_TO))
    plt.plot(phi, vLand, color='green', label='Czmax_Land =' + str(Czmax_Land))
    plt.plot(phi, vlisse, color='red', label='Czmax_lisse =' + str(Czmax_lisse))
    plt.title("Vitesse de décrochage en fonction de l'inclinaison")
    plt.xlabel("Inclinaison en degrés")
    plt.ylabel("Vitesse de décrochage en m/s")
    plt.grid()
    plt.legend()
    plt.show()


def domainedevol():
    r = rho(10000)
    Vslisse = V_decrochage(ro=r, Czmax=Czmax_lisse)
    VsTO = V_decrochage(ro=r, Czmax=Czmax_TO)
    nmax, nmin = 2.5, -1
    V, Vm, Vs = 165, Vmax(ro=r), Vslisse
    Vitesse = np.linspace(0, Vm, 550)
    n_pos, n_neg, n_volets = [], [], []

    for i in Vitesse:
        m1 = i**2 / Vs**2
        m2 = -i**2 / Vs**2
        m3 = 0.5 * r * S * i**2 * Czmax_lisse / (m * g)

        if m1 < nmax:
            n_pos.append(i**2 / Vs**2)
        else:
            n_pos.append(nmax)
        if m2 > nmin:
            n_neg.append(-i**2 / Vs ** 2)
        elif i >= V:
            n_neg.append((-(Vm - i) / (Vm - V)))
        else:
            n_neg.append(nmin)
        if m3 < nmax:
            n_volets.append(0.5 * r * S * i**2 * Czmax_lisse / (m * g))
        else:
            if i < 112:
                n_volets.append(2)
            else:
                if m1 < nmax:
                    n_volets.append(i**2 / Vs**2)
                else:
                    n_volets.append(nmax)
    n_pos[-1], n_volets[-1] = 0, 0
    plt.plot(Vitesse, n_volets, color="red", label="Domaine de vol avec volets sortis")
    plt.plot(Vitesse, n_neg, color="pink", label="Domaine de vol avec n<0")
    plt.plot(Vitesse, n_pos, color="purple", label="Domaine de vol avec volets rentrés")
    plt.xlabel("Vitesse en m/s")
    plt.ylabel("Facteur de charge n")
    plt.title("Domaine de vol")
    plt.grid()
    plt.legend()
    plt.show()


print("Partie 1 :")
print("Distance de décollage :", int(distance_decollage(r=0.1, ro=1.1588, Czmax=Czmax_TO)), "m", '\n'*2)
print(rho(ft_to_m_or_m_to_ft(1900, choice='fttom')))
print("Partie 2 :")
print("Poussée requise en montée FL350 :", int(poussée_requise(ro=0.7708, V=88, Vz=15.24)), "N")
print("Poussée max en montée FL350 :", int(poussée_max(delta_m=1, ro=0.7708)), "N")
print("Pente & vitesse de montée :", pentewVz(ro=0.7708, V=88, delta_m=1))
print("Poussée necessaire avec une panne moteur :", int(poussée_nécessaire(ro=0.7523, V=65)), "N")
print("Pente & vitesse de montée avec une panne moteur", pentewVz(ro=0.7523, V=65, delta_m=0.5))
enveloppe_vol()
pousséereq_altitude()
domainedevol()
endurance_vitesse()
rayon_vitesse()
poussée_inclinaison()
vdecrovirage_phi()
