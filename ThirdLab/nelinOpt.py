from Matrica import Matrica, JedinicnaMatrica
from functools import reduce
import math
import numpy as np
import pandas as pd



# - UNIMODALNI -

# Ulazne velicine:
# - tocka: pocetna tocka pretrazivanja
# - h: pomak pretrazivanja
# - f: ciljna funkcija

def unimodalni(x0, func, h_value = 1.0):
    
    h = Matrica(1, 1, [[h_value]])
    tocka = x0.copy()
    left = tocka - h 
    right = tocka + h
    m = tocka
    step = 1
    
    # print("unimodalni point: ", x0)
    fm = func.vrijednost(tocka)
    # print("f(point): ", fm)
    fl = func.vrijednost(left)
    fr = func.vrijednost(right)
    
    if fm < fr and fm < fl: 
        return left, right
    elif fm > fr:
        while True:
            left = m
            m = right
            fm = fr
            step *= 2
            right = tocka + (h * step)
            fr = func.vrijednost(right)
            if (fm <= fr): 
                break
    else: 
        while True:
            right = m
            m = left
            fm = fl
            step *= 2
            left = tocka - (h * step)
            fl = func.vrijednost(left)
            if fm <= fl: break
    return left, right



# - ZLATNI REZ -

# Ulazne velicine:
# - a, b: pocetne granice unimodalnog intervala
# - e: preciznost

def zlatni_rez(a, b, func, epsilon = 1e-6, max_iter = 20000):
    
    k = 0.5 * (math.sqrt(5) - 1)
        
    c = b - k * (b - a)
    d = a + k * (b - a)
    fc = func.vrijednost(c)
    fd = func.vrijednost(d)
    
    br_iter = 0
    while((b[0][0] - a[0][0]) > epsilon and br_iter <= max_iter):
        if fc < fd:
            b = d 
            d = c
            c = b - k * (b - a)
            fd = fc
            fc = func.vrijednost(c)
        else:
            a = c
            c = d
            d = a + k * (b - a)
            fc = fd
            fd = func.vrijednost(d)
        br_iter += 1
    return (a + b) / 2



# - KOORDINATNO PRETRAZIVANJE -

def koord_pretr(x0, func, epsilon = 1e-6, ispis = False):
    n = x0.br_stup
    xs = x0.copy()
    x = x0.copy()

    if ispis:
        print("n: ", n)
        print("x: ", x)

    flag = True
    while flag:
        xs = x.copy()
        for i in range(n):
            xi = Matrica(1, 1, [[x[0][i]]])

            #jedinicni vektor
            e = Matrica(1, n)
            e[0][i] = 1.0

            if ispis:
                print("iteracija ", i, "-----------------------------------")
                print("x: ", x)
                print("e: ", e)

            func_i = Funkcija(lambda l: func.vrijednost(x + l*e))
            gamma_min = minimum(xi, func_i, e, ispis)

            if ispis:
                print("l_min: ", gamma_min)
                print("x: ", x)
            x = x + gamma_min * e
            if ispis: 
                print("minimum u koordinati {0}: {1}".format(i, x))
                # print("x = ", x)
                # print("\n-----------------\n")
            
        # print("x: ", x, "xs: ", xs)

        for i in range(n):
            if not abs(xs[0][i] - x[0][i]) > epsilon: 
                # print("-----------")
                # print("STOP")
                # print("-----------")
                flag = False
                break
    return x

def minimum(x0, func, e, ispis = False):
    if ispis: print("point: ", x0)
    a, b = unimodalni(x0, func)
    if ispis: print("a: {0}, b: {1}".format(a, b))
    gamma_min = zlatni_rez(a, b, func)[0][0]
    if ispis: print("golden ratio: ", gamma_min)
    return gamma_min



# - SIMPLEX -

def simplex(x0, func, pomak = 1, ispis = False, alfa = 1.0, beta = 0.5, gama = 2.0, sigma = 0.5, epsilon=1e-6):
    
    #tocke simplexa
    pocetna_tocka = x0.copy()
#     print("pocetna tocka simplexa: ", x0)
    X = [pocetna_tocka]
    for i in range(x0.br_stup):
        nova_tocka = x0.copy()
        nova_tocka[0][i] += pomak
        X.append(nova_tocka)

    #ispis simpleks tocaka
    if ispis: 
        print("X: ", end="\t")
        for i in range(len(X)):
            print(X[i], end=" ")
        print("\n")
    
    # - h, l, centroid -
    while True:
        h = index_max(X, func, ispis = ispis)
        l = index_min(X, func, ispis = ispis)
        xc = centroid(X, h)
        fc = func.vrijednost(xc)
            
        xr = refleksija(xc, X[h], alfa)
        fr = func.vrijednost(xr)
        
        #ispis rjesenja svakog pojednog koraka
        if ispis:
            print("h: ", h, " |l: ", l, " |xc: ", xc)
            print("fc: ", fc)
            print("xr, fr: ", xr, fr)
        
        fl = func.vrijednost(X[l])
#         print("fl: ", fl)
        if fr < fl:
            xe = ekspanzija(xr, xc, gama)
#             print("fe: ", fe)
            
            if func.vrijednost(xe) < fl:
                X[h] = xe.copy()
            else:
                X[h] = xr.copy()

        else:
            # provjera
            flag = True
            for j in range(len(X)): 
                if j != h:
                    if func.vrijednost(xr) < func.vrijednost(X[j]):
                        flag = False
                        break
            if flag:
                if func.vrijednost(xr) < func.vrijednost(X[h]):
                    X[h] = xr.copy()

                xk = kontrakcija(X[h], xc, beta)

                if func.vrijednost(xk) < func.vrijednost(X[h]):
                    X[h] = xk.copy()
                else:
                    pomak_prema_xl(X, l, sigma)
            else:
                X[h] = xr.copy()
                
        # Uvjet zaustavljanja
        s = 0
        for i in range(len(X)):
            s += (func.vrijednost(X[i]) - func.vrijednost(xc))**2
        if ispis:
            print("X[l], X[h], xc:")
            print(X[l], X[h], xc)
        
        if math.sqrt(s / float(len(X) - 1)) <= epsilon:
            return xc
            
        
# - MINIMUM I MAKSIMUM FJE CILJA        
def index_min(X, func, ispis = False):
    if ispis:
        print("ulaz u min")
        print("X: ", end = " ")
        for i in range(len(X)):
            print(X[i], end = " ")
        print("\n")
    i_min = 0
    for i in range(len(X)):
        fi = func.vrijednost(X[i])
        fi_min =  func.vrijednost(X[i_min])
        if ispis: print("{0} < {1}".format(fi, fi_min), end = " ")
        if fi < fi_min: 
            if ispis: print("DA")
            i_min = i
        else: 
            if ispis: print("NE")
#         if func.vrijednost(X[i]) < func.vrijednost(X[i_min]): i_min = i
    return i_min
 
def index_max(X, func, ispis = False):
    if ispis:
        print("ulaz u max")
        print("X: ", end=" ")
        for i in range(len(X)):
            print(X[i], end = " ")
        print("\n")
    i_max = 0
    for i in range(len(X)):
        fi = func.vrijednost(X[i])
        fi_max =  func.vrijednost(X[i_max])
        if ispis: print("{0} > {1}".format(fi, fi_max), end = " ")
        if fi > fi_max: 
            if ispis: print("DA")
            i_max = i
        else:
            if ispis: print("NE")
    return i_max

# - CENTROID -
def centroid(X, h):
#     print("X: ", X)
    tocke = X.copy()
#     print("X[h]: ", tocke[h])
    del tocke[h]
#     print("tocke: ", tocke)
    
    if len(tocke) == 1:
        centroid = tocke[0].elementi
#         print("1: ", centroid)
    else: 
        red = reduce(myadd, tocke)
        centroid = (1.0 / len(tocke)) * red
#         print("2: ", centroid)
    return Matrica(1, X[0].br_stup, [centroid[0]])

def myadd(a, b):
    return Matrica(a.br_red, a.br_stup, a.elementi + b.elementi)

# - REFLEKSIJA, EKSPANZIJA, KONTRAKCIJA -
def refleksija(xc, xh, alfa = 1.0):
    return ((1.0 + alfa) * xc) - (alfa * xh)

def ekspanzija(xc, xr, gama = 2.0):
    return ((1.0 - gama) * xc) + (gama * xr)

def kontrakcija(xh, xc, beta = 0.5):
    return ((1.0 - beta) * xc) + (beta * xh)

def pomak_prema_xl(X, l, sigma = 0.5):
    for i in range(len(X)):
        X[i] = sigma * (X[i] + X[l])
    return X



# - HOOKE - JEEVES -

# x0 - početna točka
# xb - bazna točka
# xp - početna točka pretraživanja
# xn - točka dobivena pretraživanjem

def hooke_jeeves(x0, func, dx = 0.5, epsilon = 1e-6):

    xp = x0.copy()
    xb = x0.copy()
#     dx = Matrica(x0.br_red, x0.br_stup, np.full((x0.br_red, x0.br_stup), 0.5), True)
                 
    while True:
        xn = istrazi(func, xp, dx)
        
#         #ispis rjesenja svakog pojednog koraka
#         print("xb: ", xb.elementi[0][0], "xp: ", xp.elementi[0][0], "xn: ", xn.elementi[0][0], "dx: ", dx)
#         print("fxb: ", f_xb, "fxp: ", f_xp, "fxn: ", f_xn)
#         print("----------------------------------")
        
        f_xn = func.vrijednost(xn)
        f_xb = func.vrijednost(xb)
        f_xp = func.vrijednost(xp)
        
#         print(xb[0], '-' ,xp[0], '-',xn[0], "dx: ", dx)
#         print("{:4.2f} - {:4.2f} - {:4.2f}".format(f_xb, f_xp, f_xn))
#         print("----------------------------------")
        
        if f_xn < f_xb:
            xp = 2*xn - xb
            xb = xn.copy()
        else:
            dx = dx / 2
            xp = xb.copy()
        
        #uvjet zaustavljanja
        if dx < epsilon:
            return xb
        
        
def istrazi(func, xp, dx):
    x = xp.copy()
    for i in range(x.br_stup):
        p = func.vrijednost(x)
#         print("hj: ", x)
        x[0][i] += dx
#         print("hj: ", x, " += dx")
        n = func.vrijednost(x)
        if n > p:
            x[0][i] -= 2 * dx
#             print("hj: ", x, " -= 2 * dx")
            n = func.vrijednost(x)
            if n > p:
                x[0][i] += dx
#                 print("hj: ", x, " += dx")
#     print("hj: new_xn: ", x)
    return x



# - FUNKCIJA -

class Funkcija:
    
    def __init__(self, f):
        self.br_poziva = 0
        self.f = f
        self.vrijednosti = {}
    
    def reset(self):
        self.vrijednosti = {}
        sel1.br_poziva = 0
    
    def vrijednost(self, x):
        try:
#             print("self.vrijednosti: ", self.vrijednosti)
#             print("str x:", str(x))
            return self.vrijednosti[str(x)]
        except: 
            nova_vrijednost = self.f(x)
            self.br_poziva += 1
            self.vrijednosti[str(x)] = nova_vrijednost
            return nova_vrijednost

# funkcije 

def f1(vektor):
    x = vektor.elementi.flatten()
    # print("unutar fje -> x1, x2: ", x[0], x[1])
    result = 100 * (x[1] - (x[0])**2)**2 + (1 - x[0])**2
    return result

def f2(vektor):
    x = vektor.elementi.flatten()
    return (x[0] - 4)**2 + 4 * (x[1] - 2)**2

def f3(vektor):
    x = vektor.elementi.flatten()
    return (x[0] - 2)**2 + (x[1] + 3)**2

def f4(vektor):
    x = vektor.elementi.flatten()
    return (x[0] - 3)**2 + (x[1])**2


# gradijenti funkcija

def gf1(vektor):
    x = vektor.elementi.flatten()
    d1 = 400 * x[0] * (x[0]**2 - x[1]) + 2 * (x[0] - 1)
    d2 = 200 * (x[1] - x[0]**2)
    return [d1, d2]

def gf2(vektor):
    x = vektor.elementi.flatten()
    d1 = 2 * (x[0] - 4)
    d2 = 8 * (x[1] - 2)
    return [d1, d2]

def gf3(vektor):
    x = vektor.elementi.flatten()
    d1 = 2 * (x[0] - 2)
    d2 = 2 * (x[1] + 3)
    return [d1, d2]

def gf4(vektor):
    x = vektor.elementi.flatten()
    d1 = 2 * (x[0] - 3)
    d2 = 2 * x[1]
    return [d1, d2]


# Hessove matrice

def hf1(vektor):
    x = vektor.elementi.flatten()
    d11 = 400*(3*x[0]**2 - x[1]) + 2
    d12 = -400*(x[0])
    d21 = -400*(x[0])
    d22 = 200
    inverz = np.linalg.inv([[d11, d12],[d21,d22]])
    return inverz.tolist()

def hf2(vektor):
    x = vektor.elementi.flatten()
    d11 = 2
    d12 = 0
    d21 = 0
    d22 = 8
    inverz = np.linalg.inv([[d11, d12],[d21,d22]])
    return inverz.tolist()

# Ogranicenja

def o1(vektor):
    x = vektor.elementi.flatten()
    return x[1] - x[0]

def o2(vektor):
    x = vektor.elementi.flatten()
    return 2 - x[0]

def o31(vektor):
    x = vektor.elementi.flatten()
    return x[0] + 100

def o32(vektor):
    x = vektor.elementi.flatten()
    return x[1] + 100

def o41(vektor):
    x = vektor.elementi.flatten()
    return 100 - x[0]

def o42(vektor):
    x = vektor.elementi.flatten()
    return 100 - x[1]

def o5(vektor):
    x = vektor.elementi.flatten()
    return 3 - x[0] - x[1]

def o6(vektor):
    x = vektor.elementi.flatten()
    return 3 + 1.5 * x[0] - x[1]


# - FUNKCIJE 2. LABORATORIJSKE VJEZBE -

# def f3(vektor):
#     x = vektor.elementi.flatten()
#     result = 0
#     for i in range(len(x)):
#         result += (x[i] - (i + 1))**2
#     return result

# def f4(vektor):
#     x = vektor.elementi.flatten()
#     return abs((x[0] - x[1]) * (x[0] + x[1])) + math.sqrt(x[0]**2 + x[1]**2)

# def f6(v):
#     x = v.elementi.flatten()
#     suma = 0
#     try:
#         for i in range(len(x)):
#             suma += x[i]**2
#     except TypeError: 
#         suma = v**2
#     sinus = math.sin(math.sqrt(suma))
#     razlomak = ((sinus**2) - 0.5) / ((1 + 0.001 * suma)**2)
#     return 0.5 + razlomak