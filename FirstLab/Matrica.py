import numpy as np

class Matrica:
    
    def __init__ (self, brRed, brStup, identity = False, elementi = None):
        self.brRed = brRed
        self.brStup = brStup
        self.epsilon = 1e-7
        
        if elementi is None :
            self.elementi = np.zeros((brRed, brStup))
        else:
            self.elementi = np.array(elementi)
        
        if identity == True:
            for i in range(brRed):
                self.elementi[i][i] = 1.0
            
    def copy(self):
        return Matrica(self.brRed, self.brStup, self.elementi.deepCopy())
    
    # Overloading functions
    
    def __str__(self):
        return self.elementi
    
    def __eq__(self, other):
        if self.brRed != other.brRed or self.brStup != other.brStup: return False
        for i in range(self.brRed):
            for j in range(self.brStup):
                if abs(self.elementi[i][j] - other.elementi[i][j]) > self.epsilon: return False
        return True
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __str__(self):
        return f'{self.elementi}'
    
    def __add__(self, other):
        result = Matrica(self.brRed, self.brStup, False)
        for i in range(self.brRed):
            z = list(zip(self.elementi[i], other.elementi[i]))
            m = list(map(sum, z))
            result.elementi[i] = m
#         matrica = Matrica(self.brRed, self.brStup, False, result)
        return result
    
    def __sub__(self, other):
        result = Matrica(self.brRed, self.brStup, False)
        for i in range(self.brRed):
            z = list(zip(self.elementi[i], other.elementi[i]))
            m = [pair[0] - pair[1] for pair in z]
            result.elementi[i] = m
        return result
    
    def __mul__(self, other):
        t = type(other)
        if (t == int or t == float):
            result = Matrica(self.brRed, self.brStup, False)
            for i in range(self.brRed):
                m = [x * other for x in self.elementi[i]]
                result.elementi[i] = m
            return result
        else:
            result = Matrica(self.brRed, other.brStup, False)
            for i in range(self.brRed):
                resultRow = []
                for j in range(other.brStup):
                    stupac = other.dohvatiStupac(j)
                    z = list(zip(self.elementi[i], stupac))
                    m = [pair[0]*pair[1] for pair in z]
                    resultRow.append(sum(m))
                result.elementi[i] = resultRow
        return result
    
    def __iadd__(self, other):
        return self + other
        
    def __isub__(self, other):
        return self - other
    
    def transpose(self):
        result = Matrica(self.brRed, self.brStup, False)
        for i in range(self.brRed):
            for j in range(self.brStup):
                result.elementi[i][j] = self.elementi[j][i]
        return result
    
    def dohvatiStupac(self, i):
        return [self.elementi[j][i] for j in range(self.brRed)]
    
    
    ## SUPSTITUCIJA
    
    def supst_unaprijed(self, b):
        y = b.elementi.copy()
        for i in range(self.brRed-1):
            for j in range(i+1, self.brRed):
                y[j] -= self.elementi[j][i] * y[i]
        return Matrica(1, b.brStup, False, y)
    
    def supst_unazad(self, y):
        x = y.elementi.copy()
        for i in reversed(range(self.brRed)):
            x[i] = x[i] / self.elementi[i][i]
            if(i <= 0): break
            for j in range(i):
                x[j] -= self.elementi[j][i] * x[i]
        return x
        