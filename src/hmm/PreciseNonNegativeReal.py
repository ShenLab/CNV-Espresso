from math import *
from numpy import *

class PreciseNonNegativeReal:
    REAL_ZERO = 0.0
    REAL_ONE = 1.0
    REAL_TEN = 10.0
    REAL_EPSILON = 1e-6
    REAL_INFINITY = float('inf')

    def __init__(self, d, isLog10 = False):
        if isinstance(d, float) or isinstance(d, int):
            if isLog10:
                self.log10value = d
            else:
                if d < 0:
                    raise Exception("the argument muse be non-negative")
                elif d == 0:
                    self.log10value = float('-inf')
                else:
                    self.log10value = log10(d)
        else:
            self.log10value = d.log10value

    def getValue():
        return pow(10, self.log10value)

    def getLog10Value():
        return self.log10value

    # = couldn't be overloaded

    # +
    def __add__(self, other):
        tmp = PreciseNonNegativeReal(self)
        tmp += other
        return tmp

    # -
    def __sub__(self, other):
        tmp = PreciseNonNegativeReal(self)
        tmp -= other
        return tmp

    # *
    def __mul__(self, other):
        tmp = PreciseNonNegativeReal(self)
        tmp *= other
        return tmp
    
    # /
    def __div__(self, other):
        tmp = PreciseNonNegativeReal(self)
        tmp /= other
        return tmp

    # ==
    def __eq__(self, other):
        return self.compareTo(other) == 0

    # > 
    def __get__(self, other):
        return self.compareTo(other) > 0

    # <
    def __lt__(self, other):
        return self.compareTo(other) < 0

    # >=
    def __ge__(self, other):
        return self.compareTo(other) >= 0
        
    # <=
    def __le__(self, other):
        return self.compareTo(other) <= 0
        
    # +=
    def __iadd__(self, other):
        self.log10value = self.addInLogSpace(self.log10value, other.log10value)
        return self

    # -=
    def __isub__(self, other):
        self.log10value = self.absSubLog(self.log10value, other.log10value)
        return self

    # *=
    def __imul__(self, other):
        self.log10value += other.log10value
        return self

    # /=
    def __idiv__(self, other):
        self.log10value -= other.log10value
        return self

    # print
    def __repr__(self):
        return '10^('+str(self.log10value)+')'

    # if x = log(a), y = log(b), return log(a+b)
    def addInLogSpace(self, x, y):
        if x == self.REAL_INFINITY or y == self.REAL_INFINITY:
            return REAL_EPSILON

        if x == -self.REAL_INFINITY:
            return y
        if y == -self.REAL_INFINITY:
            return x

        maxVal = 0
        negDiff = 0
        if x > y:
            maxVal = x
            negDiff =  y - x
        else:
            maxVal = y
            negDiff = x - y

        #x + log(1+e^(y-x)) = log(a) + log(1+e^(log(b)-log(a))) = log(a) + log(1+b/a) = log(a+b)
        return maxVal + log10(self.REAL_ONE + pow(self.REAL_TEN, negDiff))
    
    # if x = log(a), y = log(b), return log |a-b|
    def absSubLog(self, x, y):
        if x == -self.REAL_INFINITY and y == -self.REAL_INFINITY:
            return -self.REAL_INFINITY

        maxVal = 0
        negDiff = 0
        if x > y:
            maxVal = x
            negDiff =  y - x
        else:
            maxVal = y
            negDiff = x - y
        #x + log(1-e^(y-x)) = log(a) + log(1-e^(log(b)-log(a))) = log(a) + log(1-b/a) = a - b = |a-b|, since x >= y
        return maxVal + log10(self.REAL_ONE - pow(self.REAL_TEN, negDiff))
        
    def compareTo(self, other):
        return  self.compareReal(self.log10value, other.log10value)

    def compareReal(self, r1, r2):
        try:
            diff = r1 - r2
            if abs(diff) <= self.REAL_EPSILON:
                return 0 # r1 == r2

            if diff > self.REAL_ZERO:
                return 1
            else:
                return -1
        except Exception, e:
            raise e
        
        
