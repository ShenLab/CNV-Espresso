from PreciseNonNegativeReal import *

def getMatrix(x, y, default=PreciseNonNegativeReal(0)):
    m = []
    for i in range(0, x):
        tmp = []
        for j in range(0, y):
            tmp.append(PreciseNonNegativeReal(default))
        m.append(tmp)

    return m

def mulMatrix(m1, m2):
    m1_x = len(m1)
    m1_y = len(m1[0])

    m2_x = len(m2)
    m2_y = len(m2[0])

    if m1_y != m2_x:
        raise Exception('Matrix multiplication error')

    result = getMatrix(m1_x, m2_y)

    for i in range(m1_x):
        for j in range(m2_y):
            tmp = PreciseNonNegativeReal(0)
            for a in range(m1_y):
                tmp += m1[i][a] * m2[a][j]
            result[i][j] = PreciseNonNegativeReal(tmp)

    return result

def transposeMatirx(m):
    x = len(m)
    y = len(m[0])
    m_t = getMatrix(y, x)

    for i in range(y):
        for j in range(x):
            m_t[i][j] = m[j][i]

    return m_t

def maxIndex(l):
    index = 0

    for i in range(len(l)):
        if l[i] > l[index]:
            index = i

    return index
