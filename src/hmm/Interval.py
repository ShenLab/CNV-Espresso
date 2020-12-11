class Interval(object):
    KB = 1000
    def __init__(self, _chr, _bp1, _bp2):
        self._chr = _chr
        self._bp1 = int(_bp1)
        self._bp2 = int(_bp2)


    def distance(self, another):
        if self._chr != another._chr:
            return float("inf")
        if self._bp2 < another._bp1:
            return another._bp1 - self._bp2
        elif another._bp2 < self._bp1:
            return self._bp1 - another._bp2

        return 0

    def getInfo(self):
        return self._chr+':'+str(self._bp1)+'-'+str(self._bp2)
