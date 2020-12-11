from numpy import *
from Interval import *
import math
from scipy import stats
import Util
from PreciseNonNegativeReal import *
import pdb

class ModelParams(object):
    state = {'DEL':0, 'DIPLOID':1, 'DUP':2}
    statestr = ['DEL', 'DIPLOID', 'DUP'] 

    def __init__(self, mode, params, intervals, het_nums, tagsnp):
        self.mode = mode
        self.het_nums = het_nums
        self.tagsnp = tagsnp
        self._p = float(params[0])
        self._e = int(params[1])
        self._d = int(params[2]) * Interval.KB

        self._q = 1 / float(self._e)
        self._lambda = 1 / float(self._d) 
        
        num_states = self.getNumHiddenStates()
        self._transitionMat = zeros((num_states, num_states), dtype=float)
        DEL = self.state['DEL']
        DIPLOID = self.state['DIPLOID']
        DUP = self.state['DUP']

        self._transitionMat[DEL, DEL] = 1 - self._q
        self._transitionMat[DEL, DIPLOID] = self._q
        self._transitionMat[DEL, DUP] = 0

        self._transitionMat[DIPLOID, DEL] = self._p
        self._transitionMat[DIPLOID, DIPLOID] = 1 - 2 * self._p
        self._transitionMat[DIPLOID, DUP] = self._p

        self._transitionMat[DUP, DEL] = 0 
        self._transitionMat[DUP, DIPLOID] = self._q 
        self._transitionMat[DUP, DUP] = 1 - self._q 
        
        self._mean = []
        self._sd = []

        self._mu = []
        self._fi = []

        for i in range(3, 9, 2) :
            m = float(params[i])
            sd = float(params[i+1])
            self._mean.append(m)
            self._sd.append(sd)

        # For Negative binomial distribution
        if self.mode == 'baseline' or mode == 'reference' or mode == 'single' or mode == 'single-sample':
            self._mu.append(params[-2]/float(3))
            self._mu.append(params[-2])
            self._mu.append(params[-2]*float(3))
            self._fi.append(params[-1]/float(3))
            self._fi.append(params[-1])
            self._fi.append(params[-1]/float(3))

        #calculate the distances between targets which will be used to
        #calculate transition probalities
        self.interTargetDistances = []
        for i in range(len(intervals) - 1):
            self.interTargetDistances.append(intervals[i].distance(intervals[i+1]))

        return None    

    def transitionProb(self, i, j, t1):
        distInBases = self.interTargetDistances[t1] 

        f = 0
        if distInBases < float("inf") :
            f = math.exp(-(distInBases * self._lambda))

        return f * self._transitionMat[i, j] + (1-f) * self.getDiploidTransProb()[j]

    #return the transition matrix at t1->t1+1
    def transitionMatirx(self, t1) :
        tm = Util.getMatrix(self.getNumHiddenStates(), self.getNumHiddenStates())

        for i in range(self.getNumHiddenStates()) :
            for j in range(self.getNumHiddenStates()):
                tm[i][j] = PreciseNonNegativeReal(self.transitionProb(i, j, t1))

        return tm
    
    #value: read depth
    def getEmissProbs(self, mode, value) :
        if isinstance(value, str):
            value = float(value)
        emissProbs = []
        for i in range(self.getNumHiddenStates()) :
            if mode == 'svd' or mode == 'SVD' or mode == 'pooled' or mode == 'pooled-sample':
                emissProbs.append(PreciseNonNegativeReal(stats.norm.pdf(value, self._mean[i], self._sd[i])))
            # emissProbs.append(PreciseNonNegativeReal(stats.norm.pdf(value, self._mu[i], self._std[i])))
            elif mode == 'baseline' or mode == 'reference' or mode == 'single' or mode == 'single-sample':
                emissProbs.append(PreciseNonNegativeReal(self.calculate_NB_Prob(value,self._mu[i],self._fi[i])))
        return emissProbs

    def getEmissProbsWithSNP(self, mode, value, t) :
        if isinstance(value, str):
            value = float(value)
        emissProbs = []

        for i in range(self.getNumHiddenStates()) :
            if i == 0:
                if mode == 'svd' or mode == 'SVD' or mode == 'pooled' or mode == 'pooled-sample':
                    emissProbs.append(PreciseNonNegativeReal(self.getCH(0, t) * self.getCT(0, t) * stats.norm.pdf(value, self._mean[i], self._sd[i])))
                elif mode == 'baseline' or mode == 'reference' or mode == 'single' or mode == 'single-sample':
                    emissProbs.append(PreciseNonNegativeReal(self.getCH(0, t) * self.getCT(0, t) * self.calculate_NB_Prob(value,self._mu[i],self._fi[i])))
            else:
                if mode == 'svd' or mode == 'SVD' or mode == 'pooled' or mode == 'pooled-sample':
                    emissProbs.append(PreciseNonNegativeReal(self.getCH(i, t) * self.getCT(i, t) * stats.norm.pdf(value, self._mean[i], self._sd[i])))
                elif mode == 'baseline' or mode == 'reference' or mode == 'single' or mode == 'single-sample':
                    emissProbs.append(PreciseNonNegativeReal(self.getCH(i, t) * self.getCT(i, t) * self.calculate_NB_Prob(value,self._mu[i],self._fi[i])))

        return emissProbs

    # use these probabilities to penalize the heterozgosity at the regions of deletions
    # Jixuan's original
    # def getCH(self, i, t):
    #     if i ==0 or i == 1:
    #         return 10 ** (i - 3 - self.het_nums[t])
    #     elif i >= 2:
    #         return (1 - self.getCH(0, t) - self.getCH(1, t)) / 10

    # microtan modified
    def getCH(self, i, t):
        if self.het_nums[t] > 0:
            if i == 0:
                # return 10 ** (-30*self.het_nums[t])
                return 10 ** (-3 - self.het_nums[t])
            elif i >= 1:
                # return (1 - self.getCH(0, t)) / 10
                return 1
        else:
            return 1

    def getCT(self, i, t):
        if self.tagsnp[t] > 0.8:
            if i == 0:
                return 10 ** (3 + self.tagsnp[t])
            elif i >= 1:
                return 1
        else:
            return 1
    
    def getInitProbs(self):
        v = self._transitionMat[self.state['DIPLOID']]
        return v
        norm = linalg.norm(v)
        if norm==0: 
            return v
        return v/norm

    def getDiploidTransProb(self):
        return self._transitionMat[self.state['DIPLOID'], :]

    def calculate_NB_Prob(self,x,mu,fi):
        prob_nb = (math.e)**(math.lgamma(x+1/fi)-math.lgamma(x+1)-math.lgamma(1/fi)+x*math.log(mu/(mu+1/fi))-1/fi*math.log(1+mu*fi))
        return prob_nb

    @staticmethod
    def getNumHiddenStates():
        return len(ModelParams.state)

