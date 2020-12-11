from numpy import *
from ModelParams import *
import Util
from PreciseNonNegativeReal import *
import copy
import time

class Model (object) :
    def __init__(self, _mode, _model_params, _observations):
        self._mode = _mode
        self._model_params = _model_params
        self._observations = _observations
    
    # IMPORTANT: the elements of alpha, beta, gamma are all
    # PreciseNonNegativeReal type
    def forwardBackward_Viterbi(self, mode, if_snp):
        s = ModelParams.getNumHiddenStates()
        o = len(self._observations)
        alpha = Util.getMatrix(o, s) 
        beta = Util.getMatrix(o, s)
        gamma = Util.getMatrix(o, s)
        
        # compute forward probabilities
        for t in range(o):
             #print 'Calculating forward message ' + str(t),
            self.calcForwardMessage(mode, t, alpha, if_snp)
            #print alpha[t]
        
        #compute backward probabilities
        for t in range(o)[::-1] :
            #print 'Calculating backward message ' + str(t),
            self.calcBackwardMessage(mode, t, beta, if_snp)
            #print beta[t]
        
        for x in range(o) :
            if type(self._observations[x]) in (int,float) :
                for y in range(s):
                    gamma[x][y] = PreciseNonNegativeReal(alpha[x][y] * beta[x][y])
                    # print x,self._observations[x],alpha[x][y],beta[x][y],gamma[x][y]
                    # time.sleep(2)
            else:
                print x,self._observations[x]


        _vpath = []
        for t in range(o) :
            _vpath.append(Util.maxIndex(gamma[t]))

        vpathlist = []
        for v in _vpath :
                vpathlist.append(ModelParams.statestr[v])

        return vpathlist

    def calcForwardMessage(self, mode, t, alpha, if_snp):
        forwardMess = alpha[t]
        s = ModelParams.getNumHiddenStates()

        if t == 0:
            for i in range(s):
                forwardMess[i] = PreciseNonNegativeReal(self._model_params.getInitProbs()[i])
        else:
            incomingForwardMess = alpha[t-1]
            transMat = self._model_params.transitionMatirx(t-1)

            # convert 'incomingForwardMess' to matirx
            # convert 'forwardMess' to list
            forwardMess = Util.mulMatrix([incomingForwardMess], transMat)[0]
        
        self.multiplyMessageTimesEmissProbs(mode, forwardMess, t, if_snp)
        alpha[t] = forwardMess

    def calcBackwardMessage(self, mode, t, beta, if_snp) :
        backwardMess = beta[t]
        o = len(self._observations)
        s = ModelParams.getNumHiddenStates()

        if t == o - 1 :
            backwardMess = Util.getMatrix(1, s, 1)[0]
        else :
            incomingBackwardMess = copy.deepcopy(beta[t+1])
            self.multiplyMessageTimesEmissProbs(mode, incomingBackwardMess, t+1, if_snp)
            transMatTranspose = Util.transposeMatirx(self._model_params.transitionMatirx(t))

            # convert 'incomingBackwardMess' to matirx
            # convert 'backwardMess' to list
            backwardMess = Util.mulMatrix([incomingBackwardMess], transMatTranspose)[0]

        beta[t] = backwardMess
            
    def multiplyMessageTimesEmissProbs(self, mode, message, t, if_snp):
        emissProbs = self._model_params.getEmissProbs(mode, self._observations[t])
        if if_snp:
            try:
                emissProbs = self._model_params.getEmissProbsWithSNP(mode, self._observations[t], t)
            except:
                #Renjie modified at Jun 17th, 2016
                pass
                # print 'NO SNP information!, then keep emission probability same as No SNP mode.'
            
        for i in range(ModelParams.getNumHiddenStates()) :
            message[i] *= emissProbs[i]
