import scipy.special as special
import scipy.optimize as optimize
import numpy as np
import mpmath
import pdb

class ParameterEstimation(object):
	def __init__(self, p = 0.1, r = 10):
		nbin_mpmath = lambda k, p, r: mpmath.gamma(k + r)/(mpmath.gamma(k+1)*mpmath.gamma(r))*np.power(1-p, r)*np.power(p, k)
		self.nbin = np.frompyfunc(nbin_mpmath, 3, 1)
		self.p = p
		self.r = r

	def mleFun(self, par, data, sm):
		'''
		Objective function for MLE estimate according to
		https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation

		Keywords:
		data -- the points to be fit
		sm -- \sum data / len(data)
		'''
		p = par[0]
		r = par[1]
		n = len(data)
		f0 = sm/(r+sm)-p
		f1 = np.sum(special.psi(data+r)) - n*special.psi(r) + n*np.log(r/(r+sm))
		return np.array([f0, f1])

	def fit(self, data, remove_lower, remove_upper):
		p = None
		r = None
		data = [float(i) for i in data]
		data.sort()
		data_len = len(data)
		data_selected = data[int(data_len*remove_lower):int(data_len*remove_upper)]

		if p is None or r is None:
			av = np.average(data_selected)
			va = np.var(data_selected)
			r = (av*av)/(va-av)
			p = (va-av)/(va)
		sm = np.sum(data_selected)/len(data_selected)
		x = optimize.fsolve(self.mleFun, np.array([p, r]), args=(data_selected, sm))
		self.p = x[0]
		self.r = x[1]
		mu = self.p*self.r/(1-self.p)
		fi = 1/self.r
		return np.array([mu,fi])

	def gaussian_fit(self,data,remove_lower,remove_upper):
		data = [float(i) for i in data]
		data.sort()
		data_len = len(data)
		data_selected = data[int(data_len*remove_lower):int(data_len*remove_upper)]
		mu = np.average(data_selected)
		sd = np.std(data_selected)
		return np.array([mu, sd])

	def pdf(self, k):
		return self.nbin(k, self.p, self.r).astype('float64')
