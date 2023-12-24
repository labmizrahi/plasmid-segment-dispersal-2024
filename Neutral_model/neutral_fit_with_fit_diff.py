#!/usr/bin/env python

'''
@author: jonathanfriedman
@date: 2017-07-01

Module for estimating the average coverage per position for plasmids, based on the coverage distribution along the plasmid.

Dependencies: numpy, scipy, pandas
The easiest way to istall these are to use a python distribution such as Anaconda.
'''

import numpy as np
from numpy import (unravel_index, argmax, ones, corrcoef, cov, r_, 
                   diag, sqrt, where, nan, log)
import pandas as pd
from scipy import optimize as spo
from scipy import stats as stats

def objective(M, *args):
	'''
	Objective function for fitting M.
	Returns the sum of squared errors between predicted and observed
	detection frequencies of all OTUs.
	'''
	p, det_lim, F_obs  = args
	F_expt, F_var = expected_detection_freq(M, p, det_lim)
	error = np.sum((F_obs-F_expt)**2)
	return error


def R2(obs,fit):
	'''
	Calculate the R^2 (coefficient of determination) of a fit.
	Define: 
	    SS_fit = sum(obs-fit)^2
	    SS_tot = sum(obs-<obs>)^2
	    R2 = 1 - SS_fit/SS_tot
	See: 
	    http://en.wikipedia.org/wiki/Coefficient_of_determination        
	'''
	SS_fit  = np.sum((obs-fit)**2)
	obs_avg = np.mean(obs)
	SS_tot  = np.sum((obs-obs_avg)**2)
	R2 = 1 - SS_fit/SS_tot
	return R2
    
    

def expected_detection_freq(M, p, det_lim):
	'''
	The expected detection frequency of each OTU, given the parameter M = mN_T.
	This function assumes neutrality, larger population size (continuous limit), and unbiased sampling from
	the local community.
	The parameter p (fraction in the meta community) is estimated from the counts.

	Beta model:
	    The distribution of fractions of each OTU is a binomial with parameters a_i, b_i as given above.
	    There's a detection threshold in terms of a minimal fraction. 
	    The expected detection frequency of each OTU is the probability of it's fraction exceeding the detection threshold.
	    This is given by
	        <Detection Freq> = 1 - cdf_Beta(th | a_i,b_i)
	'''
	otus    = p.index
	F_expt = pd.Series(np.zeros(len(otus)), otus)
	F_var  = pd.Series(np.zeros(len(otus)), otus)
	a = p * M
	b = (1-p)/p * a
	F_expt = pd.Series(stats.beta(a,b).sf(det_lim), otus) 
	return (F_expt, F_var)


class DetFreq(object):
	'''
	Class for fitting different community assembly models based on the detection frequency of all OTUs.
	'''
	def __init__(self,p, F_obs, det_lim, **kwargs):
		'''
		Constructor
		'''
		self.p       = p         #estimated proportion in metacommunity
		self.F_obs   = F_obs     #observed detection frequency
		self.det_lim = det_lim
		self.M       = None
		self.error   = None
		self.F_expt  = None
		self.R2      = None

	    
	def set_Fexpt(self, model, Fexpt=None):
		if Fexpt is None:
			M = self.M
			p = self.p
			det_lim = self.det_lim
			Fexpt = expected_detection_freq(M, p, det_lim)
		else:
			Fexpt = Fexpt 
		self.Fexpt[model] = Fexpt
	              
	def set_R2(self, model):
		self.R2 = R2(self.F_obs, self.F_expt)

	def estimate_M(self,  **kwargs):
		'''
		Main function used to find best fitting values of M (=mN_T).
		'''

		import scipy.optimize as spo
		#set up optimization problem
		tol = kwargs.get('tol', 1e-5) 
		lb  = kwargs.get('lb', tol)
		M0  = kwargs.get('M0', 1e1)
		iprint = kwargs.get('iprint', 0)
		#set up problem
		args = (self.p, self.det_lim, self.F_obs)
		sol = spo.minimize_scalar(objective, bounds=(lb, 1e5), method='bounded', args=args)    
		M_opt = sol.x
		error = sol.fun

		F_expt,F_var = expected_detection_freq(M_opt, self.p, self.det_lim)

		self.F_expt = F_expt
		self.M = M_opt
		self.error = error
		self.R2 = R2(self.F_obs, self.F_expt)

	def plot_fit(self, **kwargs):
		import pylab
		## parse input args
		kwargs_obs  = kwargs.get('kwargs_obs',{})
		kwargs_obs.setdefault('marker','s')
		kwargs_obs.setdefault('linestyle','')
		kwargs_obs.setdefault('label','Observed')
		kwargs_expt = kwargs.get('kwargs_expt',{})
		kwargs_expt.setdefault('linewidth',2)
		xlog = kwargs.get('xlog',True)
		ylog = kwargs.get('ylog',False)
		label_fs = kwargs.get('label_fs',16)
		legend     = kwargs.get('legend',True)
		legend_fmt = kwargs.get('legend_fmt','%.2f')
		legend_loc = kwargs.get('legend_loc','best')
		show_R2 = kwargs.get('show_R2',True)
		R2_fmt  = kwargs.get('R2_fmt','%.2f')
		new_fig = kwargs.get('new_fig',True)
		file = kwargs.get('file',None)
		show = kwargs.get('show',True)

		## plot observed freqs
		p = self.p.sort_values()
		x = p.values
		if xlog: x = np.log10(x)
		obs  = self.F_obs.reindex_like(p)
		y = obs.values
		if ylog: y = np.log10(y)
		if new_fig: pylab.figure()
		pylab.plot(x, y, **kwargs_obs)

		## plot expected freqs 
		expt = self.F_expt.reindex_like(p)
		y = expt.values
		R2 = self.R2
		if ylog: y = np.log10(y)
		label = ('Neutral, mN=' + legend_fmt) %(self.M)                
		if show_R2: label += ('\n $R^2$ = ' + R2_fmt) %R2
		pylab.plot(x, y, label=label, **kwargs_expt)
		if xlog: xlabel = 'log10(Mean relative abundance)'
		else:    xlabel = 'Mean relative abundance'
		if ylog: ylabel = 'log10(Detection Frequency)'
		else:    ylabel = 'Detection Frequency'
		pylab.xlabel(xlabel, fontsize=label_fs)
		pylab.ylabel(ylabel, fontsize=label_fs)
		if legend: pylab.legend(loc='lower right', fontsize=label_fs)
		if file is not None: pylab.savefig(file)
		if show: pylab.show()



def data_parser(x):
	x = str(x)
	x = x.replace(',', '.')
	return float(x)


def read_data(file_name, presence_th=1e-6):
	'''
	Reads the plasmid coverage profile.

	Inputs
	------
	file_name: str
		The name of the file containing the abundances of plasmid (rows) across samples (columns).

	presence_th: float
		Abundance threshold for saying that a plasmid is present in a given sample. 
		This is imposed since the prcedure for finding plasmid abundances returns small positive values even when no reads are found for a given plasmid.

	Returns
	-------
	fracs: DataFrame
		A frame countaing the relative abundances.
		Rows and columns represent plasmids and samples, correspondignly. 
	'''
	data_raw = pd.read_csv(file_name, header=0, sep='\t', index_col=0).squeeze().applymap(data_parser).dropna()
	temp = data_raw.applymap(lambda x:x if x>presence_th else 0)
	data_fitered = temp.loc[temp.sum(axis=1)>0]
	detection = data_fitered.apply(lambda x: (x>presence_th).sum(), axis=1) # number of samples easch plasmid is present at

	fracs = data_fitered.apply(lambda x:x/x.sum()) #normalize to make sure working with fractions
	return fracs
	
	
def fit_model(fracs, det_lim=None, plot=False, plot_kwargs={}):
	'''	
	Computes the avergae coverage per positions for plasmid in multiple samples..

	Inputs
	------
	fracs: DataFrame
		A frame countaing the relative abundances.
		Rows and columns represent plasmids and samples, correspondignly. 
	det_lim: float (default=None) 
		Detection limit of plasmids. If None given, infer from observed plasmid abundances.
	plot: bool (default=Fasle) 
		Flag for making plot of mean fraction vs detection frequency of observed data and fit line.
	plot_kwargs: dict (default={})
		Additional arguments to pass to plotting function.

	Returns
	-------
	M: float
		Fitted migration parameter value.
	R2: float
		R-squared value of fit.  
	'''
	n = float(fracs.shape[1]) # number of samples
	p = fracs.mean(axis=1) #mean fraction
	detection = fracs.apply(lambda x: (x>0).sum(), axis=1) # number of samples easch plasmid is present at
	if det_lim is None:
		det_lim = fracs[fracs>0].min().min()/10
	fit = DetFreq(p, detection/n, det_lim)
	fit.estimate_M()
	R2_obs = fit.R2
	if plot:
		fit.plot_fit(**plot_kwargs)

	return fit.M, fit.R2, fit


def write_fit_file(ofile, M, R2):
	o = pd.Series([M, R2], index=['M', 'R2'])
	o.to_csv(ofile)

if __name__ == '__main__':
	import argparse
	kwargs = {}
	description  = """Compute the average coverage per position of a given plasmid across multiple samples. 
					The average coverage is computed by finding a values that minimizes the log-deviation of the coverage of all positions from this average.
					Using the log of the deviation makes this average less sensitive to outliers (region with unusually high/low coverage). 
					This is an intermediate between the regular mean, which minimizes the squared deviation, and the median, which minimizes the deviaiton sign differences (deviation^0).
					Mathematically, the average coverage is given by the solution to argmin_c \sum_i log(1+abs(c-c_i)).



					"""
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("ifile",
	                  help="File containnig the abundances of plasmid (rows) across samples (columns).")
	parser.add_argument("-p", "--plot", action="store_true", default=False,
	                  help="Flag for making plot of mean fraction vs detection frequency of observed data and fit line (default=False).")
	parser.add_argument('-o', "--ofile", 
						help="File to which the output would be written.\nIf not provided, the input file name will be used with a '.fit.csv' suffix.")
	parser.add_argument('-d', "--dfile", 
						help="File to which the difference between the fit and observed detection frequencies would be written.\nIf not provided, the input file name will be used with a '.fit.diff.xlsx' suffix.")
	args = parser.parse_args()
	fracs = read_data(args.ifile)
	M, R2, fit = fit_model(fracs, plot=args.plot)
	print(M, R2)

	if args.ofile is None:
		ofile = args.ifile + '.fit.csv'
	else:
		ofile = args.ofile
	write_fit_file(ofile, M, R2)

	## 
	if args.dfile is None:
#line 297 wasn't working because a module was missing so I changed it and the file extension as well
#	dfile = args.ifile + '.fit.diff.xlsx'
		dfile = args.ifile + '.fit.diff.csv'
	else:
		dfile = args.dfile
	x = fit.p.sort_values()
	y = fit.F_obs.reindex_like(x)
	y_fit = fit.F_expt.reindex_like(x)
	d = (y-y_fit)
#this next line wasn't working because a module was missing so I changed it
#	d.to_excel(dfile)
	d.to_csv(dfile)
     

