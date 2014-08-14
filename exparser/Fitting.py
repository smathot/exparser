#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file is part of exparser.

exparser is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

exparser is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with exparser.  If not, see <http://www.gnu.org/licenses/>.
"""

from exparser.TangoPalette import *
import numpy as np
from matplotlib import pyplot as plt

def exgauss(x, m=0., s=1., l=1.):

	"""
	desc:
		Fits an exponentially modified Gaussian distribution, suitable for
		modelling skewed RT distributions.

	arguments:
		x:
			desc:	The X data.
			type:	ndarray

	keywords:
		m:
			desc:	The mean of the Gaussian.
			type:	[int, float]
		s:
			desc:	The standard deviation of the Gaussian.
			type:	[int, float]
		l:
			desc:	The inverse of the exponential paramater.
			type:	[int, float]

	returns:
		desc:	The Y data.
		type:	ndarray
	"""

	from scipy.special import erfc
	m = float(m)
	s = float(s)
	l = float(l)
	return l/2*np.e**( (l/2) * (2*m+l*s**2-2*x) ) * \
		erfc( (m+l*s**2-x) / (np.sqrt(2)*s) )

def sigmoid(x, x0=0, k=1):

	"""
	A logistic sigmoid function, useful for fitting cumulative distributions.

	Arguments:
	x		--	The X data.

	Keyword arguments:
	x0		--	The .5 point. (default=0)
	k		--	The steepness. (default=1)

	Returns:
	A cumulative distribition from 0 to 1
	"""

	y = 1 / (1 + np.exp(-k*(x-x0)))
	return y

def HL1993(t, n=10.1, tMax=930, pMin=0, pMax=1):

	"""
	An exponential pupil model for a transient pupillary response, as
	described by Hoeks and Levelt (1993). Default parameters have been estimated
	empirically by Hoeks and Levelt (1993).

	Arguments:
	t		--	Time series.

	Keyword arguments:
	n			--	The number of layers (actually corresponds to n+1).
					(default=10.1)
	tMax		--	The time at which the response is maximal. (default=930.0)
	pMin		--	Indicates the minimum (initial) pupil size value.
					(default=0)
	pMax		--	Indicates the maximum (final) pupil size value.
					(default=1)

	Returns:
	An array with pupil size values.
	"""

	a = t**n * np.exp(-n*t/tMax)
	a /= a.max()
	a = a * (pMax-pMin) + pMin
	return a

def HL1993D(t, n=10.1, tMax=930, pMin=1., vMin=0, vMax=1):

	"""
	Applying the HL algorithm to the velocity profile of the signal.

	Arguments:
	t		--	Time series.

	Keyword arguments:
	n			--	The number of layers (actually corresponds to n+1).
					(default=10.1)
	tMax		--	The time at which the response is maximal. (default=930.0)
	pMin		--	Indicates the minimum (initial) pupil size value.
					(default=0)
	pMax		--	Indicates the maximum (final) pupil size value.
					(default=1)

	Returns:
	An array with pupil size values.
	"""

	p = pMin
	l = []
	for v in HL1993(t, n=n, tMax=tMax, pMin=vMin, pMax=vMax):
		l.append(p)
		p -= v
	return np.array(l)

def PLRExp(t, n=10.1, t0=250, tMax=930, ps=1, pd=.6):

	"""
	An exponential PLR function.

	Arguments:
	t			--	Time series.

	Keyword arguments:
	n			--	The number of layers (actually corresponds to n+1).
					(default=10.1)
	t0			--	The latency of the response. (default=250)
	tMax		--	The time at which the response is maximal. (default=930.0)
	pMin		--	Indicates the minimum (initial) pupil size value.
					(default=0)
	pMax		--	Indicates the maximum (final) pupil size value.
					(default=1)

	Returns:
	An array with pupil size values.
	"""

	a = np.exp(-n*(t-t0)/tMax)
	a[np.where(t < t0)] = 1
	a -= 1
	a *= pd
	a += ps
	return a

def PLRWeibull(t, t0=150, k=2, L=200, ps=1., pd=.6):

	"""
	Keyword arguments:
	t0		--	The latency of the PLR.
	k		--	Shape parameter.
	L		--	Decay parameter.
	"""

	a = np.exp( -((t-t0)/L)**k * np.log(2) )
	a[np.where(t<t0)] = 1
	a -= 1
	a *= pd
	a += ps
	return a

def fit(x, y, func, p0=None, plot=False, color=blue[2]):

	"""
	Uses `scipy.optimize.curve_fit()` to fit a set of data using any function.
	Optionally, the data is plotted to the currently active axis.

	Arguments:
	x		--	The X data.
	y		--	The Y data.
	func	--	The function to fit.

	Keyword arguments:
	p0		--	Initial guess for the parameters for `func`. (default=None)
	plot	--	Indicates whether the fit should be plotted. (default=False)
	color	--	Indicates the line color, in case plot==True. (default=blue[2])

	Returns:
	An (residuals, params) tuple if the fit was successful or (np.nan, np.nan)
	if it was not.
	"""

	from scipy.optimize import curve_fit
	try:
		popt, pcov = curve_fit(func, x, y, p0=p0)
	except:
		print 'Failed to fit!'
		return None, None
	plt.plot(x, func(x, *popt), color=color)
	return np.sqrt(np.sum((y - func(x, *popt))**2)) / len(x), popt

if __name__ == '__main__':

	i = 1
	t = np.linspace(0, 2000)
	for func in HL1993, HL1993D, PLRExp, PLRWeibull:
		plt.subplot(6,1,i)
		plt.plot(t, func(t))
		i += 1
	plt.show()


