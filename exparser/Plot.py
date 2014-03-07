#-*- coding:utf-8 -*-

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

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from exparser.TangoPalette import *

plotFolder = 'plot'
if '--clear-plot' in sys.argv and os.path.exists(cacheFolder):
	print 'Removing plot folder (%s)' % plotFolder
	shutil.rmtree(plotFolder)
plt.rc('font', family='Arial', size=10)

# Some pre-defined sizes
xs = 4, 4
s = 6, 6
r = 8, 8
w = 12, 8
h = 8, 12
l = 12, 12
xl = 16, 16

def new(size=r):

	"""
	Creates a new figure.

	Keyword arguments:
	size	--	The figure size. (default=r)

	Returns:
	A matplotlib figure.
	"""

	return plt.figure(figsize=size)

def regress(x, y, annotate=True, symbol='.', linestyle='--', symbolColor= \
	blue[1], lineColor=blue[1]):

	"""
	Creates a regression plot.

	Arguments:
	x				--	An array or list for the X data.
	y			--	An array or list for the Y data.

	Keyword arguments:
	annotate	--	Indicates whether the correlation and p-value should be
					marked in the plot. (default=True)
	"""

	from scipy.stats import linregress
	s, i, r, p, se = linregress(x, y)
	plt.plot(x, y, symbol, color=symbolColor)
	xData = np.array([min(x), max(x)])
	yData = i + s*xData
	plt.plot(xData, yData, linestyle, color=lineColor)
	if annotate:
		plt.text(0.05, 0.95, 'r = %.3f, p = %.3f' % (r, p), ha='left', \
			va='top', transform=plt.gca().transAxes)

def save(name, show=False):

	"""
	Saves the current figure to the correct folder, depending on the active
	experiment.

	Arguments:
	name	--	The name for the figure.

	Keyword arguments:
	show	--	Indicates whether the figure should be shown as well.
				(default=False)
	"""

	try:
		os.makedirs(os.path.join(plotFolder, 'svg'))
	except:
		pass
	try:
		os.makedirs(os.path.join(plotFolder, 'png'))
	except:
		pass
	pathSvg = os.path.join(plotFolder, 'svg', '%s.svg' % name)
	pathPng = os.path.join(plotFolder, 'png', '%s.png' % name)
	plt.savefig(pathSvg)
	plt.savefig(pathPng)
	if show:
		plt.show()
	else:
		plt.clf()
