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
if '--clear-plot' in sys.argv and os.path.exists(plotFolder):
	print 'Removing plot folder (%s)' % plotFolder
	import shutil
	shutil.rmtree(plotFolder)
plt.rc('font', family='liberation sans', size=10)

# Some pre-defined sizes
xs = 4, 4
s = 6, 6
ws = 6, 3
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

	fig = plt.figure(figsize=size)
	plt.subplots_adjust(left=.15, right=.9, bottom=.15, top=.9, wspace=.3,
		hspace=.3)
	return fig

def regress(x, y, annotate=True, symbol='.', linestyle='--', symbolColor= \
	blue[1], lineColor=blue[1], label=None):

	"""
	Creates a regression plot.

	Arguments:
	x				--	An array or list for the X data.
	y			--	An array or list for the Y data.

	Keyword arguments:
	annotate	--	Indicates whether the correlation and p-value should be
					marked in the plot. (default=True)

	Returns:
	The regression parameters as a (slope, intercept, correlation, p-value,
	standard error) tuple
	"""

	from scipy.stats import linregress
	s, i, r, p, se = linregress(x, y)
	plt.plot(x, y, symbol, color=symbolColor)
	xData = np.array([min(x), max(x)])
	yData = i + s*xData
	plt.plot(xData, yData, linestyle, color=lineColor, label=label)
	if annotate:
		plt.text(0.05, 0.95, 'r = %.3f, p = %.3f' % (r, p), ha='left', \
			va='top', transform=plt.gca().transAxes)
	return s, i, r, p, se

def save(name, folder=None, show=False, dpi=200):

	"""
	Saves the current figure to the correct folder, depending on the active
	experiment.

	Arguments:
	name	--	The name for the figure.

	Keyword arguments:
	folder	--	A name for a subfolder to save the plot or None to save directly
				in the plotFolder. (default=None)
	show	--	Indicates whether the figure should be shown as well.
				(default=False)
	dpi		--	The dots per inch to use for the png export, or None to use
				the default. (default=None)
	"""

	if folder != None:
		_plotFolder = os.path.join(plotFolder, folder)
	else:
		_plotFolder = plotFolder
	try:
		os.makedirs(os.path.join(_plotFolder, 'svg'))
	except:
		pass
	try:
		os.makedirs(os.path.join(_plotFolder, 'png'))
	except:
		pass
	pathSvg = os.path.join(_plotFolder, 'svg', '%s.svg' % name)
	pathPng = os.path.join(_plotFolder, 'png', '%s.png' % name)
	plt.savefig(pathSvg)
	plt.savefig(pathPng, dpi=dpi)
	if show or '--show' in sys.argv:
		plt.show()
	else:
		plt.clf()
