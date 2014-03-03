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
from matplotlib import pyplot as plt

plotFolder = 'plot'
if '--clear-plot' in sys.argv and os.path.exists(cacheFolder):
	print 'Removing plot folder (%s)' % plotFolder
	shutil.rmtree(plotFolder)
plt.rc('font', family='Arial', size=10)

r = 8, 8
l = 12, 12
xl = 16, 16

def new(size=r):

	"""
	Creates a new figure.

	Keyword arguments:
	size	--	The figure size. (default=(12,8))

	Returns:
	A matplotlib figure.
	"""

	return plt.figure(figsize=size)

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
