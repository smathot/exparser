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
import time
import numpy as np
from exparser.DataMatrix import DataMatrix

skipCache = '--no-cache' in sys.argv
cacheFolder = '.cache'
if '--clear-cache' in sys.argv and os.path.exists(cacheFolder):
	print 'Removing cache folder (%s)' % cacheFolder
	shutil.rmtree(cacheFolder)
if not os.path.exists(cacheFolder):
	print 'Creating cache folder (%s)' % cacheFolder
	os.mkdir(cacheFolder)

def cachedArray(func):

	"""
	A decorator function that provides a cache for functions that return
	numpy arrays.

	Arguments:
	func		--	A function.
	"""

	def inner(*args, **kwargs):

		if 'cacheId' in kwargs:
			cachePath = os.path.join(cacheFolder, kwargs['cacheId']) + '.npy'
			del kwargs['cacheId']
		else:
			cachePath = None
		if skipCache or cachePath == None or not os.path.exists(cachePath):
			a = func(*args, **kwargs)
			if cachePath != None:
				print '@cachedArray: saving %s' % cachePath
				np.save(cachePath, a)
		else:
			cTime = time.ctime(os.path.getctime(cachePath))
			print '@cachedArray: loading %s (created %s)' % (cachePath, cTime)
			a = np.load(cachePath)
		return a

	return inner

def cachedDataMatrix(func):

	"""
	A decorator function that provides a cache for functions that return
	DataMatrices.

	Arguments:
	func		--	A function.
	"""

	def inner(*args, **kwargs):

		if 'cacheId' in kwargs:
			cachePath = os.path.join(cacheFolder, kwargs['cacheId']) + '.npy'
			del kwargs['cacheId']
		else:
			cachePath = None
		if skipCache or cachePath == None or not os.path.exists(cachePath):
			dm = func(*args, **kwargs)
			if cachePath != None:
				print '@cachedDataMatrix: saving %s' % cachePath
				dm.save(cachePath)
		else:
			cTime = time.ctime(os.path.getctime(cachePath))
			print '@cachedDataMatrix: loading %s (created %s)' % (cachePath, \
				cTime)
			dm = DataMatrix(cachePath)
		return dm

	return inner
