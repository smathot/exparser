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

import numpy
import copy

class BaseMatrix:

	"""Base class for the *Matrix classes"""

	def __init__(self, m):

		"""
		Constructor

		Arguments:
		m -- a numpy array
		"""

		self.m = m

	def __str__(self):

		"""
		Returns:
		A string representation
		"""

		return '%s' % self.asArray()

	def _print(self, title=None, sign=2, maxLen=20, hSep='=', vSep='|'):

		"""
		Prints a pivot matrix to the standard output

		Keyword arguments:
		title -- a title to be printed above the table (default=None)
		m -- the pivot matrix. If none, the current DataMatrix is printed
			 (default=None)
		sign -- the maximum number of decimals (default=2)
		maxLen -- the maximum column length (default=20)
		hSep -- the horizontal separator character (default='+')
		vSep -- the vertical separator character (default='|')
		"""

		m = self.asArray()
		l = 0
		for row in m:
			for col in row:
				try:
					s = str(('%%.%df' % sign) % float(col))
				except:
					s = str(col)
				l = min(max(l, len(str(s))), maxLen)

		totalL = (1+(l+2)*row.size)*len(vSep)

		print hSep*totalL
		if title != None:
			j = (totalL-len(title))/2-1
			print vSep + (' '*j + title + ' '*j).ljust(totalL-2) + vSep
			print hSep*totalL

		for row in m:
			_l = []
			for col in row:
				try:
					_l.append(str(('%%.%df' % sign) % float(col)).rjust(l+1))
				except:
					_l.append(str(col)[:maxLen].ljust(l+1))
			print '%s%s%s' % (vSep, vSep.join(_l), vSep)
		print hSep*totalL

	def asArray(self):

		"""
		Returns:
		An array representation
		"""

		return self.m

	def clone(self):

		"""
		Creates a clone (deep copy) of the current matrix

		Returns:
		A BaseMatrix
		"""

		return copy.deepcopy(self)

	def save(self, path='DataMatrix.csv', delimiter=','):

		"""
		Write the current DataMatrix to a plaintext csv file

		Keyword arguments:
		path -- the path of the file to write (default='DataMatrix.csv')
		delimiter -- the character used to separate columns (default=',')
		"""

		numpy.savetxt(path, self.asArray(), fmt='%s', delimiter=delimiter)
