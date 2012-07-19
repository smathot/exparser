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

from exparser.BaseReader import BaseReader
from exparser.DataMatrix import DataMatrix
import os
import os.path
import numpy as np

class CsvFolderReader(BaseReader):

	def __init__(self, path='data', ext='.csv', delimiter=',', quote=None):

		"""
		Constructor. Reads all csv files from a specific folder.

		Keyword arguments:
		path -- the folder containing the csv files (default='data')
		ext -- the extension of the csv files (default='.csv')
		quote -- the character used for quoting columns (default=None)
		"""

		print 'Scanning \'%s\'' % path
		self.m = None
		for fname in os.listdir(path):
			if os.path.splitext(fname)[1] == ext:
				print 'Reading %s ...' % fname,
				a = np.loadtxt(os.path.join(path, fname), dtype=str, \
					delimiter=delimiter)
				if self.m == None:
					self.m = a
				else:
					self.m = np.concatenate( (self.m, a[1:]) )
				print '(%d rows)' % (a[:,0].size-1)

		if quote != None:
			self.m = np.char.lstrip(self.m, quote)
			self.m = np.char.rstrip(self.m, quote)

	def dataMatrix(self):

		"""
		Returns:
		A DataMatrix
		"""

		return DataMatrix(self.m)
