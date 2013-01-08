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
import csv

class CsvReader(BaseReader):

	def __init__(self, path='data.csv', delimiter=',', quote='"'):

		"""
		Constructor. Reads a single csv file.

		Keyword arguments:
		path -- the folder containing the csv files (default='data')
		delimiter -- the column delimiter (default=',')
		quote -- the character used for quoting columns (default='"')
		"""

		reader = csv.reader(open(path), delimiter=delimiter, quotechar=quote)
		l = []
		for row in reader:
			l.append(row)
		self.m = np.array(l, dtype='|S128')

	def dataMatrix(self):

		"""
		Returns:
		A DataMatrix
		"""

		return DataMatrix(self.m)
