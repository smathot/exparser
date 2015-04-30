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
from exparser.CsvReader import CsvReader
import os
import os.path
import numpy as np
import csv
import warnings

class CsvFolderReader(BaseReader):

	def __init__(self, path='data', ext='.csv', delimiter=',', quote='"',
		maxN=None):

		"""
		desc:
			Constructor. Reads all csv files from a specific folder.

		keywords:
			path:
				desc:	The folder containing the csv files.
				type:	[str, unicode]
			ext:
				desc:	The extension of the csv files, or None to parse all
						files.
				type:	[str, unicode, NoneType]
			delimiter:
				desc:	The column delimiter.
				type:	[str, unicode]
			quote:
				desc:	The character used for quoting columns.
				type:	[str, unicode, NoneType]				
			maxN:
				desc:	The maximum number of files to process, or `None` to
						process all.
				type:	[int, NoneType]
		"""

		print 'Scanning \'%s\'' % path
		self.dm = None
		l = sorted(os.listdir(path))
		if maxN != None:
			l = l[:maxN]
		for fname in l:
			if ext == None or os.path.splitext(fname)[1] == ext:
				print 'Reading %s ...' % fname,
				try:
					cr = CsvReader(os.path.join(path, fname),
						delimiter=delimiter, quote=quote)
				except:
					warnings.warn('Failed to read %s' % fname)
					continue
				dm = cr.dataMatrix()
				while '__src__' not in dm.columns():
					try:
						# Apparently this can go wrong sometimes, presumably
						# due to a low-level bug. It appears to be safe to just
						# try again.
						dm = dm.addField('__src__', dtype=str)
					except:
						print 'Trying again ...'
				dm['__src__'] = fname
				if self.dm == None:
					self.dm = dm
				else:
					self.dm += dm
				print '(%d rows)' % len(dm)

	def dataMatrix(self):

		"""
		returns:
			type:	DataMatrix
		"""

		return self.dm
