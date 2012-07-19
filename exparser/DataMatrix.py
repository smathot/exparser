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

import sys
import warnings
import numpy as np
from numpy import ma
from numpy.lib import recfunctions
from scipy import stats
from copy import deepcopy
from itertools import product
from exparser.BaseMatrix import BaseMatrix
from exparser.PivotMatrix import PivotMatrix

class DataMatrix(BaseMatrix):

	"""Provides functionality for convenient processing of experimental data"""

	def __init__(self, a, structured=False):

		"""
		Constructor. Creates a DataMatrix from a np array. If this array is
		not stuctured, it is assumed to have column names on the first row and
		values in the other rows, and to be of a string type

		Arguments:
		a -- a np array or list

		Keyword arguments:
		structured -- indicates whether the passed array is structured. If not,
					  a structured array will be created, assuming that the
					  first row contains the column names. (default=False)
		"""

		# If a structured array was passed, we still need to convert it, to
		# choose the optimal data type for each column.
		if structured:
			dtype = []
			for i in range(len(a.dtype.names)):
				vName = a.dtype.names[i]
				try:
					np.array(a[vName], dtype=np.float64)
					dtype.append( (vName, np.float64) )
				except:
					dtype.append( (vName, '|S128') )
			self.m = np.array(a, dtype=dtype)
			return

		# Try to convert lists to arrays
		if type(a) == list:
			a = np.array(a)

		# Extract the variable names (first row) and values (the rest)
		vNames = a[0]
		vVals = a[1:]
		_vVals = vVals.swapaxes(0,1)

		# Determine the appropriate dtype for each column
		dtype = []
		for i in range(len(vNames)):
			try:
				np.array(_vVals[i], dtype=np.int32)
				dtype.append( (vNames[i], np.int32) )
			except:
				try:
					np.array(_vVals[i], dtype=np.float64)
					dtype.append( (vNames[i], np.float64) )
				except:
					dtype.append( (vNames[i], '|S128') )

		# Fill the structured array row by row
		self.m = np.zeros( len(vVals), dtype=dtype)
		for i in range(len(vVals)):
			try:
				self.m[i] = tuple(vVals[i])
			except:
				for j in range(len(vVals[i])):
					try:
						self.m[i][j] = vVals[i][j]
					except:
						self.m[i][j] = 0
						warnings.warn( \
							"Column %d is not numeric but '%s'. Falling back to 0" \
							% (j, vVals[i][j]))

	def __getitem__(self, vName):

		"""
		Return a column

		Arguments:
		vName -- the name of the variable
		"""

		return self.m[vName]

	def __len__(self):

		"""
		Returns:
		The number of rows
		"""

		return len(self.m)

	def __setitem__(self, vName, vVal):

		"""
		Set a certain variable

		Arguments:
		vName -- the name of the variable
		vVal -- an array with the new values
		"""

		self.m[vName] = vVal

	def addField(self, vName, dtype=np.int32):

		"""
		Modified from: <http://stackoverflow.com/questions/1201817/\
			adding-a-field-to-a-structured-numpy-array>

		Return a new array that is like "a", but has additional fields. The
		contents of "a" are copied over to the appropriate fields in the new
		array, whereas the new fields are uninitialized.  The arguments are not
		modified.

		Arguments:
		descr -- a numpy type description of the new fields

		Returns:
		A DataMatrix
		"""

		a = np.zeros(self.m.shape, dtype=self.m.dtype.descr + [(vName, dtype)])
		for name in self.m.dtype.names:
			a[name] = self.m[name]
		return DataMatrix(a, structured=True)

	def asArray(self):

		"""
		Returns:
		An array representation
		"""

		l = [list(self.m.dtype.names)]
		for row in self.m:
			l.append(list(row))
		return np.array(l, dtype='|S128')

	def calcPerc(self, vName, targetVName, keys=None, nBin=None):

		"""
		Calculates percentile scores for a variable

		Arguments:
		vName -- the variable to calculate percentile scores for
		targetVName -- the variable to store the percentile scores in. This
					   variable must exist, it is not created.

		Keyword arguments:
		keys -- keys to split the data by, before calculating percentile
				scores, so you can calculate scores individually per subject,
				condition, etc. (default=None)
		nBin -- the number of bins or None for continues scores (default=None)

		Returns:
		A DataMatrix
		"""

		if keys != None:
			l = [list(self.m.dtype.names)]
			for dm in self.group(keys):
				dm = dm.calcPerc(vName, targetVName, nBin=nBin)
				for row in dm.m:
					l.append(list(row))
			return DataMatrix(np.array(l))

		dm = DataMatrix(self.m, structured=True)
		for row in dm.m:
			p = stats.percentileofscore(self.m[vName], row[vName])
			if nBin != None:
				binSize = 100./nBin
				p = int(p-p%(binSize+.00000000001))
			row[targetVName] = p
		return dm

	def collapse(self, keys, vName):

		"""
		Collapse the data by a (list of) keys and get statistics on a dependent
		variable.

		Arguments:
		keys -- a list of key names

		Returns:
		A DataMatrix with the collapsed data, with the following descriptives on
		the vName variable
		"""

		m = [keys + ['mean', 'median', 'std', 'se', '95ci', 'count']]
		for g in self.group(keys):
			l = []
			for key in keys:
				l.append(g[key][0])
			a = g[vName]
			l.append(a.mean())
			l.append(np.median(a))
			l.append(a.std())
			l.append(a.std()/np.sqrt(a.size))
			l.append(1.96*a.std()/np.sqrt(a.size))
			l.append(a.size)
			m.append(l)
		return DataMatrix(m)

	def group(self, keys):

		"""
		Split the data into different groups based on unique values for the
		key variables.

		Arguments:
		keys -- a list of variable names

		Returns:
		A list of DataMatrices
		"""

		if len(keys) == 0:
			return [self]
		vName = keys[0]
		vVals = np.unique(self.m[vName])
		l = []
		for vVal in vVals:
			a = self.m[self.m[vName] == vVal]
			dm = DataMatrix(a, structured=True)
			l += dm.group(keys[1:])
		return l

	def select(self, query, verbose=True):

		"""
		Select a subset of the data

		Arguments:
		query -- a query, e.g. 'rt > 1000'

		Keyword arguments:
		verbose -- indicates if a summary should be printed (default=True)

		Returns:
		A DataMatrix
		"""

		l = query.split(' ')
		vName = l[0]
		op = l[1]
		test = query[len(vName)+len(op)+1:]
		flt = eval('self.m[vName] %s %s' % (op, test))
		if type(flt) == bool: # The whole array is selected
			dm = self
		else:
			dm = DataMatrix(self.m[flt], structured=True)
		if verbose:
			a = np.empty((4,3), dtype='|S128')
			a[0,0] = '[QUERY]'
			a[0,1] = query
			a[0,2] = '-'
			a[1,0] = '[N BEFORE]'
			a[1,1] = self.m.size #len(self.m)
			a[1,2] = 100.
			a[2,0] = '[N DISCARDED]'
			a[2,1] = self.m.size -dm.m.size #len(self.m) - len(dm.m)
			a[2,2] = 100.*(self.m.size-dm.m.size)/self.m.size #100.*(len(self.m) - len(dm.m))/len(self.m)
			a[3,0] = '[N AFTER]'
			a[3,1] = dm.m.size #len(dm.m)
			a[3,2] = 100.*dm.m.size/self.m.size #100.*len(dm.m)/len(self.m)
			BaseMatrix(a)._print()
		return dm
		
	def withinize(self, dv, key, verbose=True):
	
		"""
		Removes the between factor variance for a given key (such as subject or
		file) for a given depedent variable.
		
		Arguments:
		dv -- the dependent variable to withinize
		key -- the key that defines the within group
		
		Keyword arguments:
		verbose -- indicates whether the results should be printed (default=
				   True)		
		"""
		
		gAvg = self.m[dv].mean()
		if verbose:
			print "Grand avg =", gAvg
		for f in np.unique(self.m[key]):				
			i = np.where(self.m[key] == f)
			fAvg = self.m[dv][i].mean()
			if verbose:
				print "Avg(%s) = %f" % (f, fAvg)
			self.m[dv][i] += gAvg - fAvg
		return self				

def fromMySQL(query, user, passwd, db, charset='utf8', use_unicode=True):

	"""
	Returns the results from a MySQL query as a DataMatrix. A connection is
	automatically created and closed.

	Arguments:
	query -- the SQL query
	user -- the MySQL user
	passwd -- the MySQL password
	db -- the MySQL database

	Keyword arguments:
	charset -- the character set to be used for MySQL (default='utf8')
	use_unicode -- indicated whether Unicode strings should be used
				   (default=True)

	Returns:
	A DataMatrix
	"""

	import MySQLdb
	from MySQLdb.cursors import DictCursor

	db = db = MySQLdb.connect(user=user, passwd=passwd, db=db, charset=charset,
		use_unicode=use_unicode, cursorclass=DictCursor)
	cur = db.cursor()
	cur.execute(query)
	l = []
	for row in cur.fetchall():
		if len(l) == 0:
			l.append(list(row.keys()))
		l.append(list(row.values()))
	db.close()
	return DataMatrix(np.array(l, dtype='|S256'))
