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
import types
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
				for v in _vVals[i]: int(v)
				dtype.append( (vNames[i], np.int32) )
			except:
				try:
					for v in _vVals[i]: float(v)
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

	def __add__(self, dm):

		"""
		Concatenates two DataMatrices. Implements the + operator.

		Arguments:
		dm -- the DataMatrix to be appended

		Returns:
		The concatenation of the current and the passed DataMatrix
		"""

		if self.columns() != dm.columns():
			# Allow concatenation of non-matching datamatrices, but want the
			# user to avoid trouble
			warnings.warn( \
				'Adding two DataMatrices that do not have matching columns. All non-matching columns will be ignored.')
			cols = np.intersect1d(self.columns(), dm.columns())
			a = np.zeros( (1+len(self)+len(dm), len(cols)), dtype='|S128')
			i = 0
			for col in cols:
				a[0, i] = str(col)
				a[1:, i] = np.concatenate( (self[col], dm[col]) )
				i += 1
			_dm = DataMatrix(a)
		else:
			# The easy, more efficient way of concatenating
			_dm = self.clone()
			_dm.m = np.concatenate( (self.m, dm.m) )
		return _dm

	def __setitem__(self, vName, vVal):

		"""
		Set a certain variable

		Arguments:
		vName -- the name of the variable
		vVal -- an array with the new values
		"""

		self.m[vName] = vVal

	def __iter__(self):

		"""Implements an iterator for 'for' loops"""

		return DataMatrixIterator(self)

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

		if vName in self.columns():
			raise Exception('field "%s" already exists' % vName)
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
				try:
					p = int(p-p%(binSize+.00000000001))
				except:
					p = np.nan
			row[targetVName] = p
		return dm

	def collapse(self, keys, vName):

		"""
		Collapse the data by a (list of) keys and get statistics on a dependent
		variable.

		Arguments:
		keys -- a list of key names
		vName -- the dependent variable to collapse. Alternative, you can
				 specifiy a function, in which case the error will be 0.

		Returns:
		A DataMatrix with the collapsed data, with the following descriptives on
		the vName variable
		"""

		m = [keys + ['mean', 'median', 'std', 'se', '95ci', 'count']]
		for g in self.group(keys):
			l = []
			for key in keys:
				l.append(g[key][0])
			if type(vName) == types.FunctionType:
				l.append(vName(g))
				l.append(np.nan)
				l.append(np.nan)
				l.append(np.nan)
				l.append(np.nan)
				l.append(len(g))
			else:
				a = g[vName]
				l.append(a.mean())
				l.append(np.median(a))
				l.append(a.std())
				l.append(a.std()/np.sqrt(a.size))
				l.append(1.96*a.std()/np.sqrt(a.size))
				l.append(a.size)
			m.append(l)
		return DataMatrix(m)

	def columns(self, dtype=False):

		"""
		Returns a description of the columns

		Keyword arguments:
		dtype -- indicates if the datatype for each column should be returned as
				 well (default=False)

		Returns:
		If dtype == False: A list of names
		If dtype == True: A list of (name, dtype) tuples
		"""

		if dtype:
			return self.m.dtype.descr
		return list(self.m.dtype.names)

	def intertrialer(self, keys, dv, _range=[1]):

		"""
		Adds columns that contain values from the previous or next trial. These
		columns are called '[dv]_p1' for the next value, '[dv]_m1' for the
		previous one, etc.

		Arguments:
		keys -- a list of keys that define the trial order.
		dv -- the dependent variable

		Keyword argument:
		_range -- a list of integers that specifies the range for which the
				  operation should be executed (default=[1])

		Returns:
		A new DataMatrix
		"""

		dm = self.clone()
		dm.sort(keys)
		for i in _range:
			if i == 0:
				continue
			if i < 0:
				v = '%s_m%d' % (dv, -1*i)
			else:
				v = '%s_p%d' % (dv, i)
			dm = dm.addField(v, dtype=type(dm[dv][0]))
			if i < 0:
				dm[v][:i] = dm[dv][-i:]
			else:
				dm[v][i:] = dm[dv][:-i]
		return dm

	def recode(self, key, coding):

		"""
		Recodes values (i.e. changes one value into another for a given set of
		columns).

		Arguments:
		key -- the name of the variable to recode, or a list of names to recode
			   multiple variables in one go
		coding -- an (oldValue, newValue) tuple, or a list of tuples to handle
				  multiple recodings in one go
		"""

		if type(key) == list:
			for _key in key:
				self.recode(_key, coding)
			return

		if type(coding) == list:
			for _coding in coding:
				self.recode(key, _coding)
			return

		oldValue, newValue = coding
		i = np.where(self.m[key] == oldValue)
		self.m[key][i] = newValue

	def rename(self, oldKey, newKey):

		"""
		Renames a column. This function operates in place, so it modifies the
		current dataMatrix.

		Arguments:
		oldKey	--	The old name of the column
		newKey	--	The new name of the column
		"""

		dtype = []
		for key, _type in self.m.dtype:
			if key == oldKey:
				key = newKey
			dtype.append( (key, _type) )
		self.m.dtype = dtype
		return self

	def group(self, keys, _sort=True):

		"""
		Split the data into different groups based on unique values for the
		key variables.

		Arguments:
		keys -- a list of variable names, or a single variable name

		Keyword arguments:
		_sort -- indicates whether the groups should be sorted by values
				 (default=True)

		Returns:
		A list of DataMatrices
		"""

		if type(keys) == str:
			keys = [keys]
		if len(keys) == 0:
			return [self]
		vName = keys[0]
		vVals = np.unique(self.m[vName])
		if _sort:
			vVals = sorted(list(vVals))
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

	def selectByStdDev(self, keys, dv, thr=2.5, verbose=False):

		"""
		Select only those rows where the value of a given column is within a
		certain distance from the mean

		Arguments:
		keys -- a list of column names
		dv -- the dependent variable

		Keyword arguments:
		thr -- the stddev threshold (default=2.5)
		verbose -- indicates whether detailed output should be provided
				   (default=False)

		Returns:
		A selection of the current DataMatrix
		"""

		print '======================================='
		print '| Selecting based on Standard Deviation'
		print '| Threshold: %.2f' % thr
		print '| Keys: %s' % ','.join(keys)
		print '|'
		dm = self.clone()

		# Create a dummy field that combines the keys, so we can simply group
		# based on one key
		dm = dm.addField('__dummyCond__', dtype=str)
		for key in keys:
			for i in range(len(dm)):
				dm['__dummyCond__'][i] += str(dm[key][i]) + '__'

		# Add a field to store outliers
		dm = dm.addField('__stdOutlier__', dtype=int)
		dm['__stdOutlier__'] = 0

		for cond in np.unique(dm['__dummyCond__']):
			if verbose:
				print '| Condition (recoded):', cond
			# Get mean and standard deviation of one particular condition
			iCond = np.where(dm['__dummyCond__'] == cond)[0] # Indices of trials in condition
			m = dm[dv][iCond].mean()
			sd = dm[dv][iCond].std()
			if verbose:
				print '| M = %f, SD = %f' % (m, sd)
			# Get indices of trials that are too fast or too slow given the boundaries based
			# on this particular condition
			iTooSlow = np.where(dm[dv] > m+thr*sd)[0]
			iTooFast = np.where(dm[dv] < m-thr*sd)[0]
			# Get the indices of trials only in this condition, that are too fast or too
			# slow: 'cond and (tooslow or toofast)'
			i = np.intersect1d(iCond, np.concatenate( (iTooSlow, iTooFast) ))
			if verbose:
				print '| # outliers: %d of %d' % (len(i), len(iCond))
			dm['__stdOutlier__'][i] = 1

		dm = dm.select('__stdOutlier__ == 0')
		print
		return dm

	def sort(self, keys):

		"""
		Sorts the matrix

		Arguments:
		keys -- a list of keys to use for sorting. The first key is dominant,
				the second key is next-to-dominant, etc. A single string can
				also be specified.
		"""

		self.m.sort(order=keys)

	def unique(self, dv):

		"""
		Gives all unique values for a particular variable

		Arguments:
		dv -- the dependent variable

		Returns:
		An array of unique values for dv
		"""

		return np.unique(self[dv])

	def withinize(self, vName, targetVName, key, verbose=True, whiten=False):

		"""
		Removes the between factor variance for a given key (such as subject or
		file) for a given dependent variable.

		Arguments:
		vName -- the dependent variable to withinize
		targetVName -- the target variable
		key -- the key that defines the within group

		Keyword arguments:
		verbose -- indicates whether the results should be printed (default=
				   True)
		whiten -- indicates whether the data should be whitened so that the
				  standard deviation is 1 and the mean 0 (default=False)
		"""

		self.m[targetVName] = self.m[vName]
		gAvg = self.m[targetVName].mean()
		if verbose:
			print "Grand avg =", gAvg
		for f in np.unique(self.m[key]):
			i = np.where(self.m[key] == f)
			fAvg = self.m[targetVName][i].mean()
			fStd = self.m[targetVName][i].std()
			if verbose:
				print "Avg(%s) = %f" % (f, fAvg)
			if whiten:
				self.m[targetVName][i] -= fAvg
				self.m[targetVName][i] /= fStd
			else:
				self.m[targetVName][i] += gAvg - fAvg
		return self

class DataMatrixIterator(object):

	"""Implements an iterator for the dataMatrix"""

	def __init__(self, dm):

		"""
		Constructor

		Arguments:
		dm -- a DataMatrix
		"""

		self.dm = dm
		self.i = 0

	def __iter__(self):

		"""
		Returns:
		self
		"""

		return self

	def next(self):

		"""
		Return the current row and advance by one step

		Returns:
		A DataMatrix with only the current row

		Raises:
		A StopIteration exception when the DataMatrix is finished
		"""

		if self.i >= len(self.dm)-1:
			raise StopIteration
		self.i += 1
		a = self.dm.m[[self.i]]
		return DataMatrix(a, structured=True)


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
