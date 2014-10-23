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
from scipy.spatial.distance import cdist
from scipy import stats
from scipy.stats.stats import nanmean, nanstd, nanmedian
from copy import deepcopy
from itertools import product
from exparser.BaseMatrix import BaseMatrix
from exparser.PivotMatrix import PivotMatrix
from exparser import Constants

class DataMatrix(BaseMatrix):

	"""
	desc:
		Provides functionality for convenient processing of experimental data.
	"""

	def __init__(self, a, structured=False):

		"""
		desc: |
			Creates a new DataMatrix object.

		arguments:
			a:
				desc:	A NumPy array, list, or filename. For unstructured
						NumPy arrays or lists, the first row is assumed to
						contain column names. Filenames are assumed to refer to
						a `.npy` file.
				type:	[ndarray, list, str, unicode]


		keywords:
			structured:
				desc:	Indicates whether `a` is a structured NumPy array.
				type:	bool
		"""

		# Load from disk
		if isinstance(a, basestring):
			a = np.load(a)

		# Try to convert lists to arrays
		if isinstance(a, list):
			a = np.array(a)

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
					dtype.append( (vName, Constants.strDType) )
			try:
				self.m = np.array(a, dtype=dtype)
			except Exception as e:
				raise Exception('Failed to convert %s (original exception %s)' \
					% (a, e))
			return

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
					dtype.append( (vNames[i], Constants.strDType) )

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

	def __getitem__(self, key):

		"""
		desc:
			Returns a column, index, or slice. Note that some operations return
			a copy of the DataMatrix, so they cannot be used to modify the
			contents of the DataMatrix.

		example: |
			dm['rt'] # returns column 'rt' as numpy array
			dm[0] # returns first row as DataMatrix
			dm[0:2] # returns first two rows as DataMatrix
			dm[0]['rt'] = 1 # This doesn't alter the original DataMatrix
			dm['rt'][0] = 1 # This does!

		arguments:
			key:
				desc:	A column name or index.
				type:	[int, str, unicode]

		returns:
			desc:	A DataMatrix or NumPy array corresponding to a slice or
					column from this DataMatrix.
			type:	[DataMatrix, ndarray]
		"""

		if type(key) == int or type(key) == np.int64:
			a = np.array(self.m[key])
			a.shape = (1,)
			return DataMatrix(a, structured=True)
		elif isinstance(key, basestring):
			return self.m[key]
		elif isinstance(key, slice) or isinstance(key, list) \
			or isinstance(key, np.ndarray):
			return DataMatrix(self.m[key], structured=True)
		else:
			raise Exception('Cannot get %s (%s)' % (key, type(key)))

	def __len__(self):

		"""
		returns:
			desc:	The number of rows.
			type:	int
		"""

		return len(self.m)

	def __add__(self, dm, cautious=False):

		"""
		desc:
			Concatenates two DataMatrices. Implements the + operator.

		arguments:
			dm:
				desc:	The DataMatrix to be appended.
				type:	DataMatrix

		keywords:
			cautious:	DEPRECATED

		example: |
			dm3 = dm1 + dm2

		returns:
			desc:	The concatenation of the current and the passed DataMatrix.
			type:	DataMatrix
		"""

		cols = np.intersect1d(self.columns(), dm.columns())
		a = np.zeros( (1+len(self)+len(dm), len(cols)), dtype='|S128')
		i = 0
		for col in cols:
			a[0, i] = str(col)
			if self[col].dtype != dm[col].dtype:
				warnings.warn( \
					'%s has non-matching types (%s and %s)' \
					% (col, self[col].dtype, dm[col].dtype))
				a1 = np.array(self[col], dtype='|S128')
				a2 = np.array(dm[col], dtype='|S128')
				a[1:, i] = np.concatenate( (a1, a2) )
			else:
				a[1:, i] = np.concatenate( (self[col], dm[col]) )
			i += 1
		return DataMatrix(a)

	def __setitem__(self, key, val):

		"""
		desc:
			Set a certain variable. Implements assigment.

		arguments:
			key:
				desc:	The name of a key.
				type:	[str, unicode]
			val:
				desc:	An array with the new values, or a single new value to
						use for the entire column.

		example: |
			dm['rt'] = 100
		"""

		self.m[key] = val

	def __iter__(self):

		"""
		desc:
			Implements an iterator for 'for' loops to walk through a DataMatrix
			row by row.

		example: |
			for rowDm in dm:
				print rowDm

		returns:
			type:	DataMatrixIterator
		"""

		return DataMatrixIterator(self)

	def addField(self, key, dtype=np.int32, default=None):

		"""

		desc:
			Creates a new DataMatrix that is a copy of the current DataMatrix
			with an additional field.

		source:
			http://stackoverflow.com/questions/1201817/\
			adding-a-field-to-a-structured-numpy-array>

		arguments:
			key:
				desc:	The name of the new field.
				type:	[str, unicode]

		keywords:
			dtype:		The dtype for the new field.
			default:	The default value or `None` for no default.

		example: |
			dm = dm.addField('rt', dtype=float, default=-1000)

		returns:
			type:	DataMatrix.
		"""

		if key in self.columns():
			raise Exception('field "%s" already exists' % key)
		a = np.zeros(self.m.shape, dtype=self.m.dtype.descr + [(key, dtype)])
		for name in self.m.dtype.names:
			a[name] = self.m[name]
		dm = DataMatrix(a, structured=True)
		if default != None:
			dm[key] = default
		return dm

	def asArray(self):

		"""
		returns:
			desc:	An array representation of the current DataMatrix.
			type:	ndarray
		"""

		l = [list(self.m.dtype.names)]
		for row in self.m:
			l.append(list(row))
		return np.array(l, dtype='|S128')

	def balance(self, key, maxErr, ref=0, verbose=False):

		"""
		desc:
			Filters the data such that a given column is on average close to a
			reference value, and is symetrically distributed.

		arguments:
			key:
				desc:	The key to balance. It must have a numeric dtype.
				type:	[str, unicode]
			maxErr:
				desc:	The maximum mean error relative to the reference value.
				type:	[int, float]

		keywords:
			ref:
				desc:	The reference value.
				type:	[int, float]
			verbose:
				desc:	Indicates whether verbose output is printed.
				type:	bool

		returns:
			desc:	A balanced copy of the current DataMatrix.
			type:	DataMatrix
		"""

		dm = self.clone()

		# Create a distance matrix for the values, with the diagonal set to nan
		a = dm[key] - ref
		a.shape = len(a), 1 # cdist requires 2D array
		d = cdist(a, -a)
		np.fill_diagonal(d, np.nan)

		# Create the best matching matrix
		pairs = []
		while False in np.isnan(d):
			c, r = np.where(d == np.nanmin(d))
			i1 = c[0]
			i2 = r[0]
			if verbose:
				print '%.5d\t->\t%.5d' % (i1, i2)
			err = abs(a[i1][0]+a[i2][0])
			pairs.append( (i1, i2, err) )
			d[i1] = np.nan
			d[:,i1] = np.nan
			d[i2] = np.nan
			d[:,i2] = np.nan
		# Mark rows in a pairwise fashion until the overall error is low enough.
		toRemove = []
		_dm = dm.clone()
		while True:
			err = _dm[key].mean()
			if verbose:
				print 'Error = %f' % err
			#if abs(err) <= maxErr:
				#print 'Done!'
				#break
			if len(pairs) == 1:
				#print 'Set exhausted' (default=False)
				break
			i1, i2, err = pairs.pop()
			v1 = dm[key][i1]
			v2 = dm[key][i2]
			pErr = v2+v1

			if abs(pErr) <= maxErr:
				break

			#print 'Marking rows %d and %d for removal' % (i1, i2)
			toRemove += [i1, i2]
			_dm.m = np.delete(dm.m, toRemove)
			#print '%f - %f = %f' % (v1, v2, v2+v1)

		# Remove rows
		dm = dm.addField('__unbalanced__', dtype=int, default=0)
		dm['__unbalanced__'][toRemove] = 1
		return dm

	def calcPerc(self, vName, targetVName, keys=None, nBin=None):

		"""
		desc:
			Calculates percentile scores for a variable.

		arguments:
			vName:
				desc:	The variable to calculate percentile scores for.
				type:	[str, unicode]
			targetVName:
				desc:	The variable to store the percentile scores in. This
						variable must exist, it is not created.
				type:	[str, unicode]

		keywords:
			keys:
				desc:	A key or list of keys to split the data by, before
						calculating percentile scores, so you can calculate
						scores individually per subject, condition, etc.
				type:	[list, str, unicode]
			nBin:
				desc:	The number of bins or None for continuous scores.
				type:	[int, NoneType]

		returns:
			type:	DataMatrix
		"""

		print "vName = ", vName
		print "targetVName = ", targetVName

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
		desc:
			Collapse the data by a (list of) keys and get statistics on a
			dependent variable.

		arguments:
			keys:
				desc:	A key or list of keys to collapse the data on.
				type:	[list, str, unicode]
			vName:
				desc:	The dependent variable to collapse. Alternative, you can
						specifiy a function, in which case the error will be 0.
				type:	[str, unicode, function]

		returns:
			desc:	A DataMatrix with the collapsed data, with the descriptives
					statistics on `vName`.
			type:	DataMatrix
		"""

		if isinstance(keys, basestring):
			keys = [keys]

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
				l.append(nanmean(a))
				l.append(nanmedian(a))
				l.append(nanstd(a))
				l.append(nanstd(a)/np.sqrt(a.size))
				l.append(1.96*nanstd(a)/np.sqrt(a.size))
				l.append(a.size)
			m.append(l)
		return DataMatrix(m)

	def columns(self, dtype=False):

		"""
		desc:
			Returns a list of the columns.

		keywords:
			dtype:
				desc:	Indicates whether the dtype for each column should be
						returned as well.
				type:	bool

		returns:
			desc: |
				If dtype == False: A list of names
				If dtype == True: A list of (name, dtype) tuples
			type:	list
		"""

		if dtype:
			return self.m.dtype.descr
		return list(self.m.dtype.names)

	def count(self, key):

		"""
		desc:
			Returns the number of different values for a given variable.

		arguments:
			key:
				desc:	The variable to count.
				type:	[str, unicode]

		returns:
			desc:	The number of different values for 'key'.
			type:	int
		"""

		return len(self.unique(dv))

	def explode(self, N):

		"""
		desc:
			Break up the DataMatrix in N smaller DataMatrices. For splitting a
			DataMatrix based on column values, see [DataMatrix.split].

		arguments:
			N:
				desc:	The number of DataMatrices to explode in.
				type:	int

		returns:
			desc:	A list of DataMatrices.
			type:	list
		"""

		l = []
		for a in np.array_split(self.m, N):
			dm = DataMatrix(a, structured=True)
			l.append(dm)
		return l

	def intertrialer(self, keys, dv, _range=[1]):

		"""
		desc:
			Adds columns that contain values from the previous or next trial.
			These columns are called '[dv]_p1' for the next value, '[dv]_m1' for
			the previous one, etc.

		arguments:
			keys:
				desc:	A key or list of keys that define the trial order.
				type:	[list, str, unicode]
			dv:
				desc:	The dependent variable.
				type:	[str, unicode]

		keywords:
			_range:
				desc:	A list of integers that specifies the range for which
						the operation should be executed.
				type:	list

		returns:
			type:	DataMatrix
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

	def range(self):

		"""
		desc:
			Gives a list of indices to walk through the current DataMatrix.

		returns:
			A list of indices.
		"""

		return range(len(self))

	def recode(self, key, coding):

		"""
		desc:
			Recodes values (i.e. changes one value into another for a given set
			of columns).

		arguments:
			key:
				desc:	The name of the variable to recode, or a list of names
						to recode multiple variables in one go.
				type:	[str, unicode, list]
			coding:
				desc:	An (oldValue, newValue) tuple, a list of tuples to
						handle multiple recodings in one go, or a function that
						takes a value and returns the recoded value.
				type:	[tuple, list, function]
		"""

		if type(key) == list:
			for _key in key:
				self.recode(_key, coding)
			return

		if type(coding) == list:
			for _coding in coding:
				self.recode(key, _coding)
			return
		elif hasattr(coding, '__call__'):
			for i in range(len(self.m)):
				self.m[key][i] = coding(self.m[key][i])
		else:
			oldValue, newValue = coding
			i = np.where(self.m[key] == oldValue)
			self.m[key][i] = newValue

	def removeNan(self, key):

		"""
		desc:
			Remove all rows where the specified key is nan.

		arguments:
			key:
				desc:	A key that should not have any nan values.
				type:	[str, unicode]

		returns:
			type:	DataMatrix
		"""

		i = np.where(~np.isnan(self.m[key]))[0]
		print 'Removing %d rows with nans' % len(i)
		dm = DataMatrix(self.m[i], structured=True)
		return dm

	def removeField(self, key):

		"""
		desc:
			Return a DataMatrix that is a copy of the current DataMatrix without
			the specified field.

		arguments:
			key:
				desc:	The name of the field to be removed.
				type:	[str, unicode]

		returns:
			type:	DataMatrix
		"""

		newDtype = []
		for i in range(len(self.m.dtype)):
			_key = self.m.dtype.names[i]
			_dtype = self.m.dtype[i]
			if _key != key:
				newDtype.append( (_key, _dtype) )
		a = np.zeros(self.m.shape, dtype=newDtype)
		for name in self.m.dtype.names:
			if name != key:
				a[name] = self.m[name]
		return DataMatrix(a, structured=True)

	def rename(self, oldKey, newKey):

		"""
		desc:
			Renames a column. This function operates in place, so it modifies
			the current dataMatrix.

		arguments:
			oldKey:
				desc:	The old name.
				type:	[str, unicode]
			newKey:
				desc:	The new name.
				type:	[str, unicode]

		returns:
			type:	DataMatrix
		"""

		dtype = []
		for key in self.m.dtype.names:
			_type = self.m.dtype[key]
			if key == oldKey:
				key = newKey
			dtype.append( (key, _type) )
		self.m.dtype = dtype
		return self

	def group(self, keys, _sort=True):

		"""
		desc:
			Split the data into different groups based on unique values for the
			key variables.

		arguments:
			keys:
				desc:	A key or list of keys to split the data on.
				type:	[str, unicode, list]

		keywords:
			"_sort":
				desc:	Indicates whether the groups should be sorted by values.
				type:	bool

		returns:
			desc:	A list of DataMatrices.
			type:	list
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
		desc:
			Select a subset of the data.

		arguments:
			query:
				desc:	A query, e.g. 'rt > 1000'.
				type:	[str, unicode]

		keywords:
			verbose:
				desc:	Indicates if a summary should be printed.
				type:	bool

		returns:
			type:	DataMatrix
		"""

		if len(self) == 0:
			warnings.warn('Selecting from empty DataMatrix (%s)' % query)
			return self.clone()
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

	def selectColumns(self, keys):

		"""
		desc:
			Creates a new DataMatrix with only the specified columns.

		arguments:
			keys:
				desc:	A column or list of columns to select.
				type:	[list, str, unicode]

		returns:
			type:	DataMatrix
		"""

		if isinstance(keys, basestring):
			keys = [keys]
		for key in keys:
			if key not in self.columns():
				raise Exception('The column "%s" does not exist' % key)
		return DataMatrix(self.m[keys], structured=True)

	def selectByStdDev(self, keys, dv, thr=2.5, verbose=False):

		"""
		desc:
			Select only those rows where the value of a given column is within a
			certain distance from the mean.

		arguments:
			keys:
				desc:	A list of keys to create groups for which the
						deviation is calculated seperately.
				type:	list
			dv:
				desc:	The dependent variable.
				type:	[str, unicode]

		keywords:
			thr:
				desc:	The stddev threshold.
				type:	[float, int]
			verbose:
				desc:	Indicates whether detailed output should be provided.
				type:	bool

		returns:
			type:	DataMatrix
		"""
		if verbose:
			print '======================================='
			print '| Selecting based on Standard Deviation'
			print '| Threshold: %.2f' % thr
			print '| Keys: %s' % ','.join(keys)
			print '|'
		dm = self.clone()

		# Create a dummy field that combines the keys, so we can simply group
		# based on one key
		dm = dm.removeField('__dummyCond__')
		dm = dm.addField('__dummyCond__', dtype=str, default='__dummy__')
		for key in keys:
			for i in range(len(dm)):
				dm['__dummyCond__'][i] += str(dm[key][i]) + '__'

		# Add a field to store outliers
		dm = dm.removeField('__stdOutlier__')
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

		dm = dm.select('__stdOutlier__ == 0', verbose = verbose)
		print
		return dm

	def shuffle(self):

		"""
		desc:
			Shuffles the DataMatrix in place.
		"""

		# Directly shuffling the array does not preserve all items! This seems
		# to be a bug in numpy. This workaround preserves the integrity of the
		# DataMatrix.
		l = list(self.m)
		np.random.shuffle(l)
		self.m = np.array(l, dtype=self.m.dtype)

	def sort(self, keys, ascending=True):

		"""
		desc:
			Sorts the DataMatrix in place.

		arguments:
			keys:
				desc:	A key or list of keys to use for sorting. The first key
						is dominant, the second key is next-to-dominant, etc.
				type:	[str, unicode, list]

		keywords:
			ascending:
				desc:	Indicates whether the sorting should occur in ascending
						(True) or descending (False) order.
				type:	bool
		"""

		self.m.sort(order=keys)
		if not ascending:
			self.m = self.m[::-1]

	def split(self, key):

		"""
		desc:
			Splits the DataMatrix in chunks such that each chunk only has the
			same value for the specified column. For splitting a DataMatrix into
			equally sized parts, see [DataMatrix.explode].

		arguments:
			key:
				desc:	A key to split by.
				type:	[str, unicode]

		returns:
			desc:	A list of DataMatrices.
			type:	list
		"""

		l = []
		start_index = 0
		i = 0
		val = self[col][0] # Starting value
		while i < len(self)-1:
			i += 1
			_val = self[col][i]
			if _val != val:
				l.append(self[start_index:i])
				start_index = i
				val = _val
		l.append(self[start_index:i+1])
		return l

	def ttest(self, keys, dv, paired=True, collapse=None):

		"""
		desc:
			Performs t-tests between groups defined by a list of keys.

		arguments:
			keys:
				desc:	A list of keys to define the groups.
				type:	list
			dv:
				desc:	The dependent variable.
				type:	[str, unicode]

		keywords:
			paired:
				desc:	Determines whether a paired-samples t-test or an
						independent samples t-test should be conducted.
				type:	bool
			collapse:
				desc:	A key to collapse the data on, so that you can do
						t-tests on (subject) means.
				type:	[str, unicode, NoneType]

		returns:
			desc:	A list of (desc, t, p) tuples.
			type:	list
		"""

		from itertools import combinations
		if paired:
			from scipy.stats import ttest_rel as ttest
		else:
			from scipy.stats import ttest_ind as ttest

		if collapse != None:
			dm = self.collapse(collapse + keys, dv)
			dv = 'mean'
		else:
			dm = self

		_l = [['group', 'N', 'M / t', 'SE / p']]

		lDm = dm.group(keys)
		for l in combinations(lDm, 2):

			group0 = ''
			for key in keys:
				group0 += str(l[0][key][0]) + '_'
			group0 = group0[:-1]

			group1 = ''
			for key in keys:
				group1 += str(l[1][key][0]) + '_'
			group1 = group1[:-1]

			N0 = len(l[0])
			M0 = l[0][dv].mean()
			SE0 = l[0][dv].std() / np.sqrt(len(l[0]))
			_l.append( [group0, N0, M0, SE0] )

			N1 = len(l[1])
			M1 = l[1][dv].mean()
			SE1 = l[1][dv].std() / np.sqrt(len(l[1]))
			_l.append( [group1, N1, M1, SE1] )

			t, p = ttest(l[0][dv], l[1][dv])
			_l.append( [group0, group1, t, p] )

		return DataMatrix(np.array(_l))

	def unique(self, key):

		"""
		desc:
			Gives all unique values for a particular key.

		arguments:
			key:
				desc:	A column name.
				type:	[str, unicode]

		returns:
			desc:	A list of unique values.
			type:	list
		"""

		return list(np.unique(self[dv]))

	def where(self, query):

		"""
		desc:
			Return indices corresponding to the query.

		arguments:
			query:
				desc:	A query, e.g. 'rt > 1000'.
				type:	[str, unicode]

		returns:
			desc:	Indices.
			type:	ndarray
		"""

		l = query.split(' ')
		vName = l[0]
		op = l[1]
		test = query[len(vName)+len(op)+1:]
		flt = np.where(eval('self.m[vName] %s %s' % (op, test)))
		return flt[0]

	def withinize(self, vName, targetVName, key, verbose=True, whiten=False):

		"""
		desc:
			Removes the between factor variance for a given key (such as subject
			or file) for a given dependent variable. This operation acts in
			place.

		arguments:
			vName:
				desc:	The dependent variable to withinize.
				type:	[str, unicode]
			targetVName:
				desc:	The target variable to store withinized values. This
						variable should exist.
				type:	[str, unicode]
			key:
				desc:	The key that defines the group.
				type:	[str, unicode]

		keywords:
			verbose:
				desc:	Indicates whether the results should be printed.
				type:	bool
			whiten:
				desc:	Indicates whether the data should be whitened so that
						the standard deviation is 1 and the mean 0.
				type:	bool

		returns:
			type:	DataMatrix
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

		if self.i >= len(self.dm):
			raise StopIteration
		a = self.dm.m[[self.i]]
		self.i += 1
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
