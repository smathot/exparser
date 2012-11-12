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
import numpy
from scipy.spatial.distance import euclidean
import warnings

class EyelinkAscFolderReader(BaseReader):

	"""Parses Eyelink ASCII data"""
	
	def __init__(self, path='data', ext='.asc', ignoreBlinks=True, \
		startTrialKey='start_trial', endTrialKey='stop_trial', \
		variableKey='var', dtype='|S128', maxN=None, requireEndTrial=True):

		"""
		Constructor. Reads all Eyelink ASCII files from a specific folder.

		Keyword arguments:
		path -- the folder containing the csv files (default='data')
		ext -- the extension of the csv files (default='.csv')
		ignoreBlinks -- indicates whether events during a blink should be
						ignored (default=True)
		startTrialKey -- the start trial keyword (default='start_trial')
		endTrialKey -- the stop trial keyword (default='stop_trial')
		variableKey -- the variable keyword (default='var')
		dtype -- the numpy dtype to be used (default='|S128')
		maxN -- the maximum number of subjects to process (default=None)
		requireEndTrial -- indicates whether an exception should be raised if a
						  trial hasn't been neatly closed. Otherwise the trial
						  is simply disregarded. (default=True)
		"""

		self.ignoreBlinks=True
		self.startTrialKey = startTrialKey
		self.endTrialKey = endTrialKey
		self.variableKey = variableKey
		self.dtype = dtype
		self.requireEndTrial = requireEndTrial

		print '\nScanning \'%s\'' % path
		self.m = None
		nFile = 0		
		for fname in os.listdir(path):
			if os.path.splitext(fname)[1] == ext:
				print 'Reading %s ...' % fname,
				a = self.parseFile(os.path.join(path, fname))
				if self.m == None:
					self.m = a
				else:
					try:
						self.m = numpy.concatenate( (self.m, a[1:]), axis=0)
					except ValueError as e:
						print
						print 'Trying to concatenate'
						print a[0]
						print a[1:][-1]
						print 'to'
						print self.m[0]
						print self.m[-1]
						raise(e)
				print '(%d trials)' % (a[:,0].size-1)
				nFile += 1
			if maxN != None and nFile >= maxN:
				break
		print '%d files\n' % nFile

	def startBlink(self, l):

		"""
		Detects a blink start

		Returns:
		True or False
		"""

		# SBLINK R 4999128
		return l[0] == 'SBLINK'

	def endBlink(self, l):

		"""
		Detects a blink end

		Returns:
		True or False
		"""

		# EBLINK R 4999128
		return l[0] == 'EBLINK'

	def dataMatrix(self):

		"""
		Returns:
		A DataMatrix
		"""

		return DataMatrix(self.m)

	def startTrial(self, l):

		"""
		Determines whether a list corresponds to the start of a trial

		Returns:
		The trialid or None if the list is not a trial start
		"""

		if len(l) > 3 and l[0] == 'MSG' and l[2] == self.startTrialKey:
			return l[3]
		return None

	def endTrial(self, l):

		"""
		Determines whether a list corresponds to the end of a trial

		Returns:
		True or False
		"""

		return len(l) > 2 and l[0] == 'MSG' and l[2] == self.endTrialKey

	def finishTrial(self, trialDict):

		"""
		Perform some finalization after we we have parsed a trial

		Arguments:
		trialDict -- a trial dictionary
		"""

		pass

	def initTrial(self, trialDict):

		"""
		Perform some initalization before we start parsing a trial

		Arguments:
		trialDict -- a trial dictionary
		"""

		pass

	def parseFile(self, path):

		"""
		Parses a single Eyelink ASCII file

		Returns:
		An array
		"""

		# First read all trials into a list of trialDicts
		fd = open(path, 'r')
		lTrialDict = []
		while True:
			s = fd.readline()
			if s == '':
				break
			l = self.strToList(s)
			trialId = self.startTrial(l)
			if trialId != None:
				trialDict = self.parseTrial(trialId, fd)
				if trialDict != None:
					lTrialDict.append(trialDict)

		# Extract the column names
		lVar = []
		for trialDict in lTrialDict:
			for var in trialDict:
				if var not in lVar:
					lVar.append(var)
		lVar.sort()

		# Construct a numpy array with the variable names on the first row
		a = numpy.zeros( (len(lTrialDict)+1, len(lVar)), dtype=self.dtype )
		a[0] = numpy.array(lVar, dtype=self.dtype)
		i = 1
		for trialDict in lTrialDict:
			for var in lVar:
				if var in trialDict:
					a[i, lVar.index(var)] = trialDict[var]
			i += 1

		return a

	def parseLine(self, trialDict, l):

		"""
		Parse a single line (in list format) from a trial

		Arguments:
		trialDict -- the dictionary of trial variables
		l -- a list
		"""

		pass

	def parseTrial(self, trialId, fd):

		"""
		Parse a single trial

		Arguments:
		trialId -- the trial id
		fd -- a file object

		Returns:
		A dictionary like {'var1' : 'value', ...}
		"""

		trialDict = {'trialId' : trialId, 'file' : os.path.basename(fd.name), \
			'outlier' : 0, 'eye_used' : 'undefined'}
		self.initTrial(trialDict)

		inBlink = False
		while True:
			s = fd.readline()
			if s == '':
				break
			l = self.strToList(s)
			if self.endTrial(l):
				self.finishTrial(trialDict)
				return trialDict

			if self.ignoreBlinks:
				if self.startBlink(s):
					inBlink = True
				elif self.endBlink(s):
					inBlink = False

			self.parseVariables(trialDict, l)

			if not inBlink:
				self.parseLine(trialDict, l)

		if self.requireEndTrial:
			raise Exception('Trial %s was started but not ended' \
				% trialDict['trialId'])
		else:
			warnings.warn('Trial %s was started but not ended' \
				% trialDict['trialId'])
			return None

	def parseVariables(self, trialDict, l):

		"""
		Parse a single line (in list format) for variables

		Arguments:
		trialDict -- the dictionary of trial variables
		l -- a list
		"""

		# By default, variables are noted like
		# MSG [time] [variableKey] [name] [value]
		if len(l) > 4 and l[0] == 'MSG' and l[2] == self.variableKey:
			var = l[3]
			val = l[4]			
		# Sometimes, variables are notes like
		# [variableKey] [name] [value]
		elif len(l) > 2 and l[0] == self.variableKey:
			var = l[1]
			val = l[2]
		else:
			return
			
		if var in trialDict:
			warnings.warn('Variable \'%s\' occurs twice in trial %s' \
				% (var, trialDict['trialId']))
		trialDict[var] = val
			
	def strToList(self, s):

		"""
		Converts a string, corresponding to a single line from a file, to a
		list of values.

		Arguments:
		s -- a string

		Returns:
		A list
		"""

		l = []
		for v in s.split():
			try:
				l.append(int(v))
			except:
				try:
					l.append(float(v))
				except:
					l.append(v)
		return l

	def toSaccade(self, l):

		"""
		Attempts to parse a line (in list format) into a dictionary of saccade
		information.

		TODO:
		Handle other fixation formats

		Arguments:
		l -- a list

		Returns:
		None if the list isn't a saccade, otherwise a dictionary with the
		following keys: 'sx', 'sy', 'ex', 'ey', 'sTime', 'eTime', 'duration',
		'size'.
		"""

		if len(l) < 11 or l[0] != "ESACC":
			return None

		try:
			saccade = {}
			if len(l) == 15:
				saccade["sx"] = l[9]
				saccade["sy"] = l[10]
				saccade["ex"] = l[11]
				saccade["ey"] = l[12]
			else:
				saccade["sx"] = l[5]
				saccade["sy"] = l[6]
				saccade["ex"] = l[7]
				saccade["ey"] = l[8]
			saccade["size"] = numpy.sqrt( (saccade['sx']-saccade['ex'])**2 +
				(saccade['sy']-saccade['ey'])**2)
			saccade["sTime"] = l[2]
			saccade["eTime"] = l[3]
			saccade["duration"] = saccade["eTime"] - saccade["sTime"]
			return saccade
		except:
			return None

	def toSample(self, l):

		"""
		Attempts to parse a line (in list format) into a dictionary of sample
		information.

		TODO:
		Handle other fixation formats

		Arguments:
		l -- a list

		Returns:
		None if the list isn't a sample, otherwise a dictionary with the
		following keys: 'x', 'y', 'time'.
		"""

		if len(l) != 5:
			return None
		try:
			sample = {}
			sample["time"] = int(l[0])
			sample["x"] = float(l[1])
			sample["y"] = float(l[2])
			return sample
		except:
			return None

	def toFixation(self, l):

		"""
		Attempts to parse a line (in list format) into a dictionary of fixation
		information.

		TODO:
		Handle other fixation formats

		Arguments:
		l -- a list

		Returns:
		None if the list isn't a fixation, otherwise a dictionary with the
		following keys: 'x', 'y', 'sTime', 'eTime', 'duration'.
		"""

		if len(l) != 8 or l[0] != "EFIX":
			return None
		fixation = {}
		fixation["x"] = l[5]
		fixation["y"] = l[6]
		fixation["sTime"] = l[2]
		fixation["eTime"] = l[3]
		fixation["duration"] = fixation['sTime'] - fixation['eTime']
		return fixation
