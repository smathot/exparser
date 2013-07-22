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
from exparser import TraceKit
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial.distance import euclidean
import warnings

class EyelinkAscFolderReader(BaseReader):

	"""Parses Eyelink ASCII data"""

	def __init__(self, path='data', ext='.asc', startTrialKey='start_trial', \
		endTrialKey='stop_trial', variableKey='var', dtype='|S128', maxN=None, \
		maxTrialId=None, requireEndTrial=True, traceFolder='traces', \
		offlineDriftCorr=False, skipList=[], blinkReconstruct=False):

		"""
		Constructor. Reads all Eyelink ASCII files from a specific folder.

		Keyword arguments:
		path			--	the folder containing the csv files (default='data')
		ext				--	the extension of the csv files (default='.csv')
		startTrialKey	-- 	the start trial keyword (default='start_trial')
		endTrialKey		--	the stop trial keyword (default='stop_trial')
		variableKey		--	the variable keyword (default='var')
		dtype			--	the numpy dtype to be used (default='|S128')
		maxN			--	the maximum number of subjects to process
							(default=None)
		maxTrialId		--	the maximum number of trials to process
							(default=None)
		requireEndTrial	--	indicates whether an exception should be raised if a
							trial hasn't been neatly closed. Otherwise the trial
							is simply disregarded. (default=True)
		traceFolder		--	the folder to save the gaze traces. Traces are saved
							as 3d numpy arrays (x, y, pupil size) in .npy
							format. To start collecting traces, set
							`self.tracePhase` to a value. Use the value
							'__baseline__' to use an automatic baseline.							
							(default='traces')
		offlineDriftCorr	--	Indicates whether coordinates should be
								corrected based on the drift-correction check,
								in case the 'active' drift correctrion is
								disabled (as on the Eyelink 1000).
								(default=False)
		skipList		--	A list of trialIDs that should not be processed.
							(default=[])
		blinkReconstruct	--	Indicates whether pupil size should be
								interpolated during blinks. (default=False)
		"""

		self.startTrialKey = startTrialKey
		self.endTrialKey = endTrialKey
		self.variableKey = variableKey
		self.dtype = dtype
		self.requireEndTrial = requireEndTrial
		self.maxTrialId = maxTrialId
		self.tracePhase = None
		self.traceFolder = traceFolder
		self.traceSmoothParams = None
		self.offlineDriftCorr = offlineDriftCorr
		self.driftAdjust = 0,0
		self.skipList = skipList
		self.blinkReconstruct = blinkReconstruct

		print '\nScanning \'%s\'' % path
		self.dm = None
		nFile = 0
		for fname in os.listdir(path):
			if os.path.splitext(fname)[1] == ext:
				print 'Reading %s ...' % fname,
				a = self.parseFile(os.path.join(path, fname))
				dm = DataMatrix(a)
				if self.dm == None:
					self.dm = dm
				else:
					self.dm += dm
				print '(%d rows)' % len(dm)
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

		return self.dm

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
		Perform some finalization after we we have parsed a trial. This function
		should be overridden by a custom parser.

		Arguments:
		trialDict -- a trial dictionary
		"""

		pass

	def initTrial(self, trialDict):

		"""
		Perform some initalization before we start parsing a trial. This
		function should be overridden by a custom parser.

		Arguments:
		trialDict -- a trial dictionary
		"""

		pass
	
	def __finishTrial__(self, trialDict):

		"""
		Perform some finalization after we we have parsed a trial

		Arguments:
		trialDict -- a trial dictionary
		"""

		nPhase = len(self.traceDict)
		i = 1
		if '--traceplot' in sys.argv or '--traceimg' in sys.argv:
			plt.clf()
			plt.close()
			plt.figure(figsize=(12,12))
			plt.subplots_adjust(hspace=.5, wspace=.5)
		for phase, trace in self.traceDict.iteritems():
			a = np.array(trace, dtype=float)
			origA = a.copy()
			if self.blinkReconstruct:
				a[:,2] = TraceKit.blinkReconstruct(a[:,2])
			if self.traceSmoothParams != None:
				a[:,0] = TraceKit.smooth(a[:,0], **self.traceSmoothParams)
				a[:,1] = TraceKit.smooth(a[:,1], **self.traceSmoothParams)
				a[:,2] = TraceKit.smooth(a[:,2], **self.traceSmoothParams)
			self.traceDict[phase] = a
			path = os.path.join(self.traceFolder, '%s-%.5d-%s.npy' \
				% (trialDict['file'], trialDict['trialId'], phase))
			np.save(path, a)		
			trialDict['__trace_%s__' % phase] = path
			if '--traceplot' in sys.argv or '--traceimg' in sys.argv:
				plt.subplot(nPhase, 3, i)
				i += 1
				plt.title('X(%s)' % phase)
				plt.plot(a[:,0])
				if self.traceSmoothParams != None:
					plt.plot(origA[:,0])
				plt.subplot(nPhase, 3, i)
				i += 1
				plt.title('Y(%s)' % phase)
				plt.plot(a[:,1])
				if self.traceSmoothParams != None:
					plt.plot(origA[:,1])
				plt.subplot(nPhase, 3, i)
				i += 1
				plt.title('Pupil(%s)' % phase)
				plt.plot(a[:,2])	
				if self.traceSmoothParams != None or self.blinkReconstruct:
					plt.plot(origA[:,2])
		if '--traceplot' in sys.argv or '--traceimg' in sys.argv:
			plt.suptitle(path)
			if '--traceimg' in sys.argv:
				path = os.path.join(self.traceFolder, 'png', '%s-%.5d-%s.png' \
					% (trialDict['file'], trialDict['trialId'], phase))
				plt.savefig(path)
			if '--traceplot' in sys.argv:
				plt.show()

	def __initTrial__(self, trialDict):

		"""
		Perform some initalization before we start parsing a trial

		Arguments:
		trialDict -- a trial dictionary
		"""

		self.tracePhase = None
		self.traceDict = {}
		
	def parseDriftCorr(self, l):
		
		"""
		Sets the drift-correction parameters based on the drift correction data.
		This allows you to do real drift correct offline, even when this was not
		enabled during the experiment.
		
		MSG	900353 DRIFTCORRECT R RIGHT at 512,384  OFFSET 0.45 deg.  -15.5,5.3 pix.
		
		Arguments:
		l	--	a list		
		"""
				
		if 'DRIFTCORRECT' in l and len(l) > 10:			
			s = l[10].split(',')
			x = float(s[0])
			y = float(s[1])			
			self.driftAdjust = x, y			

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
			if self.offlineDriftCorr:
				self.parseDriftCorr(l)			
			trialId = self.startTrial(l)
			if self.maxTrialId != None and trialId > self.maxTrialId:
				break
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
		a = np.zeros( (len(lTrialDict)+1, len(lVar)), dtype=self.dtype )
		a[0] = np.array(lVar, dtype=self.dtype)
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
		trialId		--	the trial id
		fd			--	a file object

		Returns:
		A dictionary like {'var1' : 'value', ...}
		"""

		trialDict = {'trialId' : trialId, 'file' : os.path.basename(fd.name), \
			'outlier' : 0, 'eye_used' : 'undefined', 'n_blink' : 0}
		if trialDict['trialId'] not in self.skipList:
			self.__initTrial__(trialDict) # Internal stuff
			self.initTrial(trialDict) # To be overridden
		self.inBlink = False
		while True:
			s = fd.readline()
			if s == '':
				break
			l = self.strToList(s)
			if self.endTrial(l):
				if trialDict['trialId'] not in self.skipList:
					self.__finishTrial__(trialDict) # Internal stuff
					self.finishTrial(trialDict) # To be overridden
				return trialDict
			if trialDict['trialId'] in self.skipList:
				continue
			if self.startBlink(l):
				self.inBlink = True
				trialDict['n_blink'] += 1
			elif self.endBlink(l):
				self.inBlink = False
			self.parseVariables(trialDict, l)
			if self.offlineDriftCorr:
				self.parseDriftCorr(l)
			if self.tracePhase != None:
				self.parseTrace(l)
			self.parseLine(trialDict, l)

		if self.requireEndTrial:
			raise Exception('Trial %s was started but not ended' \
				% trialDict['trialId'])
		else:
			warnings.warn('Trial %s was started but not ended' \
				% trialDict['trialId'])
			return None
			
	def parseTrace(self, l):
		
		"""
		Adds a gaze sample to a trace

		Arguments:
		l -- a list
		"""
		
		s = self.toSample(l)
		if s == None:
			return
		if self.tracePhase not in self.traceDict:
			self.traceDict[self.tracePhase] = []			
		self.traceDict[self.tracePhase].append( (s['x'], s['y'], s['pupil']) )

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
				saccade["sx"] = l[9] - self.driftAdjust[0]
				saccade["sy"] = l[10] - self.driftAdjust[1]
				saccade["ex"] = l[11] - self.driftAdjust[0]
				saccade["ey"] = l[12] - self.driftAdjust[1]
			else:
				saccade["sx"] = l[5] - self.driftAdjust[0]
				saccade["sy"] = l[6] - self.driftAdjust[1]
				saccade["ex"] = l[7] - self.driftAdjust[0]
				saccade["ey"] = l[8] - self.driftAdjust[1]
			saccade["size"] = np.sqrt( (saccade['sx']-saccade['ex'])**2 +
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
		information. The expected format is:
		
		# Timestamp x y pupil size ...
		4815155   168.2   406.5  2141.0 ...

		or (during blinks)
		661781	   .	   .	    0.0	...


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
			if l[1] == '.':
				sample['x'] = np.nan
			else:
				sample["x"] = float(l[1]) - self.driftAdjust[0]
			if l[2] == '.':
				sample["y"] = np.nan
			else:
				sample["y"] = float(l[2]) - self.driftAdjust[1]
			if l[3] == 0:
				sample['pupil'] = np.nan
			else:
				sample['pupil'] = float(l[3])
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
		fixation["x"] = l[5] - self.driftAdjust[0]
		fixation["y"] = l[6] - self.driftAdjust[1]
		fixation["sTime"] = l[2]
		fixation["eTime"] = l[3]
		fixation["duration"] = fixation['sTime'] - fixation['eTime']
		return fixation
