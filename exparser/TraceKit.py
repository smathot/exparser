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
from scipy.stats import nanmean, nanmedian, nanstd, ttest_ind, linregress
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import mpl
import warnings
import numpy as np
from exparser.TangoPalette import *
from exparser.RBridge import RBridge
from exparser.Cache import cachedArray, cachedDataMatrix
from exparser.DataMatrix import DataMatrix

def getTrace(dm, signal=None, phase=None, traceLen=None, offset=0,
	lock='start', traceTemplate='__trace_%s__', baseline=None, baselineLen=100,
	baselineOffset=0, baselineLock='end', smoothParams=None, nanPad=True,
	transform=None, deriv=0, regress=None, **dummy):

	"""
	desc:
		Gets a trace for a single trial.

	arguments:
		dm:
			desc:	A DataMatrix with only a single trial (if more trials are in
					there, only the first will be used).
			type:	DataMatrix

	keywords:
		signal:
			desc:	'x', 'y', or 'pupil'.
			type:	[str, unicode, NoneType]
		phase:
			desc:	The name of the phase.
			type:	[str, unicode]
		traceLen:
			desc:	The length of the trace to plot.
			type:	int
		offset:
			desc:	The first (if lock == start) or last (if lock == end)
					samples to skip.
			type:	int
		lock:
			desc:	Indicates whether the trace should be locked from the
					phase start or end.
			type:	[str, unicode]
		traceTemplate:
			desc:	Used to map a phase name onto a DataMatrix key.
			type:	[str, unicode]
		baseline:
			desc:	The phase to use for the baseline.
			type:	[str, unicode, NoneType]
		baselineLen:
			desc:	The length of the baseline.
			type:	int
		baselineOffset:
			desc:	The first (if lock == start) or last (if lock == end)
					baseline samples to skip.
			type:	int
		baselineLock:
			desc:	Indicates whether the baseline should be locked from the
					phase start or end.
			type:	[str, unicode]
		smoothParams:
			desc:	A {'windowLen' : [..], 'windowType' : [..]} dictionary
					that is used to specify signal smoothing (see smooth()),
					or `None` for no smoothing.
			type:	[dict, NoneType]
		nanPad:
			desc:	If set to True, traces that are shorter than traceLen
					are padded with np.nan values.
			type:	bool
		transform:
			desc:	A function to transform the trace value.
			type:	[NoneType, function]
		deriv:
			desc:	Indicates the derivative that should be used. The 1st
					derivative corresponds to velocity, the 2nd to acceleration.
			type:	int
		regress:
			desc:	A function that regresses the X and Y coordinate out of
					pupil size.
			type:	FunctionType

	returns:
		desc:	A 1D NumPy array with the trace.
		type:	ndarray
	"""

	if len(dm) != 1:
		raise Exception('DataMatrix must have exactly one row')
	if signal == None or phase == None or traceLen == None:
		raise Exception('signal, phase, and traceLen are required keywords')
	if lock not in ['start', 'end']:
		raise Exception('lock should be start or end')
	if baselineLock not in ['start', 'end']:
		raise Exception('baselineLock should be start or end')
	if signal == 'x':
		i = 0
	elif signal == 'y':
		i = 1
	elif signal == 'pupil':
		i = 2
	else:
		raise Exception('Invalid signal!')
	# Get the trace
	npy = dm[traceTemplate % phase][0]
	if not os.path.exists(npy):
		raise Exception('Missing .npy trace file: %s (path="%s")' \
			% (traceTemplate % phase, npy))
	aFull = np.load(npy)
	if regress != None:
		_aTrace = regress(aFull, dm)
	else:
		_aTrace = aFull[:,i]
	if smoothParams != None:
		try:
			_aTrace = smooth(_aTrace, **smoothParams)
		except:
			warnings.warn('Failed to smooth trace of length %d' % len(_aTrace))
	if deriv > 0:
		_aTrace = traceDeriv(_aTrace, deriv)
	if lock == 'start':
		_aTrace = _aTrace[offset:offset+traceLen]
	elif offset > 0:
		_aTrace = _aTrace[-offset-traceLen:-offset]
	else:
		_aTrace = _aTrace[-traceLen:]
	# Paste the trace into a nan-filled trace that has exactly the desired
	# length. This is necessary to deal with traces that are shorter than the
	# specified traceLen.
	if nanPad:
		aTrace = np.empty(traceLen)
		aTrace[:] = np.nan
		if lock == 'start':
			aTrace[:len(_aTrace)] = _aTrace
		else:
			aTrace[-len(_aTrace):] = _aTrace
	else:
		aTrace = _aTrace
	# Optionally apply a transform
	if transform != None:
		aTrace = transform(aTrace)
	# If we don't apply a baseline then return right away, possible after
	# smoothing
	if baseline == None:
		if smoothParams != None:
			aTrace = smooth(aTrace, **smoothParams)
		return aTrace
	# Get the baseline
	npy = dm[traceTemplate % baseline][0]
	if not os.path.exists(npy):
		raise Exception('Missing .npy trace file: %s (path="%s")' \
			% (traceTemplate % baseline, npy))
	aFull = np.load(npy)
	if regress != None:
		aBaseline = regress(aFull, dm)
	else:
		aBaseline = aFull[:,i]
	if smoothParams != None:
		aBaseline = smooth(aBaseline, **smoothParams)
	if deriv > 0:
		aBaseline = traceDeriv(aBaseline, deriv)
	if baselineLock == 'start':
		aBaseline = aBaseline[baselineOffset:baselineOffset+baselineLen]
	elif baselineOffset == 0:
		aBaseline = aBaseline[-baselineLen:]
	else:
		aBaseline = aBaseline[-baselineOffset-baselineLen:-baselineOffset]
	if transform != None:
		aBaseline = transform(aBaseline)
	mBaseline = aBaseline.mean()
	aTrace /= mBaseline
	return aTrace

def getTraceAvg(dm, avgFunc=nanmean, **traceParams):

	"""
	Gets a single average trace

	Arguments:
	dm				--	a DataMatrix

	Keyword arguments:
	avgFunc			--	the function to use to determine the average trace. This
						function must be robust to nan values. (default=nanmean)
	*traceParams	--	see getTrace()

	Returns:
	An (xData, yData, errData) tuple, where errData contains the standard
	error.
	"""

	traceLen = traceParams['traceLen']
	mTrace = np.empty( (len(dm), traceLen) )
	mTrace[:] = np.nan
	i = 0
	for trialDm in dm:
		aTrace = getTrace(trialDm, **traceParams)
		mTrace[i, 0:len(aTrace)] = aTrace
		i += 1
	xData = np.linspace(0, traceLen, traceLen)
	yData = nanmean(mTrace, axis=0)
	errData = nanstd(mTrace, axis=0) / np.sqrt(mTrace.shape[0])
	errData = np.array( [errData, errData] )
	return xData, yData, errData


def getTracePeak(dm, peakFunc=np.nanmax, **traceParams):

	a = getTrace(dm, **traceParams)
	i = np.where(a == peakFunc(a))
	xPeak = i[0][0]
	yPeak = a[xPeak]
	return xPeak, yPeak

def getTracePeakAvg(dm, **traceParams):

	"""
	desc:
		Gets the average peak value from.

	arguments:
		dm:
			desc:	A DataMatrix.
			type:	DataMatrix

	keyword-dict:
		*traceParams:	See getTrace()

	Returns:
		An (xPeak, yPeak, xErr, yErr) tuple, where xErr and yErr contain the
		standard error.
	"""

	traceLen = traceParams['traceLen']
	mTrace = np.empty( (len(dm), traceLen) )
	mTrace[:] = np.nan
	i = 0
	aXPeak = np.empty(len(dm))
	aYPeak = np.empty(len(dm))
	for i, trialDm in enumerate(dm):
		xPeak, yPeak = getTracePeak(trialDm, **traceParams)
		aXPeak[i] = xPeak
		aYPeak[i] = yPeak
	return aXPeak.mean(), aYPeak.mean(), np.std(aXPeak)/np.sqrt(len(aXPeak)), \
		np.std(aYPeak)/np.sqrt(len(aYPeak))

def plotTraceAvg(ax, dm, avgFunc=nanmean, lineColor=blue[0], lineStyle='-',
	errColor=gray[1], errAlpha=.4, label=None, _downSample=None, aErr=None,
	orientation='horizontal', **traceParams):

	"""
	Plots a single average trace

	Arguments:
	ax				--	a Matplotlib axis
	dm				--	a DataMatrix

	Keyword arguments:
	avgFunc			--	see getTraceAvg()
	lineColor		--	the line color (default=blue[0])
	lineStyle		--	the line style (default='-')
	errColor		--	the color for the error shading (default=gray[1])
	errAlpha		--	the opacity for the error shading (default=.4)
	traceTemplate	--	is used to determine the correct key from the DataMatrix
						based on the phase (default='__trace_%s__')
	label			--	a line label (default=None)
	_downSample		--	specify a decrease in resolution, to decrease the size
						of the plot. (default=None)
	aErr			--	a 2-D array to use to draw the error shading, or None
						to use the error data calculated by getTraceAvg()
	orientation		--	'horizontal' or 'vertical'. (default='horizontal')
	*traceParams	--	see getTrace()
	"""

	xData, yData, errData = getTraceAvg(dm, avgFunc=avgFunc, **traceParams)
	if aErr != None:
		errData = aErr
	if _downSample != None:
		xData = downSample(xData, _downSample)
		yData = downSample(yData, _downSample)
		errData = downSample(errData, _downSample)
	if orientation == 'horizontal':
		ax.plot(xData, yData, color=lineColor, label=label, linestyle=lineStyle)
		if errColor != None:
			ax.fill_between(xData, yData-errData[0], yData+errData[1], \
				color=errColor, alpha=errAlpha)
	else:
		ax.plot(yData, xData, color=lineColor, label=label)
		if errColor != None:
			ax.fill_betweenx(xData, yData-errData[0], yData+errData[1],
				color=errColor, alpha=errAlpha)

def plotTraceContrast(dm, select1, select2, color1=blue[1], color2=orange[1],
	colorDiff=green[1], label1=None, label2=None, labelDiff=None,
	errAlpha=.25, model=None, showAbs=True, showDiff=False, **params):

	"""
	desc: |
		Creates a trace-contrast plot, with two lines and error shadings.

		NOTE: Passing the `cacheId` keyword will cause the lmer statistics to be
		cached.

	arguments:
		dm:
			desc:	A DataMatrix.
			type:	DataMatrix
		select1:
			desc:	A select statement for the first trace.
			type:	[str, unicode]
		select2:
			desc:	A select statement for the second trace.
			type:	[str, unicode]

	keywords:
		color1:
			desc:	A color for the first trace.
			type:	[str, unicode]
		color2:
			desc:	A color for the second trace.
			type:	[str, unicode]
		colorDiff:
			desc:	A color for the difference trace.
			type:	[str, unicode]
		label1:
			desc:	A label for the first trace.
			type:	[str, unicode, NoneType]
		label2:
			desc:	A label for the second trace.
			type:	[str, unicode, NoneType]
		labelDiff:
			desc:	A label for the difference trace.
			type:	[str, unicode, NoneType]
		errAlpha:
			desc:	Opacity level for the error bars.
			type:	[float, int]
		model:
			desc:	A statistical model to be passed onto `mixedModelTrace()` or
					`None` to skip statistics.
			type:	[str, unicode, NoneType]
		showAbs:
			desc:	Indicates whether the absolute traces (i.e. traces 1 and 2)
					should be shown.
			type:	bool
		showDiff:
			desc:	Indicates whether the difference trace should be shown.
			type:	bool
	"""

	x1, y1, err1 = getTraceAvg(dm.select(select1, verbose=False), **params)
	x2, y2, err2 = getTraceAvg(dm.select(select2, verbose=False), **params)
	y3 = y2-y1
	if model != None:
		ldm = mixedModelTrace(dm, model=model, **params)
		ldm = ldm.select('effect == "%s"' % ldm['effect'][1])
		aSe = ldm['se']
		aT = ldm['t']
		y1min = y1 - aSe/2
		y1max = y1 + aSe/2
		y2min = y2 - aSe/2
		y2max = y2 + aSe/2
		y3min = y3 - aSe
		y3max = y3 + aSe
		if showAbs:
			plt.fill_between(x1, y1min, y1max, color=color1, alpha=errAlpha)
			plt.fill_between(x2, y2min, y2max, color=color2, alpha=errAlpha)
		if showDiff:
			plt.fill_between(x1, y3min, y3max, color=colorDiff, alpha=errAlpha)
		markStats(plt.gca(), aT, **params)
	if showAbs:
		plt.plot(x1, y1, color=color1, label=label1)
		plt.plot(x2, y2, color=color2, label=label2)
	if showDiff:
		plt.plot(x1, y3, color=colorDiff, label=labelDiff)

def traceDiff(dm, select1, select2, epoch=None, **traceParams):

	"""
	Deteremines the difference in the trace between two subsets of the data
	(i.e. two groups or conditions) within a particular epoch.

	Arguments:
	dm				--	A DataMatrix.
	select1			--	A select statement for the first trace.
	select2			--	A select statement for the second trace.

	Keyword arguments:
	epoch			--	The time interval for which to estimate the trace
						difference. This can be None to use the entire trace,
						a single int to take one sample, or a (start, end) tuple
						to select a particular epoch. (default=None)
	traceParams		--	The trace parameters. (default=trialParams)

	Returns:
	A single value reflecting the trace difference.
	"""

	x1, y1, err1 = getTraceAvg(dm.select(select1, verbose=False), \
		**traceParams)
	x2, y2, err2 = getTraceAvg(dm.select(select2, verbose=False), \
		**traceParams)
	d = y1-y2
	if type(epoch) == int:
		return d[epoch]
	if type(epoch) == tuple and len(epoch) == 2:
		d = d[epoch[0]:epoch[1]]
	elif epoch != None:
		raise Exception('Epoch should be None, int, or (int, int)')
	return d.mean()

@cachedDataMatrix
def mixedModelTrace(dm, model, winSize=1, effectIndex=1, **traceParams):

	"""
	desc:
		Perform a mixed model over a single trace. The dependent variable is
		specifed through the signal and phase keywords.

	arguments:
		dm:
			desc:	A DataMatrix.
			type:	DataMatrix
		model:
			desc:	An lmer-style model. This needs to be only the fixed and
					random effects part of the model, so everything after
					the `~` sign. For example `cond + (1|subject_nr)`.
			type:	str

	keywords:
		winSize:
			desc:	Indicates the number of samples that should be skipped
					each time. For a real analysis, this should be 1, but
					for a quick look, it can be increased (default=1)
			type:	int

	keyword-dict:
		*traceParams:	See getTrace().

	returns: |
		A DataMatrix with nr-of-effects*nr-of-samples rows and the following
		columns:

		- `i`: sample number
		- `effect`: name of effect
		- `est`: estimated effect (slope/ intercept)
		- `se`: standard error of effect
		- `t`: t-value of effect
	"""

	if not model.startswith('mmdv__ ~ '):
		model = 'mmdv__ ~ ' + model
	global R
	try:
		R
	except:
		R = RBridge()
	traceLen = traceParams['traceLen']
	l = [ ['i', 'effect', 'est', 'se', 't'] ]
	for i in range(0, traceLen, winSize):
		# First calculate the mean value for the current signal slice for each
		# trial and save that in a copy of the DataMatrix
		_dm = dm.addField('mmdv__', dtype=float)
		for trialId in range(len(_dm)):
			aTrace = getTrace(_dm[trialId], **traceParams)
			if i < len(aTrace):
				sliceMean = aTrace[i:i+winSize].mean()
			else:
				sliceMean = np.nan
			_dm['mmdv__'][trialId] = sliceMean
		# Do mixed effects
		R.load(_dm)
		_dm = R.lmer(model)
		_dm._print(sign=4, title='%d - %d' % (i, i+winSize))
		for k in range(winSize):
			if i+k >= traceLen:
				break
			for j in _dm.range():
					l.append([i+k, _dm['effect'][j], _dm['est'][j], _dm['se'][j],
						_dm['t'][j]])
	return DataMatrix(l)

def statsTrace(dm, key='t', intercept=False, show=True):

	colors = brightColors[:]
	for effect in dm.unique('effect'):
		if not intercept and effect == '(Intercept)':
			continue
		_dm = dm.select('effect == "%s"' % effect)
		plt.plot(_dm['t'], label=effect, color=colors.pop())
	plt.legend(frameon=False)
	if show:
		plt.show()

def markStats(ax, aStat, below=False, _abs=True, thr=2., minSmp=200,
	color=gray[1], alpha=.2, loExt=False, hiExt=False, showSpurious=False,
	**args):

	"""
	Marks all timepoints in a figure with colored shading when the significance
	falls below an alpha threshold.

	Arguments:
	ax		--	a matplotlib axis
	aStat	--	an array with statistics, such as p-values or t-values.

	Keyword arguments:
	below		--	Indicates whether values blow the threshold (True) or above
					it (False) are considered to be signficant. (default=True)
	thr 		--	the threshold (default=.01)
	minSmp		--	the minimum number of consecutive significant samples
					(default=10)
	color		--	the color for the shading (default=gray[1])
	alpha		--	The alpha value. (default=.2)
	loExt		--	Indicates whether intervals before the signal should be
					treated as significant, such that at the start of the signal
					there is no minimum number of samples. (default=False)
	hiExt		--	Like `loExt`, but for the end of the signal. (default=False)
	showSpurious	-- 	Indicates whether 'spurious' signficant points should be
						shown as well. I.e. those that are do not meet the
						minSmp threshold.

	Returns:
	A list of (start, end) tuples with significant regions
	"""

	lRoi = []
	iFrom = None
	aSign = np.zeros(len(aStat))
	aSpurious = np.zeros(len(aStat))
	if _abs:
		aStat = np.abs(aStat)
	for i in range(len(aStat)):
		pVal = aStat[i]
		hit = (pVal < thr and below) or (pVal > thr and not below)
		if hit:
			aSpurious[i] = 1
			if iFrom == None:
				iFrom = i
		if ((not hit or (i == len(aStat)-1)) and iFrom != None):
			if i-iFrom >= minSmp-1 or (iFrom == 0 and loExt) or \
				(i == len(aStat)-1 and hiExt):
				print 'Significant region: %d - %d' % (iFrom, i-1)
				ax.axvspan(iFrom, i-1, ymax=1, color=color, zorder=-9999,
					alpha=alpha)
				lRoi.append((iFrom, i-1))
				aSign[iFrom:i-1] = 1
			iFrom = None
	if showSpurious:
		aSpurious[np.where(aSign == 1)] = 0
		for i in np.where(aSpurious == 1)[0]:
			print 'Spurious %d' % i
			ax.axvline(i, color='black')
	return lRoi

def smooth(aTrace, windowLen=11, windowType='hanning', correctLen=True):

	"""
	Source: <http://www.scipy.org/Cookbook/SignalSmooth>

	Smooth the data using a window with requested size.

	This method is based on the convolution of a scaled window with the signal.
	The signal is prepared by introducing reflected copies of the signal
	(with the window size) in both ends so that transient parts are minimized
	in the begining and end part of the output signal.

	Arguments:
	aTrace		--	an array with the input signal

	Keyword arguments:
	windowLen	--	the dimension of the smoothing window; should be an odd
					integer (default=5)
	windowType	--	the type of window from 'flat', 'hanning', 'hamming',
					'bartlett', 'blackman'. Flat window will produce a moving
					average smoothing.
	correctLen	--	indicates whether the return string should be the same
					length as the input string (default=True).

	Returns:
	An array with the smoothed signal
	"""

	if aTrace.ndim != 1:
		raise ValueError("smooth only accepts 1 dimension arrays.")
	if aTrace.size < windowLen:
		raise ValueError("Input vector needs to be bigger than window size.")
	if windowLen < 3:
		return aTrace
	if not windowType in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError( \
			"Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
	s = np.r_[aTrace[windowLen-1:0:-1], aTrace, aTrace[-1:-windowLen:-1]]
	if windowType == 'flat': #moving average
		w = np.ones(windowLen, 'd')
	else:
		func = getattr(np, windowType)
		w = func(windowLen)
	y = np.convolve(w/w.sum(), s, mode='valid')
	if correctLen:
		y = y[(windowLen/2-1):-(windowLen/2)]
		# The output array can be one shorter than the input array
		if len(y) > len(aTrace):
			y = y[:len(aTrace)]
		elif len(y) < len(aTrace):
			raise Exception('The output array is too short!')
	return y

def downSample(aTrace, i):

	"""
	Downsamples an array by skipping samples.

	Arguments:
	aTrace		--	input array
	i			--	downsampling ratio

	Returns:
	A downsampled array
	"""

	if len(aTrace.shape) == 1:
		return aTrace[::i]
	elif len(aTrace.shape) == 2:
		l = []
		for d in range(aTrace.shape[0]):
			l.append(aTrace[d,::i])
		return np.array(l)
	else:
		raise Exception('Only 1 and 2-dimensional arrays are allowed')

def traceDeriv(aTrace, deriv=0):

	"""
	desc:
		Gets the N-th derivative of a signal. Here, the 1st derivative is
		velocity, the 2nd acceleration, etc.

	arguments:
		aTrace:
			desc:	Input array.
			type:	ndarray

	keywords:
		deriv:
			desc:	The derivative.
			type:	int

	returns:
		desc:	An array with the derivative.
		type:	ndarray
	"""

	aTrace.dtype = np.float
	while deriv > 0:
		aTrace = aTrace[1:]-aTrace[:-1]
		deriv -= 1
	return aTrace

def latency(aTrace, at=None, vt=None, plot=False):

	"""
	Determines the response latency in a signal, based on an accelation and/ or
	velocity threshold.

	Arguments:
	aTrace		--	input array

	Keyword arguments:
	at			--	acceleration threshold (default=None)
	vt	 		--	velocity threshold (default=None)
	plot		--	indicates whether a plot should be shown (default=False)

	Returns:
	The first sample where the acceleration or velocity threshold is exceeded
	"""

	if at == None and vt == None:
		raise Exception( \
			'You must specify an accelation and/or velocity threshold')

	velTrace = aTrace[1:] - aTrace[:-1]
	accTrace = velTrace[1:] - velTrace[:-1]
	aLat = None
	vLat = None

	if vt != None:
		l = np.where(np.abs(velTrace) > vt)[0]
		if len(l) > 0:
			vLat = l[0]

	if at != None:
		l = np.where(np.abs(accTrace) > at)[0]
		if len(l) > 0:
			aLat = l[0]

	if aLat == None and vLat == None:
		lat = None
	elif aLat == None:
		lat = vLat
	elif vLat == None:
		lat = aLat
	else:
		lat = min(aLat, vLat)

	if plot:
		plt.subplot(311)
		plt.plot(aTrace)
		if lat != None:
			plt.axvline(lat)
		plt.subplot(312)
		plt.plot(velTrace)
		if vt != None:
			plt.axhline(vt, color='red')
		if lat != None:
			plt.axvline(lat)
		plt.axhline()
		plt.subplot(313)
		plt.plot(accTrace)
		if at != None:
			plt.axhline(at, color='red')
		if lat != None:
			plt.axvline(lat)
		plt.axhline()
		plt.show()

	return lat

def blinkReconstruct(aTrace, vt=5, maxDur=500, margin=10, plot=False):

	"""
	Reconstructs pupil size during blinks.

	Arguments:
	aTrace		--	The input trace.

	Keyword arguments:
	pvt			--	The pupil velocity threshold. Lower tresholds more easily
					trigger blinks (default=5.)
	maxDur		--	The maximum duration for a blink. Longer blinks are
					ignored. (default=500)
	plot		--	Indicates whether the algorithm should be plotted.
					(default=False)

	Returns:
	An array with the reconstructed pupil data or an (array, figure) tuple when
	plot==True.
	"""

	# Create a copy of the signal, a smoothed version, and calculate the
	# velocity profile.
	aTrace = np.copy(aTrace)
	try:
		sTrace = smooth(aTrace, windowLen=21)
	except Exception as e:
		warnings.warn(str(e))
		sTrace = aTrace

	vTrace = sTrace[1:]-sTrace[:-1]

	if plot:
		plt.clf()
		fig = plt.figure(figsize=(10,5))
		plt.rc("font", family='Liberation Sans')
		plt.rc("font", size=10)
		plt.subplots_adjust(wspace=.25, hspace=.4)
		plt.subplot(2,2,1)
		plt.title('Original signal')
		plt.plot(aTrace, color=blue[1])
		plt.xlabel('Time (ms)')
		plt.ylabel('Pupil size (arbitrary units)')
		plt.subplot(2,2,2)
		plt.title('Smoothed signal')
		plt.plot(sTrace, color=blue[1])
		plt.xlabel('Time (ms)')
		plt.ylabel('Pupil size (arbitrary units)')
		plt.subplot(2,2,3)
		plt.title('Velocity profile')
		plt.plot(vTrace, color=blue[1])
		plt.xlabel('Time (ms)')
		plt.ylabel('Velocity (arbitrary units)')

	# Start blink detection
	iFrom = 0
	lBlink = []
	while True:
		# The onset of the blink is the moment at which the pupil velocity
		# exceeds the threshold.
		l = np.where(vTrace[iFrom:] < -vt)[0]
		if len(l) == 0:
			break # No blink detected
		iStart = l[0]+iFrom
		if iFrom == iStart:
			break
		# The reversal period is the moment at which the pupil starts to dilate
		# again with a velocity above threshold.
		l = np.where(vTrace[iStart:] > vt)[0]
		if len(l) == 0:
			iFrom = iStart
			continue
		iMid = l[0]+iStart
		# The end blink period is the moment at which the pupil velocity drops
		# back to zero again.
		l = np.where(vTrace[iMid:] < 0)[0]
		if len(l) == 0:
			iFrom = iMid
			continue
		iEnd = l[0]+iMid
		iFrom = iEnd
		# We generally underestimate the blink period, so compensate for this
		if iStart-margin >= 0:
			iStart -= margin
		if iEnd+margin < len(aTrace):
			iEnd += margin
		# We don't accept blinks that are too long, because blinks are not
		# generally very long (although they can be).
		if iEnd-iStart > maxDur:
			continue
		if plot:
			plt.axvspan(iStart, iEnd, color=gray[-1], alpha=.4)
		lBlink.append( (iStart, iEnd) )

	if plot:
		plt.subplot(2,2,4)
		plt.title('Reconstructed signal')

	# Now reconstruct the trace during the blinks
	for iStart, iEnd in lBlink:
		# First create a list of (when possible) four data points that we can
		# use for interpolation.
		dur = iEnd - iStart
		l = []
		if iStart-dur >= 0:
			l += [iStart-dur]
		l += [iStart, iEnd]
		if iEnd+dur < len(sTrace):
			l += [iEnd+dur]
		x = np.array(l)
		# If the list is long enough we use cubic interpolation, otherwise we
		# use linear interpolation
		y = aTrace[x]

		if plot:
			plt.plot(x, y, 'o', color=orange[1])

		if len(x) >= 4:
			f2 = interp1d(x, y, kind='cubic')
		else:
			f2 = interp1d(x, y)
		xInt = np.arange(iStart, iEnd)
		yInt = f2(xInt)
		aTrace[xInt] = yInt

	if plot:
		plt.plot(aTrace, color=blue[1])
		plt.xlabel('Time (ms)')
		plt.ylabel('Pupil size (arbitrary units)')

	if plot:
		return aTrace, fig
	return aTrace

@cachedDataMatrix
def splitTrace(dm, splitCol, phase, phaseBefore=None, phaseAfter=None, \
	traceTemplate='__trace_%s__'):

	"""
	Splits all traces in the DataMatrix by the value in a specific column. This
	allows, for example, to split a single trace at the point at which a
	response was given.

	NOTE: This function is cached.

	Arguments:
	splitCol	--	The column that contains the values to split by.
	phase		--	The name of the phase to split.

	Keyword arguments:
	phaseBefore	--	The name of the phase that will be the first part of the
					split. None means that this part of the split will not be
					kept. (default=None)
	phaseAfter	--	The name of the phase that will be the second part of the
					split. None means that this part of the split will not be
					kept. (default=None)
	traceTemplate	--	See `getTrace()`.
	**traceParams	--	See `getTrace()`.

	Returns:
	The DataMatrix with the splitted phased added.
	"""

	if phaseBefore == None and phaseAfter == None:
		raise Exception('Either phaseBefore or phaseAfter should be specified')
	if phaseBefore != None:
		dm = dm.addField(traceTemplate % phaseBefore, dtype=str)
	if phaseAfter != None:
		dm = dm.addField(traceTemplate % phaseAfter, dtype=str)
	for i in dm.range():
		npy =  dm[traceTemplate % phase][i]
		a = np.load(npy)
		split = dm[splitCol][i]
		aBefore = a[:split]
		aAfter = a[split:]
		if phaseBefore != None:
			npyBefore = npy + '.before.npy'
			np.save(npyBefore, aBefore)
			dm[traceTemplate % phaseBefore][i] = npyBefore
		if phaseAfter != None:
			npyAfter = npy + '.after.npy'
			np.save(npyAfter, aAfter)
			dm[traceTemplate % phaseAfter][i] = npyAfter
	return dm
