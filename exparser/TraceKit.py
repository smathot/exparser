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

from scipy.stats import nanmean, nanmedian, nanstd, ttest_ind, linregress
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import mpl
import numpy as np
from exparser.TangoPalette import *

def getTrace(dm, signal=None, phase=None, traceLen=None, offset=0, \
	lock='start', traceTemplate='__trace_%s__', baseline=None, \
	baselineLen=100, baselineOffset=0, baselineLock='end', smoothParams=None):
		
	"""
	Gets a trace for a single trial
	
	Arguments:
	dm	--	a DataMatrix with only a single trial (if more trials are in there,
			only the first will be used)
	
	Keyword arguments:
	signal			--	'x', 'y', or 'pupil' (default=None)
	phase			--	the name of the phase (default=None)
	traceLen		--	the length of the trace to plot (default=None)
	offset			--	the first (if lock == start) or last (if lock == end)
						samples to skip (default=0)
	lock			--	indicates whether the trace should be locked from the
						phase start or end (default='start')
	traceTemplate	--	is used to determine the correct key from the DataMatrix
						based on the phase (default='__trace_%s__')		
	baseline		--	the phase to use for the baseline (default=None)
	baselineLen		--	the length of the baseline (default=100)
	baselineOffset	--	the first (if lock == start) or last (if lock == end)
						baseline samples to skip (default=0)
	baselineLock	--	indicates whether the baseline should be locked from the
						phase start or end (default='start')
	smoothParams	--	a {'windowLen' : [..], 'windowType' : [..]} dictionary
						that is used to specify signal smoothing (see smooth()),
						or None for no smoothing. (default=None)
	
	Returns:
	a 1D NumPy array with the trace	
	"""
	
	if signal == None or phase == None or traceLen == None:
		raise Exception('signal, phase, and traceLen are required keywords')
		
	if signal == 'x':
		i = 0
	elif signal == 'y':
		i = 1
	elif signal == 'pupil':
		i = 2
	else:
		raise Exception('Invalid signal!')
	npy = dm[traceTemplate % phase][0]
	aTrace = np.load(npy)[:,i]
	if lock == 'start':
		aTrace = aTrace[offset:offset+traceLen]
	elif offset > 0:
		aTrace = aTrace[-offset-traceLen:-offset]
	else:
		aTrace = aTrace[-traceLen:]
	if baseline == None:
		if smoothParams != None:
			aTrace = smooth(aTrace, **smoothParams)
		return aTrace
	npy = dm[traceTemplate % baseline][0]
	aBaseline = np.load(npy)[:,i]
	if baselineLock == 'start':
		aBaseline = aBaseline[baselineOffset:baselineOffset+baselineLen]
	elif baselineOffset == 0:
		aBaseline = aBaseline[-baselineLen:]
	else:
		aBaseline = aBaseline[-baselineOffset-baselineLen:-baselineOffset]
	mBaseline = aBaseline.mean()
	aTrace /= mBaseline
	if smoothParams != None:
		aTrace = smooth(aTrace, **smoothParams)
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

def plotTraceAvg(ax, dm, avgFunc=nanmean, lineColor=blue[0], errColor=gray[1], \
	errAlpha=.4, label=None, _downSample=None, aErr=None, \
	orientation='horizontal', **traceParams):
		
	"""
	Plots a single average trace
	
	Arguments:
	ax				--	a Matplotlib axis
	dm				--	a DataMatrix
	
	Keyword arguments:
	avgFunc			--	see getTraceAvg()
	lineColor		--	the line color (default=blue[0])
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
		ax.plot(xData, yData, color=lineColor, label=label)
		if errColor != None:
			ax.fill_between(xData, yData-errData[0], yData+errData[1], \
				color=errColor, alpha=errAlpha)
	else:
		ax.plot(yData, xData, color=lineColor, label=label)
		if errColor != None:
			ax.fill_betweenx(xData, yData-errData[0], yData+errData[1], \
				color=errColor, alpha=errAlpha)

def mixedModelTrace(dm, fixedEffects, randomEffects, winSize=1, nSim=1000, \
	effectIndex=0, **traceParams):

	"""
	Perform a mixed model over a single trace. The dependent variable is
	specifed through the signal and phase keywords.

	Arguments:
	dm				--	a DataMatrix
	fixedEffects	-- 	a list of fixed effects (i.e. the indendendent
						variables)
	randomEffects	--	a list of random effects, such as subject or item
	
	Keyword arguments:
	winSize			--	indicates the number of samples that should be skipped
						each time. For a real analysis, this should be 1, but
						for a quick look, it can be increased (default=1)
	nSim			--	the number of similuations. This should be increased
						for more accurate estimations (default=100)
	effectIndex		--	The index of the effect that should be read from the
						MixedEffectsMatrix. For example, if there are two
						fixed effects, than 0 and 1 are the main effects, and
						2 is the interaction. (default=0)
	*traceParams	--	see getTrace()
	
	Returns:
	A (traceLen, 3) array, where the columns are
	[p-value, 95low, 95high]
	"""

	from exparser.MixedEffectsMatrix import MixedEffectsMatrix
	
	traceLen = traceParams['traceLen']		
	aPVal = np.zeros( (traceLen, 3) )
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
		mem = MixedEffectsMatrix(_dm, 'mmdv__', fixedEffects, \
			randomEffects, nSim=nSim)
		print mem
		pVal = float(mem.asArray()[2+effectIndex][6])
		ciLow = float(mem.asArray()[2+effectIndex][3])
		ciHigh = float(mem.asArray()[2+effectIndex][4])
		aPVal[i:i+winSize,0] = pVal
		aPVal[i:i+winSize,1] = ciLow
		aPVal[i:i+winSize,2] = ciHigh
		print '%.4d: p = %.3f (%f - %f)' % (i, pVal, ciHigh, ciLow)
		
	return aPVal

def markStats(ax, aPVal, alpha=.05, minSmp=200, color=gray[1]):
	
	"""
	Marks all timepoints in a figure with colored shading when the significance
	falls below an alpha threshold.
	
	Arguments:
	ax		--	a matplitlin axis	
	aPVal	--	an array with p-values
	
	Keyword arguments:
	alpha	--	the alpha threshold (default=.01)
	minSmp	--	the minimum number of consecutive significant samples
				(default=10)
	color	--	the color for the shading (default=gray[1])
	"""
	
	iFrom = None
	for i in range(len(aPVal)):
		pVal = aPVal[i]		
		if pVal < alpha:			
			if iFrom == None:
				iFrom = i				
		if (pVal > alpha or i == len(aPVal)-1) and iFrom != None:
			if i-iFrom >= minSmp-1:
				ax.axvspan(iFrom, i-1, ymax=1, color=color, zorder=-9999)
			iFrom = None
						
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
	sTrace = smooth(aTrace, windowLen=21)
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
