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

from matplotlib import cm, pyplot as plt
import types
import numpy as np
from exparser.BaseMatrix import BaseMatrix
from exparser import Constants

class PivotMatrix(BaseMatrix):

	"""A PivotMatrix 'summary table' based on a DataMatrix"""

	def __init__(self, dm, cols, rows, dv, func='mean()', err='95ci', \
		colsWithin=False):

		"""
		Constructor. Generates a PivotMatrix that summarizes a SmartArray

		Arguments:
		dm -- a DataMatrix
		cols -- a list of keys for the colums
		rows -- a list of keys for the rows
		dv -- the variable to summarize

		Keyword arguments:
		func -- the function (np.[func]) that is applied on the dependent
				variable (default='mean')
		err -- the type of variance that should be shown: '95ci', 'std', or
			   'se' (default='95ci')
		colsWithin -- specifies whether the column variables should be treated
					  as within subject factors (default=False)
		"""

		self.colHeaders = len(rows)*2
		self.rowHeaders = len(cols)*2
		self.err = err
		self.dv = dv
		self.cols = cols
		self.rows = rows
		self.dm = dm
		self.func = func
		self.colsWithin = colsWithin

		# First create a matrix and fill it with data
		rowGroup = dm.group(rows)						
		self.nRows = len(rowGroup)				
		firstRow = True
		row = 0
		nCols = None
		for dm in rowGroup:
			colGroup = dm.group(cols)			
			self.nCols = len(colGroup)
			if firstRow:
				self.m = np.zeros( (self.rowHeaders+len(rowGroup)+2, \
					self.colHeaders+len(colGroup)+2), dtype=np.float64)
			col = 0			
			for _dm in colGroup:
				if type(func) == types.FunctionType:
					val = func(_dm[dv])
				else:
					val = eval("_dm[dv].%s" % func)	
				self.m[self.rowHeaders+row, self.colHeaders+col] = val					
				col += 1				
			if nCols != None and col != nCols:
				raise Exception('Some cells are empty!')				
			nCols = col				
			row += 1
			firstRow = False

		self.addDescriptives()

		# Optionally, remove within subject variance, assuming that the columns
		# are within subject
		if colsWithin:
			for rowIndex in range(self.nRows):
				for colIndex in range(self.nCols):
					vval = self.m[self.rowHeaders+rowIndex, self.colHeaders+colIndex]
					sMean = self.m[self.rowHeaders+rowIndex, -2]
					gMean = self.m[-2,-2]
					self.m[self.rowHeaders+rowIndex, self.colHeaders+colIndex] \
						= vval - sMean + gMean
			self.addDescriptives()

		# Fill the top-left corner of the matrix
		self.m = np.array(self.m, dtype='|S128')
		for row in range(self.rowHeaders):
			for col in range(self.colHeaders):
				if row != -1 or col != -1:
					self.m[row, col] = '-'
		self.m[-2, -1] = '-'
		self.m[-1, -2] = '-'
		self.m[1:self.rowHeaders, -2:] = '-'
		self.m[-2:, 1:self.colHeaders] = '-'
		self.m[0, -2] = '[MEAN]'
		self.m[0, -1] = '[%s]' % err
		self.m[-2, 0] = '[MEAN]'
		self.m[-1, 0] = '[%s]' % err

		# Add the row headers
		row = self.rowHeaders
		for dm in rowGroup:
			col = 0
			for vname in rows:
				self.m[row, col] = '[%s]' % vname
				self.m[row, col+1] = '[%s]' % dm[vname][0]
				col += 2
			row += 1

		# Add the column headers
		col = self.colHeaders
		for dm in colGroup:
			row = 0
			for vname in cols:
				self.m[row, col] = '[%s]' % vname
				self.m[row+1, col] = '[%s]' % dm[vname][0]
				row += 2
			col += 1

	def addDescriptives(self):

		"""Adds averages and errors to the PivotMatrix"""

		# Determine the row averages and std
		self.rowMeans = []
		self.rowStds = []
		for rowIndex in range(self.nRows):
			row = self.m[self.rowHeaders+rowIndex][self.colHeaders:-2]
			self.rowMeans.append(row.mean())
			self.rowStds.append(row.std())
			self.m[self.rowHeaders+rowIndex][-2] = row.mean()
			self.m[self.rowHeaders+rowIndex][-1] = row.std()

		# Determine the column averages and std
		_m = self.m.swapaxes(0,1)
		self.colMeans = []
		self.colErrs = []		
		for colIndex in range(self.nCols):
			col = _m[self.colHeaders+colIndex][self.rowHeaders:-2]
			_m[self.colHeaders+colIndex][-2] = col.mean()			
			if self.err == '95ci':
				e = col.std()/np.sqrt(col.size)*1.96
			elif self.err == 'se':
				e = col.std()/np.sqrt(col.size)
			elif self.err == 'std':
				e = col.std()
			else:
				raise Exception('Err keyword must be "95ci", "se", or "std"')
			_m[self.colHeaders+colIndex][-1] = e
			self.colMeans.append(col.mean())
			self.colErrs.append(e)

		# Determine the grand average and std
		self.m[-2,-2] = self.m[self.rowHeaders:-2, self.colHeaders:-2].mean()
		self.m[-1,-1] = self.m[self.rowHeaders:-2, self.colHeaders:-2].std()

	def asArray(self):

		"""
		Returns:
		An array representation
		"""

		return self.m
		
	def barPlot(self, fig=None, show=False, _dir='up', barSpacing1=1, \
		barSpacing2=.5, barWidth=.75):
		
		"""
		Draws a bar chart
		
		Keyword arguments:
		fig -- an existing matplotlib figure to draw in (default=None)
		show -- indicates whether the figure should be shown (default=False)
		_dir -- the direction of the bars ('up', 'down', 'left', 'right')
				(default='right')
		barSpacing1 -- the small spacing between bars (default=1)
		barScaping2 -- the extra spacing between groups of bars (default=.5)
		barWidth -- the width of the bars (default=.75)
		
		Returns:
		A matplotlib figure
		"""
			
		# Determine the first (v1) and second (v2) factor
		if len(self.cols) not in (1, 2):
			raise Exception('You can only plot 1 or 2 factor PivotMatrices')
		v1 = self.cols[0]			
		if len(self.cols) == 1:
			v2 = None
		else:
			v2 = self.cols[1]
		
		# Get the mean and error values from the matrix	
		aMean = np.array(self.m[-2,2:-2], dtype=float)
		aErr = np.array(self.m[-1,2:-2], dtype=float)
		
		xLabel = []
		xData = []
		colors = []
		x = 0

		for dm in self.dm.group(v1):		
			colNr = -1
			for _dm in dm.group(v2):			
				colors.append(Constants.palette[colNr])

				# Create a label for the column
				l1 = _dm[v1][0]
				if v2 == None:
					xLabel.append(l1)
				else:
					l2 = _dm[v2][0]
					xLabel.append('%s, %s' % (l1, l2))
						
				# Create an x-coordinate for the column
				x += barSpacing1
				xData.append(x)	
				
				# Advance to the next column (we are counting downwards because
				# this value is used to pick the last colours from the palette
				# list
				colNr -= 1	
				
			# Between conditions there is a gap
			x += barSpacing2	
			
		xData = np.array(xData)	
		
		# Optionally create a new figure
		if fig == None:
			fig = plt.figure()						
			
		if _dir in ('left', 'right'):
			# Draw a horizontal bar chart
			plt.barh(xData, aMean, barWidth, xerr=aErr, color=colors, \
				figure=fig)
			plt.yticks(xData+barWidth/2, xLabel, figure=fig)		
			plt.ylim(xData.min()-barSpacing2, xData.max()+barWidth+barSpacing2)			
			if _dir == 'left':
					plt.gca().invert_xaxis()				
		else:		
			# Draw a vertical (normal) bar chart
			plt.bar(xData, aMean, barWidth, yerr=aErr, color=colors, figure=fig)
			plt.xticks(xData+barWidth/2, xLabel, figure=fig)		
			plt.xlim(xData.min()-barSpacing2, xData.max()+barWidth+barSpacing2)
			if _dir == 'down':
				plt.gca().invert_yaxis()
				
		# Optionally show the figure
		if show:
			plt.show()
		
		return fig
		
	def plot(self, nLvl1=2, size=(5,4), dpi=90, show=False, errBar=True,
		lineLabels=None, xTicks=None, xLabel=None, yLabel=None, grid=False, \
		lineWidth=None, symbols=None, errCap=0, fontFamily='Arial', fontSize=9, \
		colors=None, legendPos='best', yLim=None, plotType='line', barWidth=.2,
		xMargin=.2, barEdgeColor='#000000', title=None, xData=None, fig=None, 
		legendTitle=None, dPlot=None, xTicksRot='horizontal', xTicksFmt='%d'):

		"""
		Plots the current PivotMatrix. The graph is generated such that the
		factors are assumed to be organized into columns. The nLvl1 specifies
		the number of levels of the outermost factor.

		Keyword arguments:
		nLvl1 -- the number of levels of the outer factor (default=2)
		show -- indicates whether the plot should be shown (default=False)
		errBar -- indicates whether error bars should be drawn (default=True)
		size
		dpi
		lineLabels
		xTicks
		xLabel
		yLabel
		grid
		lineWidth
		symbols
		errCap
		fontFamily
		fontSize
		colors
		legendPos
		yLim
		plotType -- 'line' or 'bar' (default=Line)
		barWidth
		barEdgeColor
		xMargin
		title
		xData
		fig
		legendTitle
		dPlot
		xTicksRot
		xTicksFmt

		Returns:
		A matplotlib.plt plot
		"""

		if colors == None:
			colors = Constants.plotLineColors + []
		if symbols == None:
			symbols = Constants.plotLineSymbols + []
		if lineWidth == None:
			lineWidth = Constants.plotLineWidth

		aMean = np.array(self.m[-2,2:-2], dtype=float)
		aMean = aMean.reshape(nLvl1, aMean.size/nLvl1)
		aStd = np.array(self.m[-1,2:-2], dtype=float)
		aStd = aStd.reshape(nLvl1, aStd.size/nLvl1)

		plt.rc("font", family=Constants.fontFamily)
		plt.rc("font", size=Constants.fontSize)

		if fig == None:
			fig = plt.figure(figsize=size, dpi=dpi)
		
		if title != None:
			plt.title(title, figure=fig)

		if grid:
			plt.grid(alpha=0.5)

		if dPlot != None:			
			if len(aMean) != 2:
				raise Exception('In order to plot a difference plot, there must be exactly two lines')
			if dPlot > 1:
				aMean = [aMean[1] - aMean[0]]
				aStd = [.5 * (aStd[1] + aStd[0])]
			else:
				aMean = [aMean[0] - aMean[1]]
				aStd = [.5 * (aStd[0] + aStd[1])]		

		i = 0
		for yData, yErr in zip(aMean, aStd):

			# Choose a label
			if lineLabels != None:
				label = lineLabels[i]
			else:
				label = str(i)

			color = colors.pop()
			symbol = symbols.pop()

			if xData == None:
				if xTicks == None:
					xData = np.linspace(0, 1,  len(yData))
					plt.xlim(-xMargin, 1.+xMargin)
				else:
					if type(xTicks[0]) not in (unicode, str):
						xData = xTicks
						m = (max(xTicks) - min(xTicks)) * .05
						plt.xlim(min(xTicks)-m, max(xTicks)+m)
					else:
						xData = np.linspace(0, len(xTicks)-1, len(yData))
						plt.xlim(-xMargin, len(xTicks)-(1.-xMargin))

			if plotType == 'line':

				# Draw a line graph
				plt.plot(xData, yData, symbol, figure=fig, label=label, \
					color=color, linewidth=lineWidth, markerfacecolor='white', \
					markeredgecolor=color, markeredgewidth=lineWidth)

			elif plotType == 'bar':

				# Draw a bar graph
				xData += i*barWidth - .5*barWidth*len(aMean)
				plt.bar(xData, yData, barWidth, color=color, \
					linewidth=lineWidth, edgecolor=barEdgeColor, label=label)
				xData += .5*barWidth


			if errBar:
				for x,y,err in zip(xData,yData,yErr):
					if plotType == 'bar':
						color = barEdgeColor
					plt.errorbar(x, y, err, fmt=None, figure=fig, ecolor=color, \
						linewidth=lineWidth, capsize=errCap)
			i += 1

		if xTicks != None:
			if type(xTicks[0]) not in (unicode, str):				
				plt.xticks(xTicks, rotation=xTicksRot)
				xTickLabels = [xTicksFmt % x for x in xTicks]
				xTickVals = xTicks
			else:
				xTickVals = np.arange(len(xTicks))
				xTickLabels = xTicks				
			plt.xticks(xTickVals, xTickLabels, rotation=xTicksRot)
			
		if xLabel != None:
			plt.xlabel(xLabel)
		if yLabel != None:
			plt.ylabel(yLabel)
		if yLim != None:
			plt.ylim(yLim)

		if legendPos != None:
			plt.legend(loc=legendPos, title=legendTitle, frameon=False)
		if show:
			plt.show()
		return fig

