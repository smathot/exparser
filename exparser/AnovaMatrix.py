#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import numpy as np
from exparser.BaseMatrix import BaseMatrix
from exparser import Constants
from rpy2 import robjects
from rpy2.robjects.packages import importr
stats = importr('stats')
base = importr('base')

class AnovaMatrix(BaseMatrix):

	"""A matrix containing the results of an Anova analysis"""

	def __init__(self, dm, factors, dv, subject=None):
	
		"""
		Constructor
		
		Arguments:
		dm -- a DataMatrix
		factors -- a list of factors
		dv -- the dependent variable
		
		Keyword arguments:
		subject -- the variable that should be treated as subject. If this is
				   specified, the Anova will be a SPSS-style repeated measures.
				   If no subject is specified, it will be a regular Anova.
				   (default=None)
		"""
		
		if subject == None:
			raise Exception( \
				'Only within-subject ANOVAs are supported at the moment. Please specify a subject-variable as keyword.')
		
		self.dm = dm
		
		# Construct a formula for R
		if subject != None:
			_factors = factors + [subject]			
			self._formula = '%s ~ %s + Error(%s/(%s))' % (dv, '*'.join(factors),
				subject, '*'.join(factors))
		else:
			_factors = factors
			self._formula = '%s ~ %s' % (dv, '*'.join(factors))
		
		# Collapse the DataMatrix
		inputDm = dm.collapse(_factors, dv)		
		
		# Register all the variables so that R can use them
		for f in _factors:
			robjects.globalenv[f] = robjects.FactorVector(inputDm[f])
		# The dependent variable is renamed to 'mean' by DataMatrix.collapse()
		robjects.globalenv[dv] = robjects.FloatVector(inputDm['mean'])

		# Perform the Anova
		self.aov = stats.aov(robjects.Formula(self._formula))

		# Convert the results into a Matrix for easy reading		
		b = base.summary(self.aov)
		l = [ ['Factor', 'Df', 'F value', 'Pr(>F)', 'sign.'] ]
		for name, vector in b.iteritems():
			name = name[len('Error: '):]
			df = vector[0][0][0]
			f = vector[0][3][0]
			p = vector[0][4][0]
			# R-style significance codes			
			# Signif. codes:  0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1 
			if p < .001:
				sign = '***'
			elif p < .01:
				sign = '**'
			elif p < .05:
				sign = '*'
			elif p < .1:
				sign = '.'
			else:
				sign = ''			
			if str(f) == 'NA':
				f = ''
			if str(p) == 'NA':
				p = ''				
			l.append( [name, df, f, p, sign] )
		self.m = np.array(l, dtype='|S128')
		
	def formula(self):
	
		"""
		Return the formula (model) that was used as input for the R aov()
		function
		
		Returns:
		A formula
		"""
		
		return self._formula

	def RSummary(self):
	
		"""
		Get an R style summary of the ANOVA
		
		Returns:
		A summary string
		"""
	
		return base.summary(self.aov)
		
		
