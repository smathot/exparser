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

import types
import numpy as np
from exparser.BaseMatrix import BaseMatrix
from exparser import Constants
from rpy2 import robjects
from rpy2.robjects.packages import importr
lme4 = importr('lme4')
utils = importr('utils')
langR = importr('languageR')

class MixedEffectsMatrix(BaseMatrix):

	"""A matrix containing the results of an Anova analysis"""

	def __init__(self, dm, dv, fixedEffects, randomEffects, continuous=[], \
		nsim=10000):

		"""
		Constructor. Note that for now all effects are assumed to be factors,
		and not continuous variables.

		Arguments:
		dm 				--	a DataMatrix
		factors			--	a list of factors
		dv				--	the dependent variable
		fixedEffects		-- 	a list of fixed effects (i.e. the indendendent
							variables)
		randomEffects	--	a list of random effects, such as subject or item

		Keyword arguments:
		continuous		--	a list of effects that should be treated as
							continuous, rather than discrete (default=[])
		nsim				--	the number of simulations to estimate the P values
							(default=10000)
		"""

		if type(dv) == types.FunctionType:
			sdv = dv.__name__
		else:
			sdv = dv

		self.dm = dm

		# Construct a formula for R, which should be of the type
		# 'dv ~ fixed1 * fixed 2  + (1|random1) + (1+random2)'
		fixedEffectsTemplate = '*'.join(fixedEffects)
		randomEffectsTemplate = '+'.join(['(1|%s)' % re for re in randomEffects])
		self._formula = '%s ~ %s + %s' % (dv, fixedEffectsTemplate, \
			randomEffectsTemplate)

		# Register all the variables so that R can use them
		for f in fixedEffects + randomEffects:
			if f in continuous:
				robjects.globalenv[f] = robjects.FloatVector(dm[f])
			else:
				robjects.globalenv[f] = robjects.FactorVector(dm[f])
		robjects.globalenv[sdv] = robjects.FloatVector(dm[dv])

		# Perform the regression
		model = lme4.lmer(robjects.Formula(self._formula), verbose=False)

		# Estimate p-values
		self.pVals = langR.pvals_fnc(model, nsim=nsim)

		## Convert the results into a Matrix for easy reading
		l = [['f0', 'slope', 'Pr(>|t|)', 'sign.']]
		r = 0
		for row in self.pVals[0].rownames:
			slope = self.pVals[0][0][r]
			p = float(self.pVals[0][5][r])
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
			r += 1
			l.append([row, slope, p, sign])
		self.m = np.array(l, dtype='|S128')

	def formula(self):

		"""
		Return the formula (model) that was used as input for the R aov()
		function

		Returns:
		A formula
		"""

		return self._formula




