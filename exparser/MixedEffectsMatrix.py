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

class MixedEffectsMatrix(BaseMatrix):

	"""A matrix containing the results of an Anova analysis"""

	def __init__(self, dm, dv, fixedEffects, randomEffects, continuous=[], \
		nSim=10000, formula=None, fixedOp='*'):

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
							continuous, rather than discrete. Instead of a list,
							you can also pass the string 'all' the treat all
							effects as continuous (default=[])
		nSim			--	the number of simulations to estimate the P values
							(default=10000)
		formula			--	use a specific formular for R, or None to
							auto-generate (default=None)
		fixedOp			--	the operator to use for the fixed effects
							(default='*')
		pVals			--	indicates whether p-values should be estimated
							(default=True)
		"""

		from rpy2 import robjects
		from rpy2.robjects.packages import importr
		lme4 = importr('lme4')
		utils = importr('utils')
		langR = importr('languageR')

		if type(dv) == types.FunctionType:
			sdv = dv.__name__
		else:
			sdv = dv

		self.dm = dm

		# Construct a formula for R, which should be of the type
		# 'dv ~ fixed1 * fixed 2 + (1|random1) + (1|random2)'
		# Random slopes example:
		# 'dv ~ fixed1 * fixed 2 + (1+fixed1|subject) + (1+random2)'
		
		if formula == None:
			fixedEffectsTemplate = fixedOp.join(fixedEffects)
			randomEffectsTemplate = '+'.join(['(1|%s)' % re for re in \
				randomEffects])
			self._formula = '%s ~ %s + %s' % (dv, fixedEffectsTemplate, \
				randomEffectsTemplate)
		else:
			self._formula = formula

		# Register all the variables so that R can use them
		for f in fixedEffects + randomEffects:
			if continuous == 'all' or f in continuous:
				robjects.globalenv[f] = robjects.FloatVector(dm[f])
			else:
				robjects.globalenv[f] = robjects.FactorVector(dm[f])
		robjects.globalenv[sdv] = robjects.FloatVector(dm[dv])

		# Perform the regression
		self.model = lme4.lmer(robjects.Formula(self._formula), verbose=False)
		
		# Get the start of the model output that has t-values and standard
		# errors
		lModel = str(self.model).split('\n')
		for i in range(len(lModel)):
			if lModel[i] == 'Fixed effects:':
				break
		i += 2				
		self.pVals = langR.pvals_fnc(self.model, nsim=nSim, ndigits=10)

		## Convert the results into a Matrix for easy reading
		l = [['f0', 'slope', 'SE', 'ci95low', 'ci95high', 't', 'Pr(>|t|)', \
			'sign.']]
		r = 0
		for row in self.pVals[0].rownames:
			slope = float(self.pVals[0][0][r])
			ci95low = float(self.pVals[0][2][r])
			ci95high = float(self.pVals[0][3][r])
			print lModel[i+r]
			t = float(lModel[i+r].split()[3])
			se = float(lModel[i+r].split()[2])
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
			l.append([row, slope, se, ci95low, ci95high, t, p, sign])
		self.m = np.array(l, dtype='|S128')

	def formula(self):

		"""
		Return the formula (model) that was used as input for the R aov()
		function

		Returns:
		A formula
		"""

		return self._formula




