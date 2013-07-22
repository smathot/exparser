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
from scipy.stats import nanmean

def cowansK(dm):

	"""
	Calculates Cowan's K, defined as

		(hit-rate + correct-rejection-rate - 1) * N

	Arguments:
	dm -- a DataMatrix

	Returns:
	Cowan's K
	"""

	changeDm = dm.select('targetChange == "yes"', verbose=False)
	noChangeDm = dm.select('targetChange == "no"', verbose=False)
	hitRate = np.mean(changeDm['correct'])
	corRejRate = np.mean(noChangeDm['correct'])
	N = np.mean(dm['N'])
	return (hitRate + corRejRate - 1) * N

def dPrime(dm):

	"""
	Calculates d', defined as

		Z(hit-rate) - Z(false-alarm-rate)

	Arguments:
	dm -- a DataMatrix

	Returns:
	d'
	"""

	changeDm = dm.select('targetChange == "yes"', verbose=False)
	noChangeDm = dm.select('targetChange == "no"', verbose=False)
	HR = np.mean(changeDm['correct'])
	FAR = np.mean(1-noChangeDm['correct'])
	return norm.ppf(HR) - norm.ppf(FAR)

def rtAdjust(dm, rtVar='response_time', correctVar='correct'):

	"""
	Calculates the adjusted RT, which is the response time of the correct trials
	divided by the accuracy.

	Arguments:
	dm			--	A DataMatrix.

	Keyword arguments:
	rtVar		--	The variable that indicates the response time. (default=
					'response_time')
	correctVar	--	The variable that indicates the correctness. Should have the
					values 0 and 1. (default='correct')

	Returns:
	The adjusted RT
	"""

	dmCor = dm.select('%s == 1' % correctVar, verbose=False)
	rt = nanmean(dmCor[rtVar])
	acc = nanmean(dm[correctVar])
	return 1.*rt / acc