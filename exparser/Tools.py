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

import sys
from exparser.DataMatrix import DataMatrix
from exparser import Cache, GUI
import time
import warnings

def callFunc(dm, mods, func, cachePrefix='autoCache.', redo=False):

	"""
	Calls a single function from a module.

	Arguments:
	dm		--	The DataMatrix to analyze.
	mods	--	A list of modules that may contain the function.
	func	--	The function name.

	Keyword arguments:
	cachePrefix	--	A prefix for the cacheId for cachable functions. The
					function name will be appended. (default='autoCache.')
	redo		--	Indicates whether functions should be redone, even if a
					cache is available. (default=False)

	Returns:
	A DataMatrix.
	"""

	if func[0] == '@':
		func = func[1:]
	found = False
	for mod in mods + [GUI]:
		if hasattr(mod, func):
			t1 = time.time()
			_func = getattr(mod, func)
			if not redo and Cache.isCached(_func):
				cacheId = cachePrefix + func
				print '-> Calling %s.%s() [cacheId=%s]' \
					% (mod.__name__, func, cacheId)
				retVal = _func(dm, cacheId=cacheId)
			else:
				print '-> Calling %s.%s() [uncached]' % (mod.__name__,
					func)
				retVal = _func(dm)
			if isinstance(retVal, DataMatrix):
				print '-> DataMatrix was modified'
				dm = retVal
			print '-> Finished %s.%s() in %.2f s' % (mod.__name__, func,
				time.time()-t1)
			found = True
			break # Break in case the same function occurs in multiple modules
	if not found:
		warnings.warn('Helper function %s does not exist' % func)
	return dm

def analysisLoop(dm, mods=[], pre=[], post=[], full=[],
	cachePrefix='autoCache.'):

	"""
	Executes an analysis loop, in which all functions that specified on the
	command are executed. A function is executed if its name is prefixed by `@`
	and if it is present in one of the helpers modules. Cachable functions are
	cached automatically. If a function returns a DataMatrix, this is used to
	replace the current DataMatrix for the following functions.

	Arguments:
	dm			--	The DataMatrix to analyze.
	mods		--	A module or list of modules that contain the analysis
					functions.
	pre			--	A list of function names that should always be executed
					before the rest.
	post		--	A list of function names that should always be executed
					after the rest.
	full		--	A list of function names that make up the full analysis
					pathway.
	cachePrefix	--	A prefix for the cacheId for cachable functions. The
					function name will be appended. (default='autoCache.')
	"""

	if not isinstance(mods, list):
		mods = [mods]
	if len(mods) == 0:
		raise Exception('No modules specified')
	print('Entering analysisLoop()')
	t0 = time.time()
	for func in pre:
		dm = callFunc(dm, mods, func, cachePrefix=cachePrefix)
	if '@full' in sys.argv:
		print('Running full analysis pathway')
		if len(full) == []:
			raise Exception('No full analysis pathway specified')
		for func in full:
			dm = callFunc(dm, mods, func, cachePrefix=cachePrefix)
	else:
		for func in sys.argv:
			if not func[0] == '@':
				continue
			if ':redo' in func:
				func = func.replace(':redo', '')
				redo = True
			else:
				redo = False
			dm = callFunc(dm, mods, func, cachePrefix=cachePrefix, redo=redo)
	for func in post:
		dm = callFunc(dm, mods, func, cachePrefix=cachePrefix)
	print 'Finished analysisLoop() (%.2f s)' % (time.time() - t0)
