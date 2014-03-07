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
import time
import subprocess
import select
import signal
from exparser.CsvReader import CsvReader
from exparser.DataMatrix import DataMatrix

class RBridge(object):

	"""Provides a very basic two-bridge to R."""

	def __init__(self):

		"""Constructor."""

		self.RProcess = subprocess.Popen(['R', '--vanilla'], \
			stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		self.poll = select.poll()
		self.poll.register(self.RProcess.stdout, select.POLLIN)

	def __del__(self):

		"""Destructor."""

		self.close()

	def aov(self, formula, aovVar='aov1'):

		"""
		Performs an ANOVA, using the `aov()` function.

		Arguments:
		formula		--	An R-style mixed-effects formula.

		Keyword arguments:
		aovVar		--	The R variable to store the `lmer()` output in.
						(default='lmer1')

		Returns:
		A DataMatrix with the `aov()` output as summarized by `summary()`.
		"""

		if os.path.exists('.rbridge-aov.txt'):
			os.remove('.rbridge-aov.txt')
		# Peform aov
		s = self.call('(%s <- aov(%s))' % (aovVar, formula))
		# Summarize the data
		self.write( \
			'capture.output(summary(%s)$"Error: Within", file=".rbridge-aov.txt")' \
			% aovVar)
		while not os.path.exists('.rbridge-aov.txt'):
			time.sleep(.1)
		# The output is saved in an ugly ad-hock format, and we need to parse it
		# line by line into a readable format for DataMatrix.
		s = open('.rbridge-aov.txt').read()
		l = [['effect', 'Df', 'SumSq', 'MeanSq', 'F', 'p', 'sign']]
		for row in s.split('\n')[1:-1]:
			# Do not include signficance ratings, which come after the `---`
			# separator.
			if row == '---':
				break
			l.append(row.split())
		dm = DataMatrix(l)
		return dm

	def anova(self, model1, model2, anovaVar='aov1'):

		"""
		Performs a model comparison using the R `anova()` function. The models
		should be fitted and given a name first, using `RBridge.lmer()`.

		Arguments:
		model1		--	The name of the first model.
		model2		--	The name of the first model.

		Keyword arguments:
		anovaVar	--	The R variable to store the `anova()` output in.
						(default='aov1')

		Returns:
		A DataMatrix with the output for the model comparison.
		"""

		if os.path.exists('.rbridge-anova.csv'):
			os.remove('.rbridge-anova.csv')
		self.write( \
			'%s <- anova(%s, %s)' % (anovaVar, model1, model2))
		self.call('write.csv(%s, ".rbridge-anova.csv")' % anovaVar)
		while not os.path.exists('.rbridge-anova.csv'):
			time.sleep(.1)
		# Try this a few times, because sometimes the csv hasn't been written
		# yet
		for i in range(10):
			try:
				dm = CsvReader('.rbridge-anova.csv').dataMatrix()
				break
			except:
				time.sleep(1)
		dm.rename('f0', 'model')
		return dm

	def call(self, cmd):

		"""
		Sends a command to R and returns the resulting output.

		Arguments:
		cmd		--	The R command to send.

		Returns:
		A string of output.
		"""

		self.write(cmd, cat=True)
		return self.read()

	def close(self):

		"""Burns the R bridge."""

		self.write('quit()')
		self.RProcess.wait()

	def lmer(self, formula, lmerVar='lmer1', pvalsVar='pv1', nsim=1000, \
		full=False, pvals=True):

		"""
		Performs a linear mixed-effects model, using the `lmer()` function from
		`lme4`. P-values are calculated to using `pvals.fnc()` from `languageR`.

		Notes for Ubuntu 13.10:

		- Make sure that you use the official R packages from the Ubuntu
		  repository, and not the (newer) packages from the unofficial R PPA.
		- Download languageR 1.4 from the website and install it manually. This
		  is necessary, because the package will fail to install through CRAN.

		Arguments:
		formula		--	An R-style mixed-effects formula.

		Keyword arguments:
		lmerVar		--	The R variable to store the `lmer()` output in.
						(default='lmer1')
		pvalsVar	--	The R variable to store the `pvals.fnc()` output in.
						(default='pv1')
		full		--	Indicates whether the full model output should be
						returned, or only a subset with renamed variable names.
						(default=False)
		pvals		--	Indicates whether the lmer should be fed into
						`pvals.fnc()` to estimate p-values and confidence
						intervals. If not, the results are summarized using
						`getSummary.mer()` from the `memisc` package.
						(default=True)
		Returns:
		A DataMatrix with the model output.
		"""

		if os.path.exists('.rbridge-lmer.csv'):
			os.remove('.rbridge-lmer.csv')

		# Peform lmer
		self.write('library(lme4)')
		s = self.call('(%s <- lmer(%s))' % (lmerVar, formula))

		# Summarize the data using `memisc`
		self.write('library(memisc)')
		self.write( \
			'write.csv(getSummary.mer(%s)$coef, ".rbridge-lmer.csv")' % \
			lmerVar)
		while not os.path.exists('.rbridge-lmer.csv'):
			time.sleep(.1)
		# Try this a few times, because sometimes the csv hasn't been written
		# yet
		for i in range(10):
			try:
				dm = CsvReader('.rbridge-lmer.csv').dataMatrix()
				break
			except:
				time.sleep(1)
		# If we don't want to use the pvals estimation, this is all we use.
		if not pvals:
			if not full:
				dm = dm.selectColumns(['f0', 'stat', 'est', 'se', 'stat', 'p', \
					'lwr', 'upr'])
				dm.rename('f0', 'effect')
				dm.rename('lwr', 'ci95lo')
				dm.rename('upr', 'ci95up')
				dm.rename('stat', 't')
			return dm
		# If we do want to use the pvals estimation, we need to extract the
		# column with t-values, so that we can add it to the pvals output, which
		# doesn't contain this information. The last two rows are dummy rows, so
		# we strip them.
		os.remove('.rbridge-lmer.csv')
		tvals = dm['stat'][:-2]
		self.write('library(languageR)')
		self.write( \
			'%s <- pvals.fnc(%s, nsim=%d, ndigits=20)' % \
			(pvalsVar, lmerVar, nsim))
		self.write('write.csv(%s$fixed, ".rbridge-lmer.csv")' % pvalsVar)
		while not os.path.exists('.rbridge-lmer.csv'):
			time.sleep(.1)
		# Try this a few times, because sometimes the csv hasn't been written
		# yet
		for i in range(10):
			try:
				dm = CsvReader('.rbridge-lmer.csv').dataMatrix()
				break
			except:
				time.sleep(1)
		if not full:
			dm = dm.selectColumns(['f0', 'Estimate', 'HPD95lower', \
				'HPD95upper', 'Pr(>|t|)'])
			dm.rename('f0', 'effect')
			dm.rename('Estimate', 'est')
			dm.rename('HPD95lower', 'ci95lo')
			dm.rename('HPD95upper', 'ci95up')
			dm.rename('Pr(>|t|)', 'p')
			dm = dm.addField('t', dtype=float)
			dm['t'] = tvals
		return dm

	def load(self, dm, frame='data'):

		"""
		Attaches a DataMatrix as a dataframe to R.

		Arguments:
		dm		--	A DataMatrix.

		Keyword arguments:
		frame	--	The variable name for the dataframe in R. (default=data)
		"""

		dm.save('.rbridge.csv')
		self.write('%s <- read.csv(".rbridge.csv")' % frame)
		self.write('attach(%s)' % frame)

	def read(self):

		"""
		Reads pending R output.

		Returns:
		A string of output.
		"""

		l = []
		signal.signal(signal.SIGALRM, self.timeout)
		while True:
			signal.alarm(5)
			try:
				s = self.RProcess.stdout.readline()
				if not s.strip().startswith('>'):
					l.append(s)
			except:
				break
			if len(l) > 0 and l[-1].strip() == 'end':
				break
		signal.alarm(0)
		return ''.join(l)

	def timeout(self):

		raise Exception('A timeout occurred')

	def write(self, cmd, cat=False):

		"""
		Sends input to R.

		Arguments:
		cmd		--	The R command to send.

		Keyword arguments:
		cat		--	Indicates whether a special 'end' tag should be appended.
					This is mostly for internal use. (default=False)
		"""

		_cmd = '\n%s\n' % cmd
		if cat:
			_cmd += 'cat("end\\n")\n'
		self.RProcess.stdin.write(_cmd)
