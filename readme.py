#!/usr/bin/env python
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

from yamldoc import DocFactory
from exparser.DataMatrix import DataMatrix
from academicmarkdown import build

readme = u"""
# Exparser

A Python library for the analysis of experimental data.

Copyright 2011-2014 Sebastiaan Mathôt

## Overview

%%--
toc:
 mindepth: 1
 exclude: [Exparser, Overview]
--%%

%(DataMatrix)s

""" % {
	'DataMatrix' : unicode(DocFactory(DataMatrix))
	}

build.MD(readme, 'readme.md')
build.PDF(readme, 'readme.pdf')
