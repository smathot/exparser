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

def inspect(dm):
	
	from PyQt4 import QtGui, QtCore
	app = QtGui.QApplication([])
	table = QtGui.QTableWidget(len(dm), len(dm.columns()))
	for row in dm.range():
		for col, name in enumerate(sorted(dm.columns())):
			if row == 0:
				item = QtGui.QTableWidgetItem(unicode(name))
				table.setHorizontalHeaderItem(col, item)
			val = dm[name][row]
			item = QtGui.QTableWidgetItem(unicode(val))
			table.setItem(row, col, item)
	table.show()
	app.exec_()
	