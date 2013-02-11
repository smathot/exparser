#-*- coding:utf-8 -*-

butter = ['#fce94f', '#edd400', '#c4a000']
orange = ['#fcaf3e', '#f57900', '#ce5c00']
chocolate = ['#e9b96e', '#c17d11', '#8f5902']
chameleon = ['#8ae234', '#73d216', '#4e9a06']
skyBlue = ['#729fcf', '#3465a4', '#204a87']
plum = ['#ad7fa8', '#75507b', '#5c3566']
scarletRed = ['#ef2929', '#cc0000', '#a40000']
aluminium = ['#eeeeec', '#d3d7cf', '#babdb6', '#888a85', '#555753', '#2e3436']

allColors = butter + orange + chocolate + chameleon + skyBlue + plum + \
	scarletRed
	
brightColors = [butter[0], orange[0], chameleon[0], skyBlue[0], plum[0], \
	scarletRed[0]]

yellow = butter
brown = chocolate
green = chameleon
blue = skyBlue
purple = plum
red = scarletRed
grey = aluminium
gray = aluminium

def rgb(col, r=1.):	
	if type(col) != str or len(col) != 7 or col[0] != '#':
		raise Exception('Expecting an HTML color, not "%s"' % col)		
	r = (int(col[1:3], 16)/255.)*r
	g = (int(col[3:5], 16)/255.)*r
	b = (int(col[5:7], 16)/255.)*r
	return r,g,b
	
	
