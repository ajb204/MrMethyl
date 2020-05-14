# ProtSkin - protein visualization tool
# Copyright (C) 2003-2004 Christophe Deprez
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation (version 2)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
import sys, os, re, string
from Tkinter import *
from BLOSUM import blosum
from Sequence import Sequence, Pattern

omitted = []

def writeout(string): sys.stdout.write(string)
def message(string=''): 
	if re.search(r'^win',sys.platform): 
		raw_input(string)
	else: 
		raise StandardError, string

def HTMLcolor(value=omitted, color=omitted):
		if value is omitted: return "#ffffff"
		halfvalue = 128 + 0.5 * value
		if color == 1: return "#%02x%02xff" % (value,value)	# blue
		if color == 2: return "#%02xff%02x" % (value,value)	# green
		if color == 3: return "#%02xffff" % value		# cyan
		if color == 4: return "#ff%02x%02x" % (value,value)	# red
		if color == 5: return "#ff%02xff" % value		# magenta
		if color == 6: return "#ffff%02x" % value		# yellow
		if color == 7: return "#ff%02x%02x" % (halfvalue,value)	# orange
		return "#%02x%02x%02x" % (halfvalue,halfvalue,halfvalue)# black
		
class TkSequence(Sequence):
	def __init__(self, sequence):
		" Create sequence object (copied from first argument) with Tk methods "
		Sequence.__init__(self)
		self.START = sequence.START
		for i in range (sequence.LENGTH):
			self.append(sequence.SEQ[i],sequence.VAL[i])
		
	def printtk(self, canvas, x, y, favorite_color=0, invert_colors=0, header=''):
		" Print sequence-value pairs in Tkinter format "
		if self.min() < 0 or self.max() > 1:
			color = self.copy().reduce()
		else :  color = self.copy().reduce(0.0,1.0)
				
		if invert_colors: 	color.flin(-1,1).flin(255)
		else:			color.flin(255)
		
		canvas.create_rectangle(x, y, (x+50), (y+20), fill='white')
		canvas.create_text(x+25,y+10, text=header)
		y += 20
		for i in range(self.LENGTH):
			canvas.create_rectangle(x, y+20*i, (x+20), (y+20)+20*i, fill=(HTMLcolor(color.VAL[i],favorite_color)))
			canvas.create_text(x+10,y+10+20*i, text=self.SEQ[i])
			canvas.create_rectangle(x+20, y+20*i, (x+50), (y+20)+20*i, fill=(HTMLcolor(color.VAL[i],favorite_color)))
			canvas.create_text(x+35,y+10+20*i, text="%.2f"%self.VAL[i])
		canvas.config(scrollregion=(0,0,x+50+10,y+20+20*i+10))
		
		#   x    20           30          
		# y +----------+---------------+
		#   |          |               |
		# 20|     +    |       +       |
		#   |          |               |
		#   +----------+---------------+


class TkPattern(Pattern):
	def __init__(self, pattern):
		" Create pattern object (copied from first argument) with Tk methods "
		self.DICT = {} # Format: { 'A':[0 0 0 ...], 'B':[0 0...] ...}
		for k in pattern.DICT.keys():
			self.DICT[k] = pattern.DICT[k]
						
	def printtk(self, canvas, x, y, favorite_color=0, invert_colors=0, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
		" Print pattern in Tkinter format "
		j = 0
		for k in range(len(keystring)):
			canvas.create_rectangle(x+30*k, y+20*j, (x+30)+30*k, (y+20)+20*j, fill='white')
			canvas.create_text((x+15)+30*k,(y+10)+20*j, text=keystring[k])
		y+=20
		for j in range(self.length()):
			for k in range(len(keystring)):
				if invert_colors: color = 255*self.DICT[keystring[k]][j]
				else: color = 255*(1.0-self.DICT[keystring[k]][j])
				canvas.create_rectangle(x+30*k, y+20*j, (x+30)+30*k, (y+20)+20*j, fill=HTMLcolor(color,favorite_color))
				canvas.create_text((x+15)+30*k,(y+10)+20*j, text="%.2f"%self.DICT[keystring[k]][j])
		canvas.config(scrollregion=(0,0,x+30+30*k+10,y+20+20*j+10))
			
		#   x         30                   30          
		# y +--------------------+--------------------+
		#   |                    |                    |
		# 20|          +         |          +         |
		#   |                    |                    |
		#   +--------------------+--------------------+
		
