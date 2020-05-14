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
from BLOSUM import blosum

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

class Sequence:
	def __init__(self):
		" Create new void sequence "
		self.START = 1
		self.SEQ = []	
		self.LENGTH = 0
		self.VAL = []
#		self.seq = ['A', 'B', 'C']
#		self.val = [1, 2, 3]
		
	def copy(self):
		copy = Sequence()
		copy.START = self.START
		for i in range (self.LENGTH):
			copy.append(self.SEQ[i],self.VAL[i])
		return copy	
	
	def start(self, start=omitted):
		if start is not omitted:
			self.START = start
		return self.START	
	
	def length(self):
		return self.LENGTH	

	def seq(self, sequence_number, one_letter_code=omitted):
		" Put/Get residue sequence "
		i = sequence_number - self.START
		if one_letter_code is not omitted:
			self.SEQ[i] = one_letter_code
		return self.SEQ[i]
		
	def val(self, sequence_number, value=omitted):
		" Put/Get value "
		i = sequence_number - self.START
		if value is not omitted:
			self.VAL[i] = value
		return self.VAL[i]

	def inc(self, sequence_number):
		" Increment value"
		i = sequence_number - self.START
		self.VAL[i] += 1
		return self.VAL[i]
		
	def append(self, residue, value=0):
		" Return new length " 
		self.SEQ += [residue]
		self.LENGTH += 1
		self.VAL += [value]
		return self.LENGTH
		
	def equalseq(self, alter):
		" Returns 1 if the sequences are equal "
		" (including numbering but not values _> use equal for full comparison) "
		if self.LENGTH != alter.LENGTH: return 0
		if self.START != alter.START: return 0
		for i in range(self.LENGTH):
			if self.SEQ[i] != alter.SEQ[i]: return 0
		return 1
		
	def equal(self, alter):
		" Returns 1 if the sequences are equal "
		" (including numbering and values) "
		if self.equalseq(alter) == 0: return 0
		for i in range(self.LENGTH):
			if self.VAL[i] != alter.VAL[i]:	return 0
		return 1

	def min(self):
		mini = self.VAL[0]
		for i in range(self.LENGTH):
			if self.VAL[i] < mini: mini = self.VAL[i]
		return mini
		
	def max(self):
		maxi = self.VAL[0]
		for i in range(self.LENGTH):
			if self.VAL[i] > maxi: maxi = self.VAL[i]
		return maxi
		
	def flin(self, a=1, b=0):
		" Submit values to a linear function "
		for i in range(self.LENGTH):
			self.VAL[i] = a * self.VAL[i] + b
		return self
		
	def string(self):
		" Return sequence as a string "
		mylist = ''
		for i in range(self.LENGTH):
			mylist += self.SEQ[i]
		return mylist

        def printseq(self, SEPARATOR=omitted, FILE=omitted):
                " Print one-letter-code sequence "
                if FILE is omitted: HANDLE = sys.stdout
                else:
                        try: HANDLE = open(FILE, 'w')
                        except IOError:
                                writeout("Can't open "+FILE+"\n")
                                message()
                for i in range(self.LENGTH):
                    if SEPARATOR is omitted:    HANDLE.write(self.SEQ[i]+' ')
                    else:                       HANDLE.write(self.SEQ[i]+SEPARATOR)
                HANDLE.close
		
	def printtext(self, FILE=omitted):
		" Print sequence-value pairs to file "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		for i in range(self.LENGTH):
			HANDLE.write("%s %5.2f\n" % (self.SEQ[i],self.VAL[i]))
		HANDLE.close
		
	def printhtml(self, HANDLE=sys.stdout, favorite_color=0, invert_colors=0, header='<BR>'):
		" Print sequence-value pairs in HTML format  "
		if self.min() < 0 or self.max() > 1:
			color = self.copy().reduce()
		else :  color = self.copy().reduce(0.0,1.0)
				
		if invert_colors: 	color.flin(-1,1).flin(255)
		else:			color.flin(255)
		
		HANDLE.write("<TABLE BORDER=1>")
		HANDLE.write("<TR><TD COLSPAN=2 align=center>"+header+"</TD></TR>")
		for i in range(self.LENGTH):
		    HANDLE.write("<TR bgcolor = \"%s\"><TD align=center>%s</TD><TD align=center>%.2f</TD></TR>" % (HTMLcolor(color.VAL[i],favorite_color),self.SEQ[i],self.VAL[i] ))
		HANDLE.write("</TABLE>")		
		
	def printvalhtml(self, HANDLE=sys.stdout, favorite_color=0, invert_colors=0, header='<BR>'):
		" Print values in HTML format  "
		if self.min() < 0 or self.max() > 1:
			color = self.copy().reduce()
		else :  color = self.copy().reduce(0.0,1.0)
				
		if invert_colors: 	color.flin(-1,1).flin(255)
		else:			color.flin(255)
		
		HANDLE.write("<TABLE BORDER=1>")
		HANDLE.write("<TR><TD COLSPAN=2 align=center>"+header+"</TD></TR>")
		for i in range(self.LENGTH):
		    HANDLE.write("<TR bgcolor = \"%s\"><TD align=center>%.2f</TD></TR>" % (HTMLcolor(color.VAL[i],favorite_color),self.VAL[i] ))
		HANDLE.write("</TABLE>")
		
	def reduce(self, mini=omitted, maxi=omitted):
		" Changes values linearly in the range [0;1] instead of [1starg;2ndarg] "
		" If no argument is given, min and max values are used to biject with 0 and 1 "
		if mini is omitted: mini = self.min()
		if maxi is omitted: maxi = self.max()
		if mini == maxi: maxi = mini + 1			
		return self.flin(1.0/(mini-maxi),(0.0+maxi)/(maxi-mini))
		
	def valfrompattern(self, pattern):
		" for each residue : retrieve val from pattern "
		if self.length() != pattern.length():
			message('Error using method : valfrompattern')
		for i in range(self.LENGTH):
			self.VAL[i] = pattern.val(self.SEQ[i],i)
						
	def valfrompatternsum(self, pattern, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
		" for each residue : get val from pattern sum "
		if self.length() != pattern.length():
			message('Error using method : valfrompatternsum')
		for i in range(self.length()):
			self.VAL[i] = 0
			for k in keystring:
				self.VAL[i] += pattern.DICT[k][i]
		return self
				
	def divide(self, divider):
		for i in range(self.length()):
			self.VAL[i] /= float(divider)
		
	def bad_spelling(self, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
		" return invalid character if sequence does not use keystring "
		" return 0 if everything's fine "
		for i in range(self.length()):
			if not re.search(self.SEQ[i],keystring):
				return self.SEQ[i]
		return 0
	
	def is_nuc(self):
		for i in range (self.LENGTH):
			if (self.SEQ[i] == "U"): return 1
		return 0	
	
class Sequence2(Sequence):
	def __init__(self, sequence=omitted):
		" Create new sequence or copy it from first argument "
		Sequence.__init__(self)
		self.VAL2 = []
		if sequence is not omitted:
			self.START = sequence.START
			for i in range (sequence.LENGTH):
				self.append(sequence.SEQ[i],sequence.VAL[i])
		
	def append(self, residue, value=0, value2=0):
		" Append optional values and return new length " 
		Sequence.append(self,residue, value)
		self.VAL2 += [value2]
		return self.LENGTH
		
	def copy(self):
		" Create a copy from same class "
		copy = Sequence2()
		copy.START = self.START
		for i in range (self.LENGTH):
			copy.append(self.SEQ[i],self.VAL[i],self.VAL2[i])
		return copy
		
	def printtext(self, FILE=omitted):
		" Print sequence-value pairs to file "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		for i in range(self.LENGTH):
			HANDLE.write("%s %5.2f %5.2f\n" % (self.SEQ[i],self.VAL[i],self.VAL2[i]))
		HANDLE.close


	def printhtml(self, HANDLE=sys.stdout, favorite_color=0, invert_colors=0, header='<BR>'):
		" Print sequence-value pairs in HTML format  "
		if self.min() < 0 or self.max() > 1:
			color = self.copy().reduce()
		else :  color = self.copy().reduce(0.0,1.0)
				
		if invert_colors: 	color.flin(-1,1).flin(255)
		else:			color.flin(255)
		
		HANDLE.write("<TABLE BORDER=1>")
		HANDLE.write("<TR><TD COLSPAN=2 align=center>"+header+"</TD></TR>")
		for i in range(self.LENGTH):
		    HANDLE.write("<TR bgcolor = \"%s\"><TD align=center>%d</TD><TD align=center>%s</TD><TD align=center>%.2f</TD></TR>" % (HTMLcolor(color.VAL[i],favorite_color),self.VAL2[i],self.SEQ[i],self.VAL[i] ))
		HANDLE.write("</TABLE>")

class Pattern:
	def __init__(self, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
		" Create new void pattern "
		self.DICT = {} # Format: { 'A':[0 0 0 ...], 'B':[0 0...] ...}
		for k in keystring:
			self.DICT[k] = []
#		print self.DICT
#		print self.DICT.keys()
		
	def fromseq(self, seq):
		" Update pattern from sequence list "
		# For each residue position in sequence(s)
		for j in range(seq[0].length()):
			# Initialize pattern
		    	for k in self.DICT.keys():
		    		self.DICT[k] += [0]
			# Scan sequences and update pattern
			for s in seq:
				self.DICT[s.seq(s.start()+j)][j] += 1
				
	def length(self):
		for k in self.DICT.keys():
			return len(self.DICT[k])	
					
	def val(self, key, index):
		return self.DICT[key][index]
			
	def ordered_keys(self):
		k = self.DICT.keys()
		k.sort()
		return k
	
	def keys(self):
		return self.DICT.keys()
	
	def percent_format(self):
		" Turn values to percentage among keys "
		for j in range(self.length()):
			somme = 0.0
			for k in self.keys(): 
				somme += self.DICT[k][j]
			for k in self.keys(): 
				self.DICT[k][j] /= somme
				
	def blosum62_weight(self, seq, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
		" Weigh pattern according to similarity to sequence given as argument"
		if self.length() != seq.length():
			message('Error using method : similarity_weight')
		for j in range(self.length()):
#			self.DICT[seq.seq(seq.start()+j)][j] -= 1
			for k in keystring:
				self.DICT[k][j] *= blosum(k,seq.seq(seq.start()+j))
		
	def sum(self, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
	# This function is not used any more!
		" Return list of sums "
		returnlist = []
		for j in range(self.length()):
			sum = 0
			for k in keystring:
				sum += self.DICT[k][j]
			returnlist.append(sum)
		return returnlist
			
	def printtext(self):
		print self.DICT
		
	def printhtml(self, HANDLE=sys.stdout, favorite_color=0, invert_colors=0, keystring = 'ACDEFGHIKLMNPQRSTVWY-'):
		" Print pattern in HTML format "
		HANDLE.write("<TABLE BORDER=1><TR>")
		for k in keystring: HANDLE.write("<TD align=center>"+k+"</TD>")	
		for j in range(self.length()):
			HANDLE.write("</TR><TR>")
			for k in keystring: 
				if invert_colors: color = 255*self.DICT[k][j]
				else: color = 255*(1.0-self.DICT[k][j])
				HANDLE.write("<TD align=center bgcolor = \"%s\">%.2f</TD>" % ( HTMLcolor(color,favorite_color), self.DICT[k][j] ))	
		HANDLE.write("</TR></TABLE>")
				
	def consensus(self, keystring = 'ACDEFGHIKLMNPQRSTVWY'):
		" Return consensus sequence from pattern "
		consensus = Sequence()
		for j in range(self.length()):
			keymax = []; omax = 0.0
			for k in keystring:
				o = self.DICT[k][j]
				if o > omax:
					omax = o
					keymax = k
			consensus.append(keymax,omax)
		return consensus


