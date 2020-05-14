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
from Sequence import Sequence, Sequence2, message
omitted = []
INITCOLOR = 0

def writeout(string): sys.stdout.write(string)

class PDB:
	def __init__(self, PDB_lines=omitted):
		" Contains main information from PDB file "
		" Features left over : "
		" 	Alternate location indicator "
		" 	Code for insertion of residues "
		self.ATOMNAME = []
		self.AMINOACID = []
		self.CHAIN = []
		self.RESIDUE = []
		self.X = []
		self.Y = []
		self.Z = []
		self.OCCUPANCY = []
		self.TEMPERATURE = []
		self.SEGMENT = []
		self.ELEMENT = []
		self.CHARGE = []
		self.COLOR = []
		self.ATOMSTART = 1
		self.ATOMSIZE = 0
		
		if PDB_lines is not omitted:
			for l in PDB_lines:
			    if re.search(r'^ATOM', l):
			    	if self.ATOMSIZE == 0: self.ATOMSTART = int(l[6:11])
				#writeout(l)
				self.ATOMNAME.append(l[12:16])
				self.AMINOACID.append(l[17:20])
				self.CHAIN.append(l[21:22])
				self.RESIDUE.append(int(l[22:26]))
				self.X.append(float(l[30:38]))
				self.Y.append(float(l[38:46]))
				self.Z.append(float(l[46:54]))
				self.OCCUPANCY.append(float(l[54:60]))
				self.TEMPERATURE.append(float(l[60:66]))
				m = re.search(r'(\w+)', l[72:76])
				if m: self.SEGMENT.append(m.group(1))
				else: self.SEGMENT.append('')
				m = re.search(r'(\w+)', l[76:78])
				if m: self.ELEMENT.append(m.group(1))
				else: self.ELEMENT.append('')
				m = re.search(r'(\w+)', l[78:80])
				if m: self.CHARGE.append(m.group(1))
				else: self.CHARGE.append('')
				self.COLOR.append(INITCOLOR)
				self.ATOMSIZE += 1					
			    if re.search(r'^ENDMDL', l): break
			    if re.search(r'^TER', l): break
			
	def printtext(self, FILE=omitted):
		" Print raw PDB structure to file given as argument (defaults to stdout) "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		for i in range(self.ATOMSIZE):
			HANDLE.write("ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%+2s%+2s\n" % (self.ATOMSTART+i,self.ATOMNAME[i],self.AMINOACID[i],self.CHAIN[i],self.RESIDUE[i],self.X[i],self.Y[i],self.Z[i],self.OCCUPANCY[i],self.TEMPERATURE[i],self.SEGMENT[i],self.ELEMENT[i],self.CHARGE[i]))
		HANDLE.close
		
	def write_color_as_temperature(self, FILE=omitted):
		" Print raw PDB structure to file given as argument (defaults to stdout), replacing the temperature factor by color "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		for i in range(self.ATOMSIZE):
			HANDLE.write("ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%+2s%+2s\n" % (self.ATOMSTART+i,self.ATOMNAME[i],self.AMINOACID[i],self.CHAIN[i],self.RESIDUE[i],self.X[i],self.Y[i],self.Z[i],self.OCCUPANCY[i],self.COLOR[i],self.SEGMENT[i],self.ELEMENT[i],self.CHARGE[i]))
		HANDLE.close
		
	def colorres(self, target=1, value=INITCOLOR):
		" Color each PDB atom of residue (first argument) with value given as second argument (else 0) "
		" Return 0 if ok, 1 if missing (first) argument, 2 if residue not found "
		return_val = 2
		for i in range(self.ATOMSIZE):
			if self.RESIDUE[i] == target:
				self.COLOR[i] = value
				return_val = 0
		return return_val
	
	def aa(self, residue_number):
		i = 0
		while i < self.ATOMSIZE and self.RESIDUE[i] != residue_number: i += 1
		if i >= self.ATOMSIZE: return "XXX"
		else: return self.AMINOACID[i]
	
	def seq(self):
		" Return sequence2 object containing sequence in 1-letter code from PDB file "
		" and color as first value, residue number as second value "
		" Beware : no guaranty that this sequence is continuous "
		seq = Sequence2()
		seq.start(self.RESIDUE[0])
		current_res = 0
		for i in range(self.ATOMSIZE):	# for each atom
			if self.RESIDUE[i] != current_res:	# next residue
				current_res = self.RESIDUE[i]
				seq.append(code(self.AMINOACID[i]),self.COLOR[i],self.RESIDUE[i])
		return seq
	
        def mismatchseq(self, seq):
                " Return 0 if sequence from PDB matches sequence from Seq object given as second argument"
                " Return the mismatched residue number otherwise"
                current_res = -1;
		start_res = seq.start()
                seq_index = -1
                for i in range(self.ATOMSIZE):  # for each atom
                        if self.RESIDUE[i] != current_res:      # next residue
                                current_res = self.RESIDUE[i]
				# print "<BR>"
				# print current_res
				# print self.AMINOACID[i]
				if start_res > 1:	# skip PDB residue before "starting" residue of seq
					start_res -= 1	# one less residue to skip !
				else:
                                	seq_index += 1
					# print seq_index
					# print seq.seq(seq.start()+seq_index)
					# print code(seq.seq(seq.start()))
                                	if seq_index > seq.length():
                                        	return self.RESIDUE[i-1]
                                	if seq.seq(seq.start()+seq_index) != code(self.AMINOACID[i]):
                                        	return self.RESIDUE[i]
                return 0

        def color(self, seq=omitted):
            " Color each PDB atom with color contained in Seq given as argument "
            " This function does not check residue concordance "
            if seq is omitted:
                for i in range(self.ATOMSIZE): self.COLOR[i] = INITCOLOR
            else:
                current_res = -1;
                start_res = seq.start()
                seq_index = -1
                for i in range(self.ATOMSIZE):  # for each atom
                        if self.RESIDUE[i] != current_res:      # next residue
                                current_res = self.RESIDUE[i]
				if start_res > 1:       # skip PDB residue before "starting" residue of seq
                                        start_res -= 1  # one less residue to skip !
					color = 0	# no color for skipped residue!
                                else:
                                	seq_index += 1
					color = seq.val(seq.start()+seq_index)
                       	self.COLOR[i] = color		# color each atom line

	def write_grasp_list(self, FILE=omitted):
		" Print GRASP property file from colors in PDB structure to file given as argument (defaults to stdout) "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		HANDLE.write("atoms=gproperty1\n")
		for i in range(self.ATOMSIZE):
			HANDLE.write("%.2f\n" % self.COLOR[i])
		HANDLE.close
		
	def list_molmol_color(self, HANDLE=sys.stdout, favorite_color=0):
		" Local macro : Print list of colored atoms to HANDLE (in color given as argument) "
		current_res = 0; 
		for i in range(self.ATOMSIZE):	# for each atom
			if self.RESIDUE[i] != current_res:	# next residue
				current_res = self.RESIDUE[i]
				HANDLE.write("SelectAtom ':%d'\n" % self.RESIDUE[i])
				HANDLE.write(coloratom(self.COLOR[i],favorite_color))

	def write_molmol_worm(self, FILE=omitted, favorite_color=0):
		" Print MOLMOL macro from colors in PDB structure to file given as argument (defaults to stdout) "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		HANDLE.write("SelectBond ''\n")
		HANDLE.write("StyleBond invisible\n")
		HANDLE.write("CalcSecondary\n")
		HANDLE.write("SelectAtom '@CA'\n")
		HANDLE.write("AddRibbon\n")
		HANDLE.write("StyleRibbon round sharp sharp const\n")
		HANDLE.write("RadiusPrim 1.0\n")
		HANDLE.write("PaintRibbon atom_smooth\n")
		HANDLE.write("SelectAtom '@CA'\n")
		HANDLE.write("ColorAtom 1 1 1\n")
		HANDLE.write("Light infinite 0 0 7\n")
		HANDLE.write("Projection perspective\n")
		self.list_molmol_color(HANDLE,favorite_color)
		HANDLE.close

	def write_molmol_surface(self, FILE=omitted, favorite_color=0):
		" Print MOLMOL macro from colors in PDB structure to file given as argument (defaults to stdout) "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()
		HANDLE.write("SelectBond ''\n")
		HANDLE.write("StyleBond invisible\n")
		self.list_molmol_color(HANDLE,favorite_color)
		HANDLE.write("Light point 5 2 15\n")
		HANDLE.write("Projection perspective\n")
		HANDLE.write("SelectAtom ''\n")
		HANDLE.write("AddSurface vdw solvent 1.4 shaded 0\n")
		HANDLE.write("PaintSurface atom 0 1.4 0 0 ''\n")
		HANDLE.close
		
	def write_pymol(self, FILE=omitted, favorite_color=0):
		" Print PYMOL macro from colors in PDB structure to file given as argument (defaults to stdout) "
		if FILE is omitted: HANDLE = sys.stdout
		else:
			try: HANDLE = open(FILE, 'w')
			except IOError:
				writeout("Can't open "+FILE+"\n")
				message()				
		
		current_res = 0; 
		for i in range(self.ATOMSIZE):	# for each atom
			if self.RESIDUE[i] != current_res:	# next residue
				current_res = self.RESIDUE[i]
				HANDLE.write("set_color color%d = " % self.RESIDUE[i])
				HANDLE.write(set_pymol_color(self.COLOR[i],favorite_color))
				HANDLE.write("color color%d, resi %d\n" % (self.RESIDUE[i],self.RESIDUE[i]))
		HANDLE.close

	def switch_proton(self, residue_number, atom_name):
		" Change proton nomenclature at residue and atom given as argument (and mark changed atom with color 1)"
		for i in range(self.ATOMSIZE):
			if ((self.RESIDUE[i] == residue_number) and (re.search(atom_name, self.ATOMNAME[i]) and (self.COLOR[i] == 0))):
				if (re.search(r'1H',atom_name)):
					self.ATOMNAME[i] = re.sub(r'1',r'2',self.ATOMNAME[i],1)
				elif (re.search(r'2H',atom_name)):
					self.ATOMNAME[i] = re.sub(r'2',r'1',self.ATOMNAME[i],1)
				self.COLOR[i] = 1 # 'atom changed' color mark

	def patch1(self):
                " change atom nomenclature from '1HG2' to 'HG21', '2HB ' to 'HB2 ', '1H  ' to 'H1  ', etc "
                for i in range(self.ATOMSIZE):
                        m = re.search(r'^(\d)(\w+)(\s*)',self.ATOMNAME[i])
                        if m: self.ATOMNAME[i] = m.group(2) + m.group(1) + m.group(3)
                                                                                                                           
	def shiftres(self, shift=0):
                " Shift residue number "
                for i in range(self.ATOMSIZE): self.RESIDUE[i] += shift
                                                                                

	
def set_pymol_color(value=omitted, color=omitted):
	if value is omitted: return "[1.0, 1.0, 1.0]\n"
	halfvalue = 0.5 + 0.5 * value
	if color == 1: return "[%.3f, %.3f, 1.0]\n" % (value,value)
	if color == 2: return "[%.3f, 1.0, %.3f]\n" % (value,value)
	if color == 3: return "[%.3f, 1.0, 1.0]\n" % value
	if color == 4: return "[1.0, %.3f, %.3f]\n" % (value,value)
	if color == 5: return "[1.0, %.3f, 1.0]\n" % value
	if color == 6: return "[1.0, 1.0, %.3f]\n" % value
	if color == 7: return "[1.0, %.3f, %.3f]\n" % (halfvalue,value)
	return "[%.3f, %.3f, %.3f]\n" % (value,value,value)

def coloratom(value=omitted, color=omitted):
	if value is omitted: return "ColorAtom 1 1 1\n"
	halfvalue = 0.5 + 0.5 * value
	if color == 1: return "ColorAtom %.3f %.3f 1\n" % (value,value)
	if color == 2: return "ColorAtom %.3f 1 %.3f\n" % (value,value)
	if color == 3: return "ColorAtom %.3f 1 1\n" % value
	if color == 4: return "ColorAtom 1 %.3f %.3f\n" % (value,value)
	if color == 5: return "ColorAtom 1 %.3f 1\n" % value
	if color == 6: return "ColorAtom 1 1 %.3f\n" % value
	if color == 7: return "ColorAtom 1 %.3f %.3f\n" % (halfvalue,value)
	return "ColorAtom %.3f %.3f %.3f\n" % (value,value,value)
	
def code(string):
	if string == "ALA": return "A"
	if string == "CYS": return "C"
	if string == "ASP": return "D"
	if string == "GLU": return "E"
	if string == "PHE": return "F"
	if string == "GLY": return "G"
	if string == "HIS": return "H"
	if string == "ILE": return "I"
	if string == "LYS": return "K"
	if string == "LEU": return "L"
	if string == "MET": return "M"
	if string == "ASN": return "N"
	if string == "PRO": return "P"
	if string == "GLN": return "Q"
	if string == "ARG": return "R"
	if string == "SER": return "S"
	if string == "THR": return "T"
	if string == "VAL": return "V"
	if string == "TRP": return "W"
	if string == "TYR": return "Y"
	# Include code for nucleic acids:
        if string == "A  ": return "A"
        if string == "C  ": return "C"
        if string == "G  ": return "G"
        if string == "U  ": return "U"
	writeout("Fatal error : unknown residue "+string+"\n")
	exit(1)	










