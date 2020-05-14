#!/usr/bin/python
#
# ProtSkin - protein visualization tool
# Copyright (C) 2003-2006 Christophe Deprez
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
#################################################################################
# 	 			PROTSKIN					#
#    CREATE GRASP PYMOL MOLMOL INPUT FROM BLAST ALIGNMENT OF PROTEIN SEQUENCE	#
#	Version : 1.82 - Date : Aug 31, 2009 - Author : Christophe Deprez 	#
#################################################################################

VERSION = '1.82'

import sys, os, re, string, cgi, time

from Sequence import Sequence, Sequence2, Pattern, message
from PDB import PDB
from BLOSUM import blosum
from BLAST import parse_BLAST

omitted = []

def writeout(string): sys.stdout.write(string)

def man(HANDLE=sys.stdout): HANDLE.write("""
	NAME :
ProtSkin - protein visualization tool

	SYNOPSIS :
ProtSkin [-b -u -s -i -1 -2 -3 -4 -5 -6 -7] filename [pdbfilename]

        DESCRIPTION :
ProtSkin is a protein visualization tool that create scripts for different programs to use with a PDB file.
The scripts allow the coloring of the protein surface or worm, colored and shaded using numerical
per residue values. Scripts for MOLMOL and PyMOL as generated, as well as a GRASP property file 
(containing one value per atom in the PDB file) and a pseudo-PDB file (where the temperature factor 
is replaced by the desired score) to use with programs such as InsightII.

ProtSkin can be used in 2 ways :
1) You can use ProtSkin to map any scalar data, such as heteronuclear NOE or chemical shift differences 
onto your protein structure. 
Prepare a plain text file (first column : residue numbers, second column : values to map)
and give the file's name as first argument to ProtSkin. Give the PDB file name as second argument.

2) You can use ProtSkin to read a BLAST, CLUSTAL or MSF alignment and create scripts to visualize residue conservation scores 
onto the structure. Scores correspond to conservation percentage of the corresponding residue in the alignment.
The BLAST alignment should be saved as "flat query anchored with letters for identities" into a file.
Use the -b option and give the file's name as first argument to ProtSkin. Give your PDB file name as second 
argument (if no PDB file name is provided here, conservation scores are calculated but no script is generated).
Insertions in the query sequence are discarded. Lowercase letters are converted to uppercase.
To allow the script to check the residue correspondance between the query sequence and the PDB file, ensure that 
the number of the first residue of the query sequence in the BLAST file corresponds to its number in the PDB file.
Two types of conservations scores are generated simultaneously, corresponding to conservation to the query sequence 
or to conservation to the consensus sequence (the latter scripts start with the letter c for consensus).
	
        OPTIONS :
-b   Input is a BLAST, CLUSTAL or MSF file (otherwise, it is a plain list of scores as described above)
-u   Redundant sequences in the alignment are discarded.
-s   Similarity to the query sequence is calculated for the score (instead of identity).
-i   Invert colors (higher score = less color)
-1   Color is blue
-2   Color is green
-3   Color is cyan
-4   Color is red
-5   Color is magenta
-6   Color is yellow
-7   Color is orange

	EXAMPLES :
1) Map scalar data (in red) : 
ProtSkin -4 mydatafile mypdbfile

2) Map BLAST (or CLUSTAL or MSF) alignment discarding redundant sequences and coloring the less conserved residues (in orange).
ProtSkin -b -u -i -7 myblastfile mypdbfile

        REFERENCES :
http://www.ncbi.nlm.nih.gov/BLAST/
http://npsa-pbil.ibcp.fr/cgi-bin/npsa_automat.pl?page=npsa_clustalw.html
http://trantor.bioc.columbia.edu/grasp/
http://www.mol.biol.ethz.ch/wuthrich/software/molmol/
http://pymol.sourceforge.net/
http://www.accelrys.com/insight/index.html

""")


NB_OPT = 6

CGI = 0 # ------ STANDALONE NO COMMENT HERE ------ #

#################################################################################
# 			PARSE AND CHECK ARGUMENTS				#
#################################################################################

if CGI: pass # ------ STANDALONE NO COMMENT HERE ------ #
else: # Standalone Version
   	try:
		# Tkinter interface
		from Tkinter import *
		from TkInterface import Interface
		from TkSequence import TkSequence, TkPattern
	
		INTERFACE = Interface()
	
		if (INTERFACE.is_ok()): ARGUMENTS = INTERFACE.getarg()
		else: sys.exit()

		if re.search(r'^win',sys.platform): 	FILE_PREFIX = 'Results\\'
		else: 					FILE_PREFIX = ''

		TKINTER = 1
   	except (ImportError, TclError):
   		# Command line interface
   		TKINTER = 0
   	
		# Get options
		if re.search(r'^win',sys.platform):
			### SIMULATED command line arguments for WINDOWS platform
			print "Please enter ProtSkin command line arguments and press RETURN :"
			ARGUMENTS = re.split(r'\s+',raw_input('ProtSkin '))	
			if ARGUMENTS[0] == '': ARGUMENTS.pop(0)
			FILE_PREFIX = ''
			FILE_PREFIX = 'Results\\'
		else:
			### TRUE command line arguments
			ARGUMENTS = sys.argv
			ARGUMENTS.pop(0)
			FILE_PREFIX = ''
	
#	print ARGUMENTS
	
	OPT = ''
	while len(ARGUMENTS)>0 and re.search(r'^-\w*',ARGUMENTS[0]):
		OPT += re.search(r'^-(\w*)',ARGUMENTS[0]).group(1)
		del ARGUMENTS[0]
		
	# Get files and check for void files
	if len(ARGUMENTS) < 2 or ARGUMENTS[1] == '':
	     PDB_EXIST = 0
	     if len(ARGUMENTS) < 1 or ARGUMENTS[0] == '': 
		man(sys.stdout)
		message('Missing argument')
	else: PDB_EXIST = 1
		
	try: HANDLE = open(ARGUMENTS[0], 'rU')
  	except IOError:
		message("Can't open file : "+ARGUMENTS[0]+"\n")
  	BLAST_lines = HANDLE.readlines()
  	HANDLE.close()

	if PDB_EXIST:
	  try: HANDLE = open(ARGUMENTS[1], 'rU')
  	  except IOError:
		message("Can't open file : "+ARGUMENTS[1]+"\n")
  	  PDB_lines = HANDLE.readlines()
  	  HANDLE.close()


 	BLAST_EXIST = 1 # The no-BLAST feature is currently implemented in CGI only
	# !!! MacOS uses /r as a newline, which is not recognized by the readlines function 
	if len(BLAST_lines) == 1: BLAST_lines = re.split('\r',BLAST_lines[0])
	if PDB_EXIST and len(PDB_lines) == 1: PDB_lines = re.split('\r',PDB_lines[0])

ALIGNMENT_FILE = FILE_PREFIX + "sequence_alignment";
CONSERVATION_FILE = FILE_PREFIX + "residue_conservation";	
GRASP_FILE = FILE_PREFIX + "grasp_file";
MOLMOLWORM_FILE = FILE_PREFIX + "molmol-worm.mac";	
MOLMOLSURFACE_FILE = FILE_PREFIX + "molmol-surface.mac";	
PYMOL_FILE = FILE_PREFIX + "pymol.pml";
TEMPERATURE_FILE = FILE_PREFIX + "temperature.pdb";	
CGRASP_FILE = FILE_PREFIX + "cgrasp_file";
CMOLMOLWORM_FILE = FILE_PREFIX + "cmolmol-worm.mac";	
CMOLMOLSURFACE_FILE = FILE_PREFIX + "cmolmol-surface.mac";	
CTEMPERATURE_FILE = FILE_PREFIX + "ctemperature.pdb";	
CPYMOL_FILE = FILE_PREFIX + "cpymol.pml";
ERROR_FILE = FILE_PREFIX + "error.log";

m = re.search(r'(\d)',OPT); 
if m: favorite_color = int(m.group(1))
else: favorite_color = 0

if re.search(r'i',OPT): invert_colors = 1
else: invert_colors = 0

if not BLAST_EXIST:
##############################################
#      JUST GET SEQUENCE FROM PDB FILE       #
##############################################
  if PDB_EXIST:
  	writeout("Here is the sequence from your PDB file:<BR><BR><B>") 
  	pdb = PDB(PDB_lines)
  	seq = pdb.seq()
  	seq.printseq('')
	writeout("</B><BR><BR>You might want to <A HREF=http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CDD_SEARCH=on&amp;FORMAT_TYPE=Text&amp;ALIGNMENT_VIEW=FlatQueryAnchoredNoIdentities target=_blank>run a BLAST search</A> on this sequence!<HR>")
  else:
  	writeout("Input data missing<BR>")

else:
  if re.search(r'b',OPT) or re.search(r'a',OPT) : 
  #########################################################################################
  #			PARSE AND ANALYZE SEQUENCE ALIGNMENT  				#
  #########################################################################################

    if re.search(r'b',OPT):
  	############################################
  	# 		PARSE BLAST FILE      	   #
  	############################################
 	# Input : BLAST_lines, OPT
  	# Output : seq = list of Sequence()
  	seq = parse_BLAST(BLAST_lines, ERROR_FILE)	

    if re.search(r'a',OPT):
  	############################################
  	# 	PARSE SEQUENCE ALIGNMENT   	   #
  	############################################
	seq = []		
	for i in range(len(BLAST_lines)):
	    seq.append(Sequence())
	    for j in range(len(BLAST_lines[i])-1): # -1 : discard carriage return
	    	seq[i].append(BLAST_lines[i][j])				

    # Get rid of redundancy if every sequence should be unique ! 
    if re.search(r'u',OPT):    
      for i in range(len(seq)-1, 0, -1):
        for j in range(i):
  	   if seq[i].equalseq(seq[j]):
	   	del seq[i]
		break

    # Write sequence alignment to file
    try: ALIGNMENT = open(ALIGNMENT_FILE, 'w')
    except IOError:
	message("Can't open file : "+ALIGNMENT_FILE+"\n")
    for i in range(len(seq)): 
  	ALIGNMENT.write(seq[i].string()+'\n')
    ALIGNMENT.close()

    # Detect Nucleic Acid Sequence
    if seq[0].is_nuc(): 
	writeout("RNA Sequence!")
	alphabet = 'ACGU-'
    else: alphabet = 'ACDEFGHIKLMNPQRSTVWXY-'
  
    # Check sequence validity
    for s in seq:
	SpellError = s.bad_spelling(alphabet)
  	if SpellError:
	   if CGI:
		writeout("Sequence Error : residue not recognized ("+SpellError+"). Check input file.<BR>")	
		writeout("Here is the <A HREF=\"/ProtSkin/tmp/"+UNIQUE_ID+"-sequence_alignment\" TARGET=\"_blank\">plain sequence alignment</A> derived from it. ");
	   else: message("Sequence Error : residue not recognized ("+SpellError+"). Check input file.<BR>")
	   sys.exit()
		
    ##########################################
    # 		ANALYZE ALIGNMENT     	   #
    ##########################################
  
    pattern = Pattern(alphabet)
    pattern.fromseq(seq)
    pattern.percent_format()
  #  print pattern.DICT
  #  pattern.printhtml(sys.stdout,favorite_color,invert_colors,alphabet)

    query = seq[0].copy()
    query.valfrompattern(pattern)
    consensus = pattern.consensus(alphabet)
    consensus.start(query.start())
  
  # query.printhtml();
  # consensus.printhtml();

  # Calculate Similarity scores :
    if re.search(r's',OPT) and (alphabet == 'ACDEFGHIKLMNPQRSTVWXY-'):  
	  no_query_pattern = Pattern(alphabet)
	  no_query_pattern.fromseq(seq[1:])
	  no_query_pattern.blosum62_weight(query,alphabet)
	  
	  query.valfrompatternsum(no_query_pattern,alphabet).divide(len(seq)-1)
  
	  full_pattern = Pattern(alphabet)
	  full_pattern.fromseq(seq)
	  full_pattern.blosum62_weight(consensus,alphabet)
	  
	  consensus.valfrompatternsum(full_pattern,alphabet).divide(len(seq))   
 		
    # Write residue conservation to file	
    try: CONSERVATION = open(CONSERVATION_FILE, 'w')
    except IOError:
	message("Can't open file : "+CONSERVATION_FILE+"\n")
    CONSERVATION.write("  Query  |Consensus\n\n")
    for j in range(query.start(), query.start()+query.length()):
  	CONSERVATION.write( "%5.2f  %s | %5.2f  %s\n" % (query.val(j), query.seq(j), consensus.val(j), consensus.seq(j)) )
    CONSERVATION.close()

    #######################################################
    # 	WRITE GRASP & MOLMOL FILES USING PDB TEMPLATE	#
    #######################################################
    if PDB_EXIST:
  	pdb = PDB(PDB_lines)

	#Check sequence
	error_res = pdb.mismatchseq(seq[0])
	if error_res:
	   if CGI:
		writeout("<B>Error at residue "+str(error_res)+" : ")
#		writeout(pdb.aa(error_res)+" versus "+seq[0].seq(error_res))
		writeout("<BR>The sequence in your PDB file does not match your BLAST query sequence.<BR></B>")
		writeout("Here is the <A HREF=\"/ProtSkin/tmp/"+UNIQUE_ID+"-sequence_alignment\" TARGET=\"_blank\">plain sequence alignment</A> derived from it. ");
	   else: message("Error at residue "+str(error_res)+" :\nThe sequence in your PDB file does not match your BLAST query sequence.\n")
	   sys.exit()		
		
	# First with query sequence :
	render = query.copy()
	
	# Write GRASP and PDB files
	pdb.color(render)
	pdb.write_grasp_list(GRASP_FILE)
	pdb.write_color_as_temperature(TEMPERATURE_FILE)
	
	# Write MOLMOL and PyMOL files
	if re.search(r's',OPT): render.reduce()
	else:			render.reduce(0,1)
	if invert_colors:	render.flin(-1,1)
	pdb.color(render)
	pdb.write_molmol_worm(MOLMOLWORM_FILE,favorite_color)
	pdb.write_molmol_surface(MOLMOLSURFACE_FILE,favorite_color)
	pdb.write_pymol(PYMOL_FILE,favorite_color)
	   
	# Again with consensus sequence :
	render = consensus.copy()
	
	# Write GRASP and PDB files
	pdb.color(render)
	pdb.write_grasp_list(CGRASP_FILE)
	pdb.write_color_as_temperature(CTEMPERATURE_FILE)
	
	# Write MOLMOL and PyMOL files
	if re.search(r's',OPT): render.reduce()
	else:			render.reduce(0,1)
	if invert_colors:	render.flin(-1,1)
	pdb.color(render)
	pdb.write_molmol_worm(CMOLMOLWORM_FILE,favorite_color)
	pdb.write_molmol_surface(CMOLMOLSURFACE_FILE,favorite_color)
	pdb.write_pymol(CPYMOL_FILE,favorite_color)
	     
	     
  else:    
  #########################################################################################
  # 	  			PARSE PLAIN LIST OF SCORES   	   			#
  #########################################################################################
    pdb = PDB(PDB_lines)

    # Read score list and color PDB object
    for l in BLAST_lines:
      if re.search(r'^#',l): continue
      m = re.search(r'(\d+)\s+(-?\d*\.?\d+)',l)
      if m: 
   	residue = int(m.group(1))
	value = float(m.group(2))
	if pdb.colorres(residue,value):
		writeout("Residue "+str(residue)+" not found<BR>")	
			
    # Get seq object from PDB object
    pdbseq = pdb.seq()	
    # Reset START of seq object, or "missing residues" will be skipped when coloring pdb object with seq (see DEV_NOTES)
    pdbseq.start(1) 
                                                                                                                                                    
    # Write GRASP file
    pdb.write_grasp_list(GRASP_FILE)
    pdb.write_color_as_temperature(TEMPERATURE_FILE)
	
    # Write MOLMOL file
    render = pdbseq.copy().reduce()
    if invert_colors: render.flin(-1,1)
    pdb.color(render)
    pdb.write_molmol_worm(MOLMOLWORM_FILE,favorite_color)
    pdb.write_molmol_surface(MOLMOLSURFACE_FILE,favorite_color)
    pdb.write_pymol(PYMOL_FILE,favorite_color)
          
	     
  #########################################################################################
  #    				 OUTPUT IF EVERYTHING WENT FINE     	 		#
  #########################################################################################

  if CGI: pass # ------ STANDALONE NO COMMENT HERE ------ #
  elif TKINTER: #Tkinter
   root = Tk()
   root.title('ProtSkin - Results')
    
   f0 = Frame(root, relief='raised', borderwidth=5)
   f0.pack(side=TOP, fill=BOTH, expand=1)

   if re.search(r'b',OPT) or re.search(r'a',OPT):   
    	Label(f0, text='The residue conservation scores for the query sequence and for the consensus sequence are displayed below :').pack()
    	
    	vbar = Scrollbar(f0)
    	vbar.pack(side=RIGHT, fill=Y)
    	hbar = Scrollbar(f0, orient=HORIZONTAL)
    	hbar.pack(side=BOTTOM, fill=X)
	canvas = Canvas(f0, yscrollcommand=vbar.set, xscrollcommand=hbar.set)
	
	TkSequence(query).printtk(canvas,1,1,favorite_color,invert_colors,'query');
 	
   	if not re.search(r's',OPT):
   		TkPattern(pattern).printtk(canvas,61,1,favorite_color,invert_colors,alphabet);
    		TkSequence(consensus).printtk(canvas,731,1,favorite_color,invert_colors,'consens');
    	else:	
		TkSequence(consensus).printtk(canvas,61,1,favorite_color,invert_colors,'consens');
	
	canvas.pack(fill=BOTH, expand=1)
   	vbar.config(command=canvas.yview)
   	hbar.config(command=canvas.xview)
	
   else:
   	if PDB_EXIST:
		Label(f0, text='The residue scores are displayed below :').pack()
		
    		vbar = Scrollbar(f0)
    		vbar.pack(side=RIGHT, fill=Y)
   		canvas = Canvas(f0, width=160, yscrollcommand=vbar.set)
 
 		TkSequence(pdbseq).printtk(canvas,1,1,favorite_color,invert_colors,'query');
 
   		canvas.pack(fill=BOTH, expand=1)
   		vbar.config(command=canvas.yview)
   
   Label(root, text='If you provided a PDB file, your macros were stored as separate files...').pack()
#   Button(root, text="OK", command=root.destroy, relief='raised', borderwidth=5).pack()
   root.mainloop()
       
   
print 'ProtSkin ' + VERSION + '(' + sys.platform + ')'

if re.search(r'^win',sys.platform): 
	raw_input('Results stored in \'Results\' directory.\nPress ENTER to exit')

	
	
     
