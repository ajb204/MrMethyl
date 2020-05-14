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

import sys, os, re, string, cgi
from Sequence import Sequence, Sequence2

def writeout(string): sys.stdout.write(string)

def parse_BLAST(BLAST_lines, ERROR_FILE):
  ####################################################################################
  # 				PARSE "BLAST" FILE      	  		     #
  ####################################################################################
  # Input : BLAST_lines
  #	    
  # Output : seq = list of Sequence()

  identificator = []
  def checkscope(index,range_length):
  	if index >= range_length:
		writeout("Parse Error : scan out of scope. Check alignment file. ");	
		# Log error without complaint
		try: 
			ERROR = open(ERROR_FILE, 'w');
			ERROR.write("Parse Error : scan out of scope. Check alignment file below.\n");
			ERROR.write("------------------------------------------------------------\n");
			for l in BLAST_lines: ERROR.write(l);
  		except IOError: pass
		sys.exit()	
	
  # Detect special formats : CDD (1), MSF (2), CLUSTAL (3)
  ALIGNMENT_FORMAT = 0
  for line in BLAST_lines:
    # Detect CDD format based on occurence of the phrase "NCBI CDD Logo" 
    if re.search(r'NCBI\sCDD\sLogo',line): 
 	ALIGNMENT_FORMAT = 1; break
    # Detect MSF format based on Keywords: MSF Type Check .. 
    if re.search(r'MSF.*Type.*P.*Check.*\.\.',line):
    	ALIGNMENT_FORMAT = 2; break
    # Detect CLUSTAL format based on Keyword CLUSTAL
    if re.search(r'CLUSTAL',line):
    	ALIGNMENT_FORMAT = 3; break
    # Detect (BLAST) XML format based on string "xml version"
    if re.search(r'xml version',line):
    	ALIGNMENT_FORMAT = 4; break
          
  if ALIGNMENT_FORMAT == 4:   
     ##########################################
     #          XML FORMAT
     ##########################################
        print "XML NOT SUPPORTED YET !!<BR>Please try another format<BR><BR>"	
	l = 0
	# Until the end of the Hits
	while not re.search(r'</Iteration_hits>',BLAST_lines[l]):
	  # Reach start of Hit
  	  while not re.search(r'<Hit>',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))

  	  while not re.search(r'<Hit_id>',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	  print BLAST_lines[l]
  	  while not re.search(r'<Hsp_qseq>',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	  print BLAST_lines[l]
  	  while not re.search(r'<Hsp_hseq>',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	  print BLAST_lines[l]

	  # Reach end of Hit
 	  while not re.search(r'</Hit>',BLAST_lines[l]): 
		#print BLAST_lines[l]
		l += 1; checkscope(l,len(BLAST_lines))
	
	  l += 1
	  print "<HR>"
	  #print BLAST_lines[l][9:12]
	sys.exit()


  if ALIGNMENT_FORMAT == 3:   
     ##########################################
     #    	CLUSTAL FORMAT
     ##########################################
	first_res = 1
     	
	# Skip CLUSTAL line
     	l = 0
     	while not re.search(r'CLUSTAL',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	l += 1
	
	# Reach first bloc
	while not re.search(r'^(\S+)\s+([A-Za-z\-\.]+)',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	l1 = l
	
	# Get identificators
     	while re.search(r'^(\S+)\s+([A-Za-z\-\.]+)',BLAST_lines[l]):
		try: 
			identificator.append(re.search(r'^(\S+)\s+([A-Za-z\-\.]+)',BLAST_lines[l]).group(1))
			l += 1
		except: writeout("CLUSTAL Parse Error : Check your CLUSTAL file."); sys.exit()
								
	# Parse from first block again for sequence parts
     	l = l1
	seq = [] # Format is : [ {'id':'SEQ', 'id':'SEQ'} {'id'...]
     	while l<len(BLAST_lines):
     		if re.search(r'^(\S+)\s+([A-Za-z\-\.]+)',BLAST_lines[l]):
	   		# Extract sequence and store, replacing dots by dashes
	   		seq = seq + [{}]
	   		while re.search(r'^(\S+)\s+([A-Za-z\-\.]+)',BLAST_lines[l]):
#	   			print BLAST_lines[l]
	   			m = re.search(r'^(\S+)\s+([A-Za-z\-\.]+)',BLAST_lines[l])
	   			seq[-1][ m.group(1) ] = re.sub(r'\.','-',m.group(2))
				l += 1
		l += 1	
				
  elif ALIGNMENT_FORMAT == 2:   
     ##########################################
     #    	MSF FORMAT
     ##########################################
	first_res = 1
     	     
	# Reach first bloc
     	l = 0
     	while not re.search(r'Name.*Len.*Check.*Weight.*',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	
	# Get identificators
     	while re.search(r'Name.*Len.*Check.*Weight.*',BLAST_lines[l]):
		try: 
			identificator.append(re.search(r'Name\:\s*(\S*)\s*.*Len.*Check.*Weight.*',BLAST_lines[l]).group(1))
			l += 1
		except: writeout("MSF Parse Error : Check your MSF file."); sys.exit()
			
	# Reach end of header
     	while not re.search(r'\/\/',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
	l += 1			
				
	# Parse end of file for sequence parts
     	seq = [] # Format is : [ {'id':'SEQ', 'id':'SEQ'} {'id'...]
     	while l<len(BLAST_lines):
     		if re.search(r'(\S+)\s+([A-Za-z\-\.\s]+)',BLAST_lines[l]):
	   		# Extract sequence and store, replacing dots by dashes
	   		seq = seq + [{}]
	   		while re.search(r'(\S+)\s+([A-Za-z\-\.\s]+)',BLAST_lines[l]):
#	   			print BLAST_lines[l]
	   			m = re.search(r'(\S+)\s+([A-Za-z\-\.\s]+)',BLAST_lines[l])
	   			seq[-1][ m.group(1) ] = re.sub(r'\.','-',m.group(2))
				l += 1
		l += 1	
							
  elif ALIGNMENT_FORMAT == 1: 
     ##########################################
     #    CONSERVED DOMAIN DATABASE INPUT
     ##########################################
     
     # Reach ruler
     l = 0
     while not re.search(r'\.\.\.\.\*\.\.\.\.',BLAST_lines[l]): l += 1;	checkscope(l,len(BLAST_lines))
     parse_shift = BLAST_lines[l].index("....*....")   
     
     # Reach first bloc
     l = 0
     while not re.search(r'^\s?query\s+(\d+)',BLAST_lines[l]): l += 1;	checkscope(l,len(BLAST_lines))
     # Get first residue number
     m = re.search(r'query\s+(\d+)',BLAST_lines[l])
     first_res = int(m.group(1))
     
     # Get identificators
     while re.search(r'\w',BLAST_lines[l]):
     	identificator.append(BLAST_lines[l][0:parse_shift-5])
	l += 1
	
     # Parse file again for sequence parts
     seq = [] # Format is : [ {'id':'SEQ', 'id':'SEQ'} {'id'...]
     l = 0
     while l<len(BLAST_lines):
     	if re.search(r'^\s?query\s+\d+',BLAST_lines[l]):
	   # Extract sequence and store
	   seq = seq + [{}]
	   while BLAST_lines[l] != '\n':
#	      print BLAST_lines[l]
	      if re.search(r'(\D+)\s+\d+',BLAST_lines[l][parse_shift:]):
	   	m = re.search(r'(\D+)\s+\d+',BLAST_lines[l][parse_shift:])
	   	seq[-1][ BLAST_lines[l][0:parse_shift-5] ] = m.group(1)
	      l += 1
	l += 1	

  else: 
     ##########################################
     # FLAT QUERY ANCHORED WITHOUT IDENTITIES
     ##########################################
     
     # Reach first bloc
     l = 0
     while not re.search(r'^ALIGNMENT',BLAST_lines[l]) and not re.search(r'^.lignment',BLAST_lines[l]): l += 1;	checkscope(l,len(BLAST_lines))
     while not re.search(r'^(\S+)\s+(\d+)\s+\D+\s\d+',BLAST_lines[l]): l += 1; checkscope(l,len(BLAST_lines))
     # Get query identificator and first residue number
     m = re.search(r'(\w+)\s+(\d+)\s+',BLAST_lines[l])
     query_identificator = m.group(1)
     identificator.append(m.group(1))
     first_res = int(m.group(2))

     # Reach sequence description list
     l = 0
     while not re.search(r'Sequences producing significant alignments',BLAST_lines[l]):	l += 1;	checkscope(l,len(BLAST_lines))
     # Get other identificators from header (NB: this is not reliable !!)  
     #l += 2
     #while re.search(r'\w+\|+(\w+)',BLAST_lines[l]):
	#identificator.append(re.search(r'\w+\|+(\w+)',BLAST_lines[l]).group(1))
	#l += 1
 
     # Parse file again to determine parse shifts
     parse_shift_list = []
     l = 0
     while not re.search(r'^ALIGNMENT',BLAST_lines[l]) and not re.search(r'^.lignment',BLAST_lines[l]): l += 1;	checkscope(l,len(BLAST_lines))
     while l<len(BLAST_lines):
     	if re.search(r'^(\S+)\s+(\d+)\s+\D+\s\d+',BLAST_lines[l]):
	   # Extract clear sequence and store in identificator hash table
	   parse_shift_list.append(0)
	   while re.search(r'^(\S+)\s+(\d+)\s+\D+\s\d+',BLAST_lines[l]):
	   	m = re.search(r'^(\S+)\s+(\d+)\s+\D+\s\d+',BLAST_lines[l])
		if len(m.group(2)) > parse_shift_list[-1]: parse_shift_list[-1] = len(m.group(2))		
		l += 1
	l += 1
 
     # Parse file again for sequence parts
     seq = [] # Format is : [ {'id':'SEQ', 'id':'SEQ'} {'id'...]
     l = 0
     while not re.search(r'^ALIGNMENT',BLAST_lines[l]) and not re.search(r'^.lignment',BLAST_lines[l]): l += 1;	checkscope(l,len(BLAST_lines))
     while l<len(BLAST_lines) and not re.search(r'^Matrix',BLAST_lines[l]):
     	if re.search(r'^(\S+)\s+(\d+\s+\D+)\s\d+',BLAST_lines[l]):
	   # Extract clear sequence and store in identificator hash table
	   seq = seq + [{}]
	   parse_shift = parse_shift_list.pop(0) + 1
	   while re.search(r'^(\S+)\s+(\d*\s+\D+)',BLAST_lines[l]): 
		# The * here is to tolerate "empty line" of insertions only (with no res number!)
	   	m = re.search(r'^(\S+)\s+(\d*\s+\D+)',BLAST_lines[l])
		seq[-1][ m.group(1) ] = m.group(2)[parse_shift:]
		l += 1
	l += 1
				
     # Scan seq to find identificators :
     identificator = []
     for part in seq:
	identificator += part.keys()
     # Remove duplicates (NB: initial order is altered)
     identificator.sort()
     for i in range(len(identificator)-1,0,-1):
        if identificator[i] == identificator[i-1]: del identificator[i]
     # put query identificator in the first place :
     identificator.remove(query_identificator)
     identificator.insert(0,query_identificator)

  ####################################################################################
  # 				PROCESS SEQUENCES     	   			     #
  ####################################################################################
  # Fill blank sequence parts, turn to uppercase and merge

  fullseq = []  # Format is : ['SEQ' 'SEQ' 'SEQ' ... ]
  for i in range(len(identificator)):  
     fullseq.append('')
     for part in seq:  
     	if not part.has_key(identificator[i]):
		part[identificator[i]] = ' ' * len(part[identificator[0]])
	fullseq[i] += part[identificator[i]].upper()
					
  # Split strings into sequence objects starting at 1
  # Skip insertions in query sequence and fill void sequence ends
  seq = []
  for i in range(len(fullseq)):
  	seq.append(Sequence())
	
  for j in range(len(fullseq[0])):
     # Skip insertions in query sequence :
     if re.search(r'\w',fullseq[0][j:j+1]): # we have a true residue in query sequence
	for i in range(len(fullseq)):
		# Replace any void by '-' (occurs at some sequence ends)
		m = re.search(r'(\w)',fullseq[i][j:j+1])
		if m: letter = m.group(1)
		else: letter = '-'
		seq[i].append(letter)
		
  seq[0].start(first_res)
  return seq

