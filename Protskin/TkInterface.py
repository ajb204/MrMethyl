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
import re
from Tkinter import *
from FileDialog import LoadFileDialog

class FilenameEntry(Frame):
    def __init__(self, master, width):
	Frame.__init__(self, master)
	self.filename = StringVar()
	Entry(self, textvariable=self.filename, width=width).pack(side=LEFT, fill=X, expand=1)
	Button(self, text="Browse...", command=self.browse).pack(side=RIGHT)

    def browse(self):
	file = LoadFileDialog(self).go(pattern='*')
	if file:
	    self.filename.set(file)

    def get(self):
	return self.filename.get()
			
    def set(self, string=''):
    	self.filename.set(string)
	
class Interface:
    def __init__(self):    	
	self.root = Tk()
	self.root.title('ProtSkin')
#	print root.keys()
	f0 = Frame(self.root, relief=RAISED, borderwidth=5)
	f0.pack(side=TOP, fill=X, expand=1)

	#---------
	f1 = Frame(f0)
	f1.pack(side=LEFT, fill=Y)
	
	Label(f1, text='Residue score from :', height=2).pack(side=TOP)
	#---------
	
	#---------
	f2 = Frame(f0)
	f2.pack(anchor=N, fill=X, expand=1)
	
	self.filename1 = FilenameEntry(f2, 50)
	self.filename1.pack(side=TOP, fill=BOTH)
	
	f20 = Frame(f2)
	f20.pack(anchor=W, fill=X, expand=1)

	self.OPT1 = StringVar()	
	f21 = Frame(f20, relief=GROOVE, borderwidth=2)
	f21.pack(anchor=W)
    	Radiobutton(f21, text="Plain list of scores", variable=self.OPT1, value='').pack(anchor=W)
	
	f22 = Frame(f20, relief=GROOVE, borderwidth=2)
	f22.pack(anchor=W, fill=X)
	########---------
	f221 = Frame(f22)
	f221.pack(side=LEFT)    
		
 	Radiobutton(f221, text="BLAST, CLUSTAL or MSF file", variable=self.OPT1, value='-b').pack(anchor=W)
    	Radiobutton(f221, text="Sequence alignment (from a previous query)", variable=self.OPT1, value='-a').pack(anchor=W)		

	
	f222 = Frame(f22, relief=GROOVE, borderwidth=2)
	f222.pack(side=RIGHT)    
	
	self.OPT2 = StringVar()
    	Checkbutton(f222, text="Discard redundant sequences", variable=self.OPT2, offvalue='', onvalue='-u').pack(anchor=W)
	
	Label(f222, text='Residue conservation score based on :').pack(side=LEFT)
	self.OPT3 = StringVar()
    	Radiobutton(f222, text="Identity", variable=self.OPT3, value='').pack(side=LEFT)
    	Radiobutton(f222, text="Similarity", variable=self.OPT3, value='-s').pack(side=LEFT)
	
	########---------
	f23 = Frame(f20, relief=GROOVE, borderwidth=2)
	f23.pack(anchor=W)
	#---------

	Label(self.root).pack(side=TOP)

	#---------
	f3 = Frame(self.root, relief=RAISED, borderwidth=5)
	f3.pack(anchor=W, fill=X, expand=1)
	
	f31 = Frame(f3)
	f31.pack(anchor=W, fill=X, expand=1)
	Label(f31, text='PDB file :', height=2).pack(side=LEFT)
	self.filename2 = FilenameEntry(f31, 50)
	self.filename2.pack(side=TOP, fill=BOTH)
	
	f32 = Frame(f3)
	f32.pack(anchor=W, fill=X, expand=1)
	Label(f32, text='Favorite Color :').pack(side=LEFT)
	self.OPT5 = StringVar()
	f321 = Frame(f32, relief=GROOVE, borderwidth=5, background='#FF0000')
	f321.pack(side=LEFT)
    	Radiobutton(f321, text="red", variable=self.OPT5, value='-4').pack(side=LEFT)
	f322 = Frame(f32, relief=GROOVE, borderwidth=5, background='#FF00FF')
	f322.pack(side=LEFT)
    	Radiobutton(f322, text="pink", variable=self.OPT5, value='-5').pack(side=LEFT)
	f323 = Frame(f32, relief=GROOVE, borderwidth=5, background='#0000FF')
	f323.pack(side=LEFT)
    	Radiobutton(f323, text="blue", variable=self.OPT5, value='-1').pack(side=LEFT)
	f324 = Frame(f32, relief=GROOVE, borderwidth=5, background='#00FF00')
	f324.pack(side=LEFT)
    	Radiobutton(f324, text="green", variable=self.OPT5, value='-2').pack(side=LEFT)
	f325 = Frame(f32, relief=GROOVE, borderwidth=5, background='#FFFF00')
	f325.pack(side=LEFT)
    	Radiobutton(f325, text="yellow", variable=self.OPT5, value='-6').pack(side=LEFT)
	f326 = Frame(f32, relief=GROOVE, borderwidth=5, background='#FF8800')
	f326.pack(side=LEFT)
    	Radiobutton(f326, text="orange", variable=self.OPT5, value='-7').pack(side=LEFT)
	f327 = Frame(f32, relief=GROOVE, borderwidth=5, background='#333333')
	f327.pack(side=LEFT)
    	Radiobutton(f327, text="black", variable=self.OPT5, value='-0').pack(side=LEFT)
	self.OPT6 = StringVar()
    	Checkbutton(f32, text="Invert colors", height=2, variable=self.OPT6, offvalue='', onvalue='-i').pack(anchor=W)	
	#---------
	
	Label(self.root).pack(side=TOP)

	#---------
    	Button(self.root, text="Default", command=self.init).pack(side=LEFT)
   	Button(self.root, text="Cancel", command=self.cancel).pack(side=RIGHT)
     	Button(self.root, text="Process", command=self.ok).pack(side=RIGHT)
	#---------

	self.init()
	self.root.mainloop()
	
    def init(self):
    	self.filename1.set('')
	self.filename2.set('')
	self.OPT1.set('-b')
	self.OPT2.set('-u')
	self.OPT3.set('')
	self.OPT5.set('-4')
   	self.OPT6.set('')
	
	self.OK = 0
	
# FOR TESTING ONLY # --------
#	self.filename1.set('blast.txt')
#	self.filename2.set('blast.pdb')
#################### --------

    def ok(self):
    	self.OK = 1
	self.root.destroy()

    def cancel(self):
    	self.OK = 0
	self.root.destroy()

    def is_ok(self):
    	return self.OK

    def getarg(self):
	OPTS = ''
	OPTS += self.OPT1.get() + ' '
	OPTS += self.OPT2.get() + ' '
	OPTS += self.OPT3.get() + ' '
	OPTS += self.OPT5.get() + ' '
	OPTS += self.OPT6.get()
	OPT = re.split(r'\s+',OPTS)	
	if OPT[0] == '': OPT.pop(0)
	if OPT[-1] == '': OPT.pop(-1)
	
	OPT.append(self.filename1.get())
	OPT.append(self.filename2.get())
	
	return OPT

    def results(self):
	top = Toplevel(self.root)
	top.title("About this application...")
	msg = Message(top, text='HELLOHELLOHELLOHELLOHELLOHELLO')
	msg.pack()
	button = Button(top, text="Dismiss", command=top.quit)
	button.pack()
	mainloop()
	self.root.destroy()


   
