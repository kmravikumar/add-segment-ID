#!/usr/bin/python
#
# ./add_segid.py -f input -o output [-p pdb_format] [-n output_index] [-s segid_info_file]
#
#    ################################
#
#    Input options
#     -f  file   reqd  Input .pdb file (PDB.org/CHARMM/Gromacs pdb files)
#     -o  file   reqd  Output .pdb file
#     -p  string opt   Output PDB format "pdb" or "charmm"
#                      default is "charmm" 
#     -n  file   opt   Output index file with segids - used in Gromacs 
#     -s  file   opt   Input segid_info_file with residue information for building segment IDs
#                      eg segid_info_file:
#
#                   1 - 13
#                   14 - 30
#                   31 - 270
#               
#               The above segid_info_file means that the 
#               pdb input file has 3 segments
#               with residue numbers 1-13, 14-30, 31-270.
#
#               In case this file is not given as input the script 
#               will automatically try to assign segment IDs. (read below)
#
#    eg usage:
#       1) ./add_segid.py -f 1BKV.pdb -o 1bkv.pdb -p pdb -n index_file.ndx
#       2) ./add_segid.py -f 1BKV.pdb -o 1bkv.pdb -s resid.dat
#
#    NOTE: If atom indices are not proper (like some PDB.org files)
#          then the index file will not be consistent with .pdb file.
#          In such cases use pdb2gmx first and get an output .pdb file.
#          Then run this script.
#
#          Double check the index file after you execute the script.
#
#
#    ################################
#
# The script will convert input_file into output_file 
# which is either Protein Data Bank type pdb file (or)
# Charmm type pdb file (key difference is in the way 
# residue numbers are read when residue number > 999).
#
# Often a pdb file may not have segment ID (like Gromacs generated pdb file)
# In such cases it might we difficult to visualize the molecule
# in VMD.In such cases we can use this script to add segment IDs to any pdb
# file.
#
# Automatically assign segment IDs 
# when segid_info_file is not given :
#
# IDEA:   New segments are created whenever the distance between
#         C=O carbon and next residue N-H nitrogen is > 3.00A.
# CAVEAT:  The above idea will work only on protein/aminoacid files.
#          (when -s flag is not used)
#          


from string import *
from sys import argv
import sys
from os import rename
from math import sqrt
import getopt

try:
    import psyco
    psyco.full()
except:
    print "Install python-psyco"
    print "It will help speed up the process"


def read_string(line,start,end=0):
    """
    read line(string) from start postition to end position
    """
    
    if end==0:
        end = start
    return line[start-1:end]


class pdb_line:
    """
    Read different fields of a single "ATOM" line
    in a pdb file.

        self.ATOM    = read_string(line,1,6)
        self.serial  = read_string(line,7,11)
        self.name    = read_string(line,13,16)
        self.altloc  = read_string(line,17)
        self.resname = read_string(line,18,21)
        self.chain   = read_string(line,22)

        self.resnum        = read_string(line,23,26)
        self.charmm_resnum = read_string(line,24,27)
        self.icode         = read_string(line,27)

        self.x = read_string(line,31,38) - X,Y,Z positions
        self.y = read_string(line,39,46)
        self.z = read_string(line,47,54)

        self.occ     = read_string(line,55,60)
        self.tempfac = read_string(line,61,66)
        self.segid   = read_string(line,73,76)
        self.elename = read_string(line,77,78)
        self.charge  = read_string(line,79,80)

        self.data = self.write_pdb()        

    
    """

    def init(self,line=""):

        # If the length of the line is 75(charmm pdb) 
        #     make it 80 which is standard pdb line length
        # If length of line != 80 then make the line empty.
        if len(line) > 1 and len(line) < 81 and read_string(line,1,4) == 'ATOM':
            line = line[:] + " "*(80-len(line)) + "\n" 

        if len(line) != 81 :
            line = " "*80+"\n"

        if read_string(line,1,6) != 'HETATM' \
           and read_string(line,1,4) != 'ATOM':
            line = " "*80+"\n"

        self.ATOM    = read_string(line,1,6)
        self.serial  = read_string(line,7,11)
        self.name    = read_string(line,13,16)
        self.altloc  = read_string(line,17)
        self.resname = read_string(line,18,21)
        self.chain   = read_string(line,22)

        self.resnum        = read_string(line,23,26)
        self.charmm_resnum = read_string(line,24,27)
        self.icode         = read_string(line,27)

        # Sometimes pdb file has residue numbers like 74A
        # In such cases use only the interger part for residue number
        if read_string(line,27).isalpha():
            self.charmm_resnum = read_string(line,24,26)+" "
            self.icode = " "

        try:
            if self.resnum != "    " and (int(self.resnum) < int(self.charmm_resnum)):
                self.resnum = self.charmm_resnum
                self.icode = " "
        except:
            if self.resnum != "    ":
                self.resnum = self.charmm_resnum
                self.icode = " "
        

        self.x = read_string(line,31,38)
        self.y = read_string(line,39,46)
        self.z = read_string(line,47,54)

        self.occ     = read_string(line,55,60)
        self.tempfac = read_string(line,61,66)
        self.segid   = read_string(line,73,76)
        self.elename = read_string(line,77,78)
        self.charge  = read_string(line,79,80)

        self.data = self.write_charmm()        

    def write_charmm(self):
        """
        Write the line in charmm pdb format
        """

        line = self.ATOM+self.serial+" "+self.name+self.altloc+self.resname
        line = line+""+self.chain+" "+self.resnum+"   "
        line = line+self.x+self.y+self.z+self.occ+self.tempfac+"      "
        line = line+self.segid+self.elename+self.charge+"\n" 
        return line            

    def write_pdb(self):    
        """
        Write the line in PDB.org pdb format
        """

        line = self.ATOM+self.serial+" "+self.name+self.altloc+self.resname
        line = line+""+self.chain+self.resnum+self.icode+"   "
        line = line+self.x+self.y+self.z+self.occ+self.tempfac+"      "
        line = line+self.segid+self.elename+self.charge+"\n" 
        return line



def check_line(line):
    """
    To check if the line is a valid line for pdb file in charmm
    valid line returns 1 while invalid line returns 0
    """

    atom = pdb_line(line)

    if atom.ATOM=='ATOM  ' or atom.ATOM=='HETATM':
        return 1
    return 0



def read_relevant_line(file):
    """
    To read the next relevant line from the raw pdb file.
    Omit all unwanted lines and read only 
    "ATOM" or "HETATM" entries
    """

    line1 = file.readline() 
    
    # loop over next lines and check if it is a relevant line
    while (not check_line(line1)) and line1:
        line1 = file.readline()
        
    return line1



def make_compatible(input_filep,out_filep,format="charmm"):
    """
    Convert a pdb file into required format("charmm" or "pdb")
    """
 
    print "  Make file %s into a %s type file %s" \
           % (input_filep.name,format,out_filep.name)

    input_filep.seek(0,0)

    print "  %s --> %s " \
        % (input_filep.name,out_filep.name)
    
 
    # initial segid = 1 & resid = 1
    segid = 1
    resnum = 1
    line1 = read_relevant_line(input_filep)
    atom1 = pdb_line(line1) 
 
    # If the first resnum is not 1
    # then keep the residue number - do not initialize it as 1
    if atom1.resnum and int(atom1.resnum)!=1:
        resnum = int(atom1.resnum)
 
    atom1.resnum = rjust(str(resnum),4)
    original_segid = atom1.segid
    original_resnum = atom1.resnum
    atom1.segid = 'E'+ zfill(str(segid),3)
    chain = atom1.chain
    atom1.chain = " "
    atom1.ATOM = 'ATOM  '
    atom1.elename = '  '

    if format == "charmm":
        out_filep.write(atom1.write_charmm())
    else:
        out_filep.write(atom1.write_pdb()) 

    line2 = read_relevant_line(input_filep)
    atom2 = pdb_line(line2)
    first_water = 0
 
    while line2:

       # If resname is not water 
       if atom2.resname not in ["HOH ","TIP3","TIP4",\
                                "TIP5","SPC ","SPCE","WAT ","SOL "]:

           if (original_resnum != atom2.resnum or \
	           atom1.resname != atom2.resname or \
	           atom2.name == ' N  ') :
	
               resnum = resnum + 1
	
           original_resnum = atom2.resnum 
           atom2.resnum = rjust(str(resnum),4)
	
           # If there is no segid then make 
           # segid of next atom same as previous atom
           if atom2.segid == "    ":
               atom2.segid = atom1.segid
	       
           # If segid of atom2 and original segid of atom1 were same
           # then make segid of atom2 same as current segid of atom1
           if atom2.segid == original_segid:
               atom2.segid = atom1.segid
	
	
           # Compare distance between coordinates of previous 
           # O atom and next N atom.
           # If the distance is > 3.00 A then increment segid
           if atom2.name == ' O  ':
               last_O_coor = float(atom2.x),float(atom2.y),float(atom2.z)
           if atom2.name == ' N  ':
               dist  = (last_O_coor[0]-float(atom2.x))**2
               dist += (last_O_coor[1]-float(atom2.y))**2
               dist += (last_O_coor[2]-float(atom2.z))**2
               dist = sqrt(dist)
               if dist > 3.00:
                   # Just make sure atom1.segid != atom2.segid
                   # rest will be taken care in the next "if" condition
                   atom2.segid = " "
    

           # If chainID or segid of atom1 and atom2 differ
           # then increment segid
           if chain != atom2.chain  or atom2.segid!=atom1.segid:
               original_segid = atom1.segid
               resnum = 1
               if int(atom1.resnum)!=1:
                   resnum = int(original_resnum)
    
               atom2.resnum = rjust(str(resnum),4)
               segid = segid + 1
               atom2.segid = 'E'+ zfill(str(segid),3)
    
               # If resname is ACY then make
               # segid as ACY
               if atom2.resname == "ACY ":
                   atom2.segid = " ACY"


       # If resname is HOH then segid="WAT"
       # For every OH increment resnum by 1
       if atom2.resname in ["HOH ","TIP3","TIP4","TIP5",\
                            "SPC ","SPCE","WAT ","SOL "]:
           atom2.segid = " WAT"

           # for first water make resnum as 1
           if first_water ==0:
               resnum = 0
               first_water = 1

           if atom2.name.strip()[0] == "O":
               resnum = resnum + 1 

           atom2.resnum = rjust(str(resnum),4)

       chain = atom2.chain
       atom2.chain=" "
       atom2.ATOM = 'ATOM  '
       atom2.elename = '  '
       if format == "charmm":
           out_filep.write(atom2.write_charmm())
       else:
           out_filep.write(atom2.write_pdb())           

       atom1 = atom2
       line2 = read_relevant_line(input_filep)
       atom2 = pdb_line(line2)

    return 1



def check_file(file_name):
    """ 
    To check if a file exists or is empty 
    return 0 if file is empty
    return 1 if file is present and non-empty
    """

    try:
        filep = open(file_name,'r')
        if filep.readline():
            filep.close()
            return 1
        else:
            print filep.name, "file is empty"
            return 0
    except IOError:
        print "File does not exist"
        return 0


def make_index(input_filep,output_filep): 
    """
    Make an index file with segids
    """

    line = read_relevant_line(input_filep)
    atom = pdb_line(line)
    segid = ""
    atom_index = 0
    entry_counter = 0

    while line:
        atom_index += 1
        entry_counter += 1
        
        # if segid != previous_segid - create new
        # index with segid as its name
        if segid != atom.segid:
            data = "\n\n"+ "[ "+ atom.segid  +" ]" + "\n"
            output_filep.write(data)
            segid = atom.segid
            entry_counter = 1
            
        # write atom_index into the file
        data = "%8d " %(atom_index)
        output_filep.write(data)
        
        # put 8 indices per line
        if entry_counter == 8:
            output_filep.write("\n")
            entry_counter = 0
            
        line = read_relevant_line(input_filep)    
        atom = pdb_line(line)            

    return 1
    

def usage():
    """ 
    Usage help for the program
    """

    print """ 
    ################################

    Input options
     -f  reqd  Input .pdb file
     -o  reqd  Output .pdb file
     -p  opt   Output PDB format "pdb" or "charmm"
               default is "charmm" 
     -n  opt   Output index file with segids (default: segid.ndx)
     -s  opt   Input residue information for building segment IDs
                   1 - 13
                   14 - 30
                   31 - 270
               
               The above segid_info_file says that the 
               pdb input file has 3 segments
               with residue numbers 1-13, 14-30, 31-270.

    eg usage:
       1) ./add_segid.py -f 1BKV.pdb -o 1bkv.pdb -p pdb -n index_file.ndx
       2) ./add_segid.py -f 1BKV.pdb -o 1bkv.pdb -s resids.dat

    NOTE: If atom indices are not proper (like some PDB.org files)
          then the index file will not be consistent with .pdb file.
          In such cases use pdb2gmx first and get an output .pdb file.
          Then run this script.

          Double check the index file after you execute the script.


    ################################
    """
    return 1


def print_input(opts,args):
    """
    Print the input arguments given by user
    """

    print "################################"
    print "  Input "
    for o,a in opts:
        print "  ",o,"  ",a
    print "################################"

    return 1

def make_segid(input_filep,out_filep,format="charmm",resid_filep=" "):
    """
    Make segment names using residue sequence information
    given in resid_filep
    """
    
    # stores maximum number of segments in -s input
    max_segid = 0

    # begin[] and end[] store start and end
    # residue numbers from the -s option
    begin = []
    end = []

    # run through the -s input file and store
    # required start and end residue ids in begin and end arrays
    # Also do checks to make sure the input file format 
    # in -s file is correct
    for res_seq in resid_filep:
        if len(res_seq.split("-")) != 2:
            print "Bad file input -s ",resid_filep.name
            usage()
            sys.exit()
        else:
            try:
                begin.append(int(res_seq.split("-")[0]))
                end.append(int(res_seq.split("-")[1]))
                max_segid += 1
            except:
                print "Bad residue input in file -s",resid_filep.name,". Exiting"
                usage()
                sys.exit()

    # For each atom in the input pdb file
    # identify which begin-end sequence it belongs to using resnum.
    # Then assign corresponding segmentID
    # If an atom doesn't belong to any of the residue sequences
    # then assign a unique segid ("E" + zfill(str(max_segid+1),3))
    line = read_relevant_line(input_filep)
    while line:
        atom = pdb_line(line)
        atom.ATOM = "ATOM  "
        for i in range(0,max_segid):
            if int(atom.resnum) in range(begin[i],end[i]+1):
                atom.segid = "E" + zfill(str(i+1),3)
                break
            else:
                atom.segid = "E" + zfill(str(max_segid+1),3)
        if format == "charmm":
            out_filep.write(atom.write_charmm())
        else:
            out_filep.write(atom.write_pdb())           
        line = read_relevant_line(input_filep)

    return 1


def main(argv):    
    """
    Main function for input 
    arguments parsing
    """
    
    try:                                
        opts, args = getopt.getopt(argv, "hf:o:p:n:s:") 

    except getopt.GetoptError,err:           
        print "\n  ",str(err)
        usage()                          
        sys.exit(2) 

    format = "charmm"
    print_input(opts,args)

    if len(opts) == 0:
        print "No arguments given \n"
        print "    Printig Help"
        usage()
        sys.exit(2)
        
    index_option = 0
    resid_option = 0

    for o, a in opts:
        if o in ("-f"):
            inp_file_name = a
            if not check_file(inp_file_name):
                print "Check",o,"filename. Exiting"
                sys.exit()
            
        elif o in ("-o"):
            out_file_name = a
            
        elif o in ("-p"):
            if a in ["pdb","charmm"]:
                format = a
            else:
                print "Unrecognized -p format option"
                sys.exit()

        elif o in ("-h"):
            print " -h Help option"
            usage()
            sys.exit()
        
        elif o in ("-n"):
            index_option = 1
            index_filep = open(a,"w")

        elif o in ("-s"):
            try:
                resid_filep = open(a)
                resid_option = 1
            except:
                print "File -s ",a, "does not exist, Exiting"
                sys.exit()
        else:
            assert False, "unhandled option"

    input_filep = open(inp_file_name,'r')
    output_filep = open(out_file_name, 'w')

    if resid_option == 0:
        make_compatible(input_filep,output_filep,format)
    else:
        make_segid(input_filep,output_filep,format,resid_filep)

    input_filep.close()
    output_filep.close()

    if index_option == 1:
        output_filep = open(out_file_name, 'r')
        make_index(output_filep,index_filep)
        output_filep.close()
        index_filep.close()

    return 1


if name == "main":
    main(argv[1:])



