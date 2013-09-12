
#
# Top level routine to import json configuration files 
#

import json
import re

import linac
import readjson
from bbf_config import BBF_Config
#from bbf_config_like_oct import BBF_Config
from noise_config import noise_config

#
# Configure everything from a list of .cfg's and return C data structures
#
def LoadConfig(files):
    #
    # Read the Configuration
    #
    confdict = ReadDict(files)
    #print confdict
    #
    # Run the configuration routines
    #
    linp_pylist,linp_arr,gun = ReadAccelerator(confdict)
    Nlinac = len(linp_pylist)
    bbfstructs = BBF_Config(confdict,gun,linp_arr,Nlinac)
    bbf = bbfstructs[0]
    nsrc = noise_config(confdict['Noise'])
    #print confdict

    return confdict, linp_pylist,linp_arr,gun, bbf, nsrc

#
# Build a dictionary from a list of files
#
def ReadDict(files):
    if(type(files)==str):
        files = [files]
    masterdict = {}
    for fname in files:
        #
        # See if the filename is actually a dictionary 
        # passed in through the command line.
        #
        if re.match("^{.*}$",fname):
            fdic = json.loads(fname)
        else:
            f = open(fname)
            fdic = json.load(f)
            f.close()
        #
        # See if the user made any references inside of the file...
        #
        try:
            includes = fdic['#include']
            print "INCLUDING THESE FILES: ",includes
            incdict = ReadDict(includes)
            masterdict.update(incdict)
        except:
            print "No #include found... continuing"
        print "Loading ",fname,"..."
        OverlayDict(masterdict,fdic)
        #print masterdict
        
    return masterdict


#
# Recursively overlay two dictionaries
#
def OverlayDict(olddict,newdict):
    "Routine to recursively overlay two dictionaries"
    for k in newdict.iterkeys():
        if type(newdict[k])==dict:
            if olddict.has_key(k):
                OverlayDict( olddict[k], newdict[k] )
            else:
                olddict[k] = newdict[k]
        else:
            olddict[k] = newdict[k]

#
# Read in the Linac Accelerators from the dictionary
#
def ReadAccelerator(confdict):
    from readjson import readgun,readlinac
    #
    # Read in the gun parameters
    #
    gun,Elast=readgun(confdict)
    
    #
    # Find the acclerator connectivty and read in parameters
    #
    net=confdict['Accelerator']["connect"]

    #
    # Go through the net and build the list of linacs
    #
    allaccel=[]
    for k in range(len(net)):
        elem=confdict[net[k]]
        t=elem['type']
        if (t=="linac"):
            linout,Elast = readlinac(confdict,net[k],Elast)
            allaccel.append(linout)
        elif (t=="chicane"):
            print "type chicane not implemented yet"
            #print t
            #readchicance()
        else:
            print "type called {0} not supported".format(t)

    # 
    # Copy the pointers in the python list 'allaccel' into a C array 'linp_arr'
    #
    linp_arr=linac.Linac_Param_Array(len(allaccel))
    for l in range(len(allaccel)):
        linp_arr[l]=allaccel[l]

    #
    # Find accelerator name else use generic name
    #
    try:
        name=confdict['Accelerator']['name']
    except:
        name="Accelerator"
    # But we don't actually care about the name anymore...
    
    return allaccel,linp_arr,gun

