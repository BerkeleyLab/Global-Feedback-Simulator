
#
# Top level routine to import json configuration files 
#

import json

import linac
import readjson_accelerator
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
    from readjson import ReadDict
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
# Read in the Linac Accelerators from the dictionary
#
def ReadAccelerator(confdict):
    from readjson_accelerator import readgun,readlinac
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