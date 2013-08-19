
#
# Top level routine to import json configuration files 
#

import json
import re

import linac
import readjson
from bbf_config import BBF_Config
from noise_config import noise_config

def LoadConfig(files):
    #
    # Read the Configuration
    #
    confdict = ReadDict(files)

    #
    # Run the configuration routines
    #
    linp_pylist,linp_arr,gun = ReadAccelerator(confdict)
    Nlinac = len(linp_pylist)
    bbfstructs = BBF_Config(confdict,gun,linp_arr,Nlinac)
    bbf = bbfstructs[0]
    nsrc = noise_config(confdict['Noise'])
    print confdict

    return confdict, linp_pylist,linp_arr,gun, bbf, nsrc

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
            pass
        print "Loading ",fname,"..."
        masterdict.update(fdic)
    return masterdict


def ReadAccelerator(confdict):
    from readjson import readgun,readlinac
    import pydot
    #start a pydot graph for display
    graph = pydot.Dot(graph_type='digraph',rankdir="LR")
    
    #read in the gun parameters
    gun,Elast=readgun(confdict)
    
    #find the acclerator connectivty and read in parameters
    net=confdict['Accelerator']["connect"]

    allaccel=[] #array to store linac points temporarily
    for k in range(len(net)):

        if (k<len(net)-1):
            edge=pydot.Edge(confdict[net[k]]['name'],confdict[net[k+1]]['name'])
            graph.add_edge(edge)

        elem=confdict[net[k]]
        t=elem['type']
        if (t=="linac"):
            #print t
            linout,Elast=readlinac(confdict,net[k],Elast,graph)
            allaccel.append(linout)
        elif (t=="chicane"):
            print t
            #readchicance()
        else:
            print ("type called {0} not supported".format(t))

    #move pointers in allaccel into c array
    linp_arr=linac.Linac_Param_Array(len(allaccel))
    for l in range(len(allaccel)):
        linp_arr[l]=allaccel[l]

    #find accelerator name else use generic name
    try:
        name=confdict['Accelerator']['name']
    except:
        name="Accelerator"

    graph.write_png("{0}_graph.png".format(name))
    return allaccel,linp_arr,gun
