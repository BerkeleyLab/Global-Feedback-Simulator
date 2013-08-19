import linac

def loadall(acceleratorfile,bbffile,noisefile="noise_test.cfg",default_acc="default.cfg",default_bbf="bbf_default.cfg",Verbose=False):
    
    bbf_conf=jsontodict(bbffile,defaultfile=default_bbf, Verbose=Verbose)
    noise_conf=jsontodict(noisefile)
    a,allaccel,linp_arr,gun=loadaccelerator(acceleratorfile,defaultfile=default_acc, Verbose=Verbose)
    
    return a,allaccel,linp_arr,gun,bbf_conf,noise_conf
    
def jsontodict(filename, defaultfile="default.cfg", Verbose=False):
    import json

    #read in file of interest
    f=open(filename)
    nondefault=json.load(f)
    f.close()

    #look for a default file
    try:
        f=open(defaultfile)
        a=json.load(f)
        f.close()
    except Exception,e:
        print "No default parameters read"
        print str(e)
        a={}
    
    #add the new to the defualt overwriteing changed defualt values
    a.update(nondefault)

    return a



def loadaccelerator(filename, defaultfile="default.cfg", Verbose=False):
    import json
    import pydot

    #read in file of interest
    f=open(filename)
    nondefault=json.load(f)
    f.close()

    #look for a default file
    try:
        f=open(defaultfile)
        a=json.load(f)
        f.close()
    except Exception,e:
        print "No default parameters read"
        print str(e)
        a={}
    
    #add the new to the defualt overwriteing changed defualt values
    a.update(nondefault)
    #print a
    #start a pydot graph for display
    graph = pydot.Dot(graph_type='digraph',rankdir="LR")
    


    #read in the gun parameters
    gun,Elast=readgun(a)
    
    #find the acclerator connectivty and read in parameters

    net=a['Accelerator']["connect"]
    allaccel=[] #array to store linac points temporarily
    for k in range(len(net)):

        if (k<len(net)-1):
            edge=pydot.Edge(a[net[k]]['name'],a[net[k+1]]['name'])
            graph.add_edge(edge)

        elem=a[net[k]]
        t=elem['type']
        if (t=="linac"):
            #print t
            linout,Elast=readlinac(a,net[k],Elast,graph)
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
        name=a['Accelerator']['name']
    except:
        name="Accelerator"

    graph.write_png("{0}_graph.png".format(name))
    return a,allaccel,linp_arr,gun


def readgun(a,Verbose=False):
    #get the name of the gun object, return defualt if not found
    GUN=a["Accelerator"].get("Gun","d_Gun")
    if(Verbose):
        print "reading linac called {0}".format(a[GUN]["name"])
    
    #update the gun object with accelerator specific values
    a[GUN]=dict(a['d_Gun'].items()+a[GUN].items())

    #read gun output properties
    E=readentry(a,a[GUN]["E"])*1e9       #[GeV]->[Ev]
    sz0=readentry(a,a[GUN]["sz0"])*1e-3  #[mm]->[m]
    sd0=readentry(a,a[GUN]["sd0"])*1e-2  #[%]->fraction 
    Q=readentry(a,a[GUN]["Q"])

    #create and load the c array
    gun=linac.Gun_Param()
    linac.Gun_Config(gun,
                     E,sz0,sd0,Q)
    return gun, E

def readlinac(a,linac_key,Elast,graph,Verbose=False):
    from numpy import pi
    if(Verbose):
        print "reading linac called {0}".format(a[linac_key]['name'])
    
    #read in default blocks to use in linac and overwrite if 
    #run specific object given
    a[linac_key]=dict(a['d_linac'].items()+a[linac_key].items()) 
    
    #pull out object names (step is for convienence)
    TRF1=a[linac_key]["TRF1"]
    CLIP=a[linac_key]["CLIP"]
    TRF2=a[linac_key]["TRF2"]
    CAV=a[linac_key]["CAV"]
    RXF=a[linac_key]["RXF"]
    CHIC=a[linac_key]["Chicane"]
    CONT=a[linac_key]["Controller"]
        
    #read in defaults and replace with run specific values
    a[TRF1]=dict(a['d_TRF1'].items()+a[TRF1].items())
    a[CLIP]=dict(a['d_CLIP'].items()+a[CLIP].items())  
    a[TRF2]=dict(a['d_TRF2'].items()+a[TRF2].items())
    a[CAV]=dict(a['d_CAV'].items()+a[CAV].items())
    a[CHIC]=dict(a['d_Chicane'].items()+a[CHIC].items()) 
    a[RXF]=dict(a['d_RXF'].items()+a[RXF].items())
    a[CONT]=dict(a['d_Controller'].items()+a[CONT].items())
    
    #read in the filter poles for TRF1,TRF2,RXF
    scale=1e6  #scale pole values to [Hz] from [MHz]
    p_TRF1 = linac.complexdouble_Array(2)
    p_TRF1[0] = readentry(a,a[TRF1]['poles'][0][0],localdic=a[TRF1])*scale
    p_TRF1[1] = readentry(a,a[TRF1]['poles'][1][0],localdic=a[TRF1])*scale
    p_TRF2 = linac.complexdouble_Array(1)
    p_TRF2[0] = readentry(a,a[TRF2]['poles'][0][0],localdic=a[TRF2])*scale
    p_RXF = linac.complexdouble_Array(3)
    p_RXF[0] = readentry(a,a[RXF]['poles'][0][0],localdic=a[RXF])*scale
    p_RXF[1] = readentry(a,a[RXF]['poles'][1][0],localdic=a[RXF])*scale
    p_RXF[2] = readentry(a,a[RXF]['poles'][2][0],localdic=a[RXF])*scale

    #CLIP parameters
    saturate_c=readentry(a,a[CLIP]['sat_c'],localdic=a[CLIP])
    kly_max_v=readentry(a,a[CLIP]['kly_max_v'],localdic=a[CLIP])

    #Cavity parameters
    E=readentry(a,a[CAV]['E'],localdic=a[CAV])*1e9 #[Gev]->[ev]
    dE=(E-Elast)
    phi=readentry(a,a[CAV]['phi'],localdic=a[CAV])*pi/180  #[deg]->[rad]
    lam=readentry(a,a[CAV]['lam'],localdic=a[CAV])         #[m]
    s0=readentry(a,a[CAV]['s0'],localdic=a[CAV])*1e-3      #[mm]->[m]
    aper=readentry(a,a[CAV]['a'],localdic=a[CAV])*1e-3  #[mm]->[m]  
    # named aper because a already used
    L=readentry(a,a[CAV]['L'],localdic=a[CAV])          #[m]
    nom_grad=readentry(a,a[CAV]['nomgrad'],localdic=a[CAV])
    psd_llrf=readentry(a,a[CAV]['psd_llrf'],localdic=a[CAV])
    w0=readentry(a,a[CAV]['w0'],localdic=a[CAV])     #[rad/s]
    bunch_rep=readentry(a,a[CAV]['bunch_rep'],localdic=a[CAV])
    Q_L=readentry(a,a[CAV]['Q_L'],localdic=a[CAV])
    R_Q=readentry(a,a[CAV]['R_Q'],localdic=a[CAV])
    beta_in=readentry(a,a[CAV]['beta_in'],localdic=a[CAV])
    beta_out=readentry(a,a[CAV]['beta_out'],localdic=a[CAV])
    beta_beam=readentry(a,a[CAV]['beta_beam'],localdic=a[CAV])
    
    #chicane parameters
    R56=readentry(a,a[CHIC]['R56'],localdic=a[CHIC])
    T566=readentry(a,a[CHIC]['T566'],localdic=a[CHIC])
    
    #linac properties
    n_cav=readentry(a,a[linac_key]['n_cav'],localdic=a[linac_key])

    #controller parameters
    stable_gbw=readentry(a,a[CONT]['stable_gbw'],localdic=a[CONT])
        
    #simparamters
    dt=readentry(a,a["Simulation"]['dt'],localdic=a["Simulation"])

    lin = linac.Linac_Param()
    linac.Linac_Config(lin,
                   dt,
                   dE,R56,T566,phi,
                   lam,s0,aper,L,
                   saturate_c,
                   p_TRF1,p_TRF2,p_RXF,
                   n_cav, nom_grad,
                   psd_llrf,w0,bunch_rep,
                   Q_L,R_Q,
                   beta_in,beta_out,beta_beam,
                   kly_max_v,stable_gbw
                   )
    
    return lin,E


def readentry(dictin,entry,localdic=None):
#Daniel Driver
#Alejandro Quiruga
#LBL 2013

#function the recursive read and evaluate the entries in a dictionaty

#inputs-
    #dictin : dicitionary to be searched if entry cannot be
    #         evaluated
    #entry : the dicitionary value you would like to get
    #localdic: optional input to pass perviously found
    #          dictionary entrys down the recursive chain
 
#outputs-
    #value of entry is return as float if is can be evaluated using the
    #       in the dicitonary and otherwise it returns the string


#
# TODO: Insert license
#



    #replace localdic with an empty dictionary in nothing passed
    if localdic is None:
        localdic={}

    #try to evluate the entry to interest
    #print localdic
    try:
        out=eval(str(entry),{},localdic) 
        #if the entry can not be evaluated look at error to get missing entry     
    except NameError as e: 
        #pull out the missing variable from the expressio
        name=str(e).split("'")[1] 
        
        try:
            #search search the dicitionary for the entry
            #and add to local variables and retry the evaluation
            localdic[name]=readentry(dictin,dictin[str(name)],localdic=localdic) 
            out=readentry(dictin,entry,localdic=localdic)

        except:
            #go down the dictionary until a value cannot be found
            #print out what value is missing and return failed entry
            print '{0} has no no numeric evaluation in dictionary'.format(name)
            out = entry   
   
    except TypeError:
        return entry

    return out 



def readfilter(a,f,graph):
    
    #
    #
    #currently unused-probably does not work
    #
    #
    
    print "in readfilter"
    #
    # Create a Filter
    #
    #fil = linac.Filter_Allocate_New(3,3) #this uses malloc, lets avoid this 
    fil = linac.Filter() #declare one in python and then set it up
    linac.Filter_Allocate_In(fil,3,3)
    
    #push in some random poles

    poles = linac.complexdouble_Array(3)
    poles[0] = 1.0
    poles[1] = 1.0-2.0j
    poles[2] = 2.5j
    linac.Filter_Append_Modes(fil,poles,3)
    
    poles = linac.complexdouble_Array(1)
    poles[0] = -1.0-1.5j
    linac.Filter_Append_Modes(fil,poles,1)
    poles[0] = -1.0-2.5j
    linac.Filter_Append_Modes(fil,poles,1)

    poles = linac.complexdouble_Array(2)
    poles[0] = -1.0-3.5j
    poles[1] = 2.1
    linac.Filter_Append_Modes(fil,poles,2)
    
     # Use the swig carrays.i to get handles to the pointer arrays
    
    fil.A_modes = linac.intArray_frompointer(fil.modes)
    fil.A_coeff_start = linac.intArray_frompointer(fil.coeff_start)
    fil.A_coeffs = linac.complexdouble_Array_frompointer(fil.coeffs)
    fil.A_poles = linac.complexdouble_Array_frompointer(fil.poles)


    return 0
