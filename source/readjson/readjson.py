   
def jsontodict(filename, defaultfile="default.cfg", Verbose=False):
    #
    #a little helper routine for reading in a json file and a optionally
    # a default file and return a python dictionary
    #
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


def readentry(dictin,entry,localdic=None, safedic={}):
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
    try:
       
        out=eval(str(entry),{},safedic)
        if isinstance(out,str) or isinstance(out,int) or isinstance(out, unicode):
            out = float(out)
        elif isinstance(out,list):
            out = [float(x) if isinstance(x,int) or isinstance(x,str) or isinstance(x,unicode) else x for x in out]
    
    #if the entry can not be evaluated look at error to get missing entry     
    except NameError as e: 
        #pull out the missing variable from the expression
        name=str(e).split("'")[1] 
        print 'Looking for {0}'.format(name)

        #search the dictionary for the entry
        #and add to local variables and retry the evaluation
        try:
            newentry=dictin[str(name)]
        except KeyError:
            try:
                newentry=localdic[str(name)]
            except KeyError:
                print '{0} has no numeric evaluation in dictionary'.format(name)
                newentry = str(name)

        safedic[name]=readentry(dictin,newentry,None,localdic) 
        out=readentry(dictin,entry,localdic,safedic)

        #except Exception:
            #go down the dictionary until a value cannot be found
            #print out what value is missing and return failed entry
            #print '{0} has no numeric evaluation in dictionary'.format(name)
            #out = entry   
   
    except TypeError as e:
        print 'Oops! TypeError: ' + str(e)
        return entry

    return out 