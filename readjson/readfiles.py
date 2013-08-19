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
