
#
# Routine to create a noise data structure from a dictionary
#

import linac

def noise_config(noisedict):
    nsrc = linac.Noise_Source()
    typeArr = linac.intArray_frompointer(nsrc.type)
    setsArr = linac.double_Array_frompointer(nsrc.settings)

    fielddict = ["dQ_Q","dtg","dE_ing","dsig_z","dsig_E","dchirp"]
    typedict = {"None":0,"White":1}

    for i in xrange(len(fielddict)):
        key = fielddict[i]
        try:
            entry = noisedict[key]
            print entry
            typeArr[i] = typedict[entry["Type"]]
            settings = entry["Settings"]
            for k in xrange(len(settings)):
                setsArr[linac.N_NOISE_SET*i+k] = float(settings[k])
        except:
            typeArr[i] = 0
            print "oh no"
        #print typeArr[i]
        #print setsArr[i]
    return nsrc
