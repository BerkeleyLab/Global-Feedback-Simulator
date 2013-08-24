
#
# Routine to create a noise data structure from a config dictionary
#
# Alejandro F Queiruga
# Daniel Scott Driver
# LBNL Summer 2013
#

import linac

def noise_config(noisedict):
    #
    # Allocate data structures and create handles to the underlying
    # C-arrays
    #
    nsrc = linac.Noise_Source()
    typeArr = linac.intArray_frompointer(nsrc.type)
    setsArr = linac.double_Array_frompointer(nsrc.settings)

    #
    # Dictionaries for parameters and indices
    #
    fielddict = ["dQ_Q","dtg","dE_ing","dsig_z","dsig_E","dchirp","adc_real","adc_imag"]
    typedict = {"None":0,"White":1,"Sine":2,"Chirp":3,"Step":4}

    #
    # Loop over all sources of noise we have defined, and 
    # see if the user specified them
    #
    for i in xrange(len(fielddict)):
        key = fielddict[i]
        # Make sure the key corresponds to a field that the user specified
        try:
            entry = noisedict[key]
            # Figure out its type
            typeArr[i] = typedict[entry["Type"]]
            # Write the settings for the noise into the C structure
            settings = entry["Settings"]
            for k in xrange(len(settings)):
                setsArr[linac.N_NOISE_SET*i+k] = float(settings[k])
        except:
            # Default to no noise and whine to the user
            typeArr[i] = 0
            print "Noise source ", key, " not found, defaulted to None"

    return nsrc
