#!/usr/bin/python

# Python parser of Accelerator Configuration parameters from JSON file.

 # Get list of JSON files and return dictionary
def loadConfig(files, Verbose=True):

    from readjson import ReadDict
    # Read the configuration file(s)
    confdict = ReadDict(files, Verbose)

    return confdict

def ParseAccelerator(file_list, Verbose=True):

    # Import accelerator-specific functions
    import readjson_accelerator as read_acc
    # Extract dictionary from JSON configuration file(s)
    conf_dictionary = loadConfig(file_list, Verbose)

    # Instantiate simulation and accelerator classes using configuration in the dictionary
    simulation, accelerator = read_acc.readConfiguration(conf_dictionary)

    return simulation, accelerator

def DoSimulation():

    file_list =  ["configfiles/LCLS-II/default_accelerator.json", "configfiles/LCLS-II/LCLS-II_accelerator.json"]

    simulation, accelerator = ParseAccelerator(file_list)

    print simulation
    print accelerator

# If this file is called as main, run it from command line arguments,
# or some default settings.
#
if __name__=="__main__":
    DoSimulation()
