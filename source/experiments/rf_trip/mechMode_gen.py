

def data_to_JSON():

    import json, copy
    data_list = [
        [59.9,   110.0,   -0.1,    0.0],
        [60.8,   90.0,    -0.04,   0.0],
        [68.3,   90.0,    -0.02,   0.0],
        [80.8,   90.0,    -0.03,   0.0],
        [93.9,   50.0,    -0.4,    0.0],
        [97.4,   300.0,   -0.2,    0.0],
        [97.7,   300.0,   -0.2,    0.0],
        [114.4,  50.0,    -0.1,    0.0],
        [117.0,  8.0,     -0.2,    0.0],
        [125.0,  15.0,    -0.6,    0.0],
        [130.2,  60.0,    0.2,     0.0],
        [131.0,  60.0,    -0.2,    0.0],
        [140.0,  40.0,    0.8,     0.0],
        [144.0,  60.0,    0.0,     2.0],
        [171.0,  60.0,    0.25,    0.25],
        [179.1,  102.0,   1.3,     1.0],
        [188.8,  60.0,    0.1,     0.1],
        [195.1,  80.0,    0.23,    0.08],
        [196.7,  80.0,    0.23,    0.09],
        [203.5,  120.0,   0.07,    0.06],
        [209.8,  70.0,    0.6,     0.2],
        [212.2,  70.0,    0.6,     0.2],
        [222.0,  120.0,   0.5,     0.4],
        [230.0,  70.0,    0.0,     0.3],
        [232.0,  70.0,    0.0,     -0.2],
        [234.5,  70.0,    0.8,     0.3],
        [246.0,  400.0,   0.07,    0.1],
        [250.4,  23.0,    3.3,     0.0],
        [500.0,  1.0,     8.0,     10.0]
    ];

    data_dict = {}

    mechMode_pattern = {
        "MechMode":{
            "type": "MechMode",
            "name": "MechMode",
            "f0": {
                "value": 30e3,
                "units": "Hz",
                "description": "Frequency"
            },
            "Q": {
                "value": 5.0,
                "units": "Unitless",
                "description": "Quality factor"
            },
            "full_scale": {
                "value": 1.13,
                "units": "sqrt(J)",
                "description": "FPGA full-scale resonator amplitude register"
            }
        }
    }

    elecMode_dict = {
        "ElecMode0": {
            "mech_couplings" : {
                "value": {"MechMode0": 0.0, "MechMode1": 0.0, "MechMode2": 0.0}
            }
        }
    }

    piezo_dict = {
        "Piezo0": {
            "mech_couplings" : {
                "value": {"MechMode0": 0.0, "MechMode1": 0.0, "MechMode2": 0.0}
            }
        }
    }

    elecMode_couplings = {}
    piezo_couplings = {}

    for i, line in enumerate(data_list):

        mechMode_dict = copy.deepcopy(mechMode_pattern)
        key, value = mechMode_dict.popitem()

        new_key = key + str(i)
        mechMode_dict[new_key] = value
        mechMode_dict[new_key]["name"] = new_key
        mechMode_dict[new_key]["f0"]["value"] = line[0]
        mechMode_dict[new_key]["Q"]["value"] = line[1]

        data_dict.update(mechMode_dict)
        elecMode_couplings[new_key] = line[2]
        piezo_couplings[new_key] = line[3]

    elecMode_dict["ElecMode0"]["mech_couplings"]["value"] = elecMode_couplings
    piezo_dict["Piezo0"]["mech_couplings"]["value"] = piezo_couplings

    data_dict.update(elecMode_dict)
    data_dict.update(piezo_dict)

    print json.dumps(data_dict, indent=4, separators=(',',':'), sort_keys=True)


if __name__=="__main__":
    data_to_JSON()
