#!/usr/bin/python

#
# Routines for pretty-printing data structures for linacs etc
#
# Alejandro F Queiruga
# Daniel Scott Driver
# 2013 LBL
#

import linac

####################################
#
# Converting a Linac_Param and children to string
#
####################################


def filtostr(fil):

    fil.A_modes = linac.intArray_frompointer(fil.modes)
    fil.A_coeff_start = linac.intArray_frompointer(fil.coeff_start)
    fil.A_coeffs = linac.complexdouble_Array_frompointer(fil.coeffs)
    fil.A_poles = linac.complexdouble_Array_frompointer(fil.poles)

    ret = "  Filter of Order: {0}\n".format(fil.order)
    for o in xrange(fil.order):
        ret+="  Modes for {0}: poles    a    b    scale\n".format(o)
        for m in xrange(fil.A_modes[o]):
            cs = fil.A_coeff_start[o]+m
            ret+="   "+str(fil.A_poles[cs])+", "+ \
                str(fil.A_coeffs[3*cs+0])+", "+\
                str(fil.A_coeffs[3*cs+1])+", "+\
                str(fil.A_coeffs[3*cs+2])+"\n"
        #ret+="\n"
    return ret

def cavtostr(cav):
    return """  psd_llrf: {0}
  w0: {1}
  bunch_rep: {2}
  Q_L: {3}
  R_Q: {4}
  beta_in: {5}
  beta_out: {6}
  beta_beam: {7}

  bandw: {8}
  noise_rms: {9}
  bw_ol: {10}
  k: {11}
  
  nom_beam_phase: {12}
  rf_phase: {13}
  design_voltage: {14}
  unity_voltage: {15}
""" .format(cav.psd_llrf,cav.w0,cav.bunch_rep, cav.Q_L,
            cav.R_Q,cav.beta_in,cav.beta_out,cav.beta_beam,
            cav.bandw,cav.noise_rms,cav.bw_ol,cav.k,
            cav.nom_beam_phase,cav.rf_phase,cav.design_voltage,
            cav.unity_voltage)

def fpgatostr(fpga):
    return """  gain: {0}
  int_gain: {1}
  set_point: {2}
""".format(fpga.gain,fpga.int_gain,fpga.set_point)

def lintostr(lin):
    return """---Beam Paramters:---
  dE: {0}
  R56: {1}
  T566: {2}
  phi: {3}
  lam: {4}
  s0: {5}
  a: {6}
  L: {7}
""".format(lin.dE,lin.R56,lin.T566,lin.phi,
                  lin.lam,lin.s0,lin.a,lin.L) + \
                  "\n---TRF1:---\n" + filtostr(lin.TRF1) + \
                  "\n-Saturate_c: {0}\n".format(lin.saturate_c) + \
                  "\n---TRF2:---\n" + filtostr(lin.TRF2) + \
                  "\n---Cavity:---\n" + cavtostr(lin.cav) + \
                  "\n---Cavity Filter:---\n" + filtostr(lin.Cav_Fil) + \
                  "\n---RXF:---\n" + filtostr(lin.RXF) + \
                  "\n---FPGA:---\n" + fpgatostr(lin.fpga)





####################################
#
# Converting a Filter into a graphviz image
#
####################################

def filter_dotify(fil,fname="test.png"):
    import pydot
    graph = pydot.Dot(graph_type='digraph',rankdir="LR")
    innode = pydot.Node("Input")
    graph.add_node(innode)
    las = innode
    for o in xrange(fil.order):
        if fil.A_modes[o]>1:
            summer = pydot.Node("Sigma{0}".format(o),label="&Sigma;",
                                shape="circle")
            graph.add_node(summer)

        for m in xrange(fil.A_modes[o]):
            nod = pydot.Node(str(fil.A_poles[fil.A_coeff_start[o]+m]))
            graph.add_node(nod)
            #these += [nod]
            if fil.A_modes[o]>1:
                graph.add_edge(pydot.Edge(nod,summer))
            edg = pydot.Edge(las,nod)
            graph.add_edge(edg)
        if fil.A_modes[o]>1:
            las=summer
        else:
            las=nod

    ounode = pydot.Node("Output")
    graph.add_node(ounode)
    graph.add_edge(pydot.Edge(las,ounode))
    graph.write_png(fname)
