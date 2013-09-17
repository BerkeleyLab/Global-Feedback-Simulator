#!/usr/bin/python

# shim for plotting without X server
import matplotlib as mtp
mtp.use('Agg')

import sys
import numpy as np
import pylab as py

def plotdata(plotcfgfile):

    import sys
    sys.path.append("../")
    import json
    from readjson.readjson import jsontodict
    from readjson.readjson import readentry
    import numpy as np
    import pylab as py

    plotdict = jsontodict(plotcfgfile) #read json and return python dictionary
    acceldict = jsontodict(plotdict["Accelerator Config"]) #read accelerator config

    datafile=plotdict["Data"]
    connect=acceldict['Accelerator']['connect']
    Nlinac=len(connect)

    closeflag=plotdict.get('closewindow','none')
    display= plotdict.get('display','True') in ['True','true']
    py.close(closeflag)

    for (typeplot,plotcont) in plotdict.iteritems():
        plotflag=True
        if type(plotcont) != dict or "type" not in plotcont:
            continue

        figurenum=int(plotcont.get('figure',1))
        output_fn=plotcont.get('filename',None)
        local_display=plotcont.get('display',None)
        if (local_display):
            display = local_display in ['True','true']

        if plotcont["type"]=='TF':
            in_col=iodict_to_column_num(plotcont["input"],connect=connect)
            out_col=iodict_to_column_num(plotcont["output"],connect=connect)

            data=np.loadtxt(datafile,dtype=np.complex,usecols=(in_col,out_col))
            indata=data[:,0]
            outdata=data[:,1]

            scale_out=float(plotcont["output"].get("scale",1.0))
            scale_in=float(plotcont["input"].get("scale",1.0))
            OL_suppression=float(plotcont.get("OL_suppression",0.0))

            dt=float(acceldict["Simulation"]["dt"])
            effective_dt=dt*float(acceldict['Simulation']['Outputfreq'])

            wintype=plotcont.get('windowtype',None)
            steadyN=int(readentry(plotcont,plotcont.get('steadyN',0),acceldict['Simulation']))

            py.figure(figurenum)
            errcode=TF_plot(indata,outdata,effective_dt,
                    scale_in=scale_in,scale_out=scale_in,
                    OL_suppression=OL_suppression,wintype=wintype,
                    steadyN=steadyN)
            if output_fn:
                py.savefig(output_fn, bbox_inches=0)
            if display:
                py.show()

        elif (plotcont["type"]=='versus'):
            py.figure(figurenum)
            errcode=versus_plot(datafile,plotcont,connect=connect)
            if output_fn:
                py.savefig(output_fn, bbox_inches=0)
            if display:
                py.show()
        else:
            print 'Case for that type {0}  is not present'.format(plotcont["type"])

    return 0

def iodict_to_column_num(iodict,connect=None):
    if connect:
        linidx={connect[i]:i for i in range(len(connect))}

    linackeys={
        "error_vol_a":0,
        "error_vol_p":1,
        "dE_E":2,
        "dtz":3,
        "sz":4,
        "dE_Ei":5,
        "dE_Ei2":6,
        "cav_voltage":7,
        "fpga_err":8,
        "fpga_set_point":9,
        "fpga_drive":10,
        "fpga_state":11
        }
    otherkeys={
        "t":0,
        "Q":1,
        "dQ_Q":2,
        "dtg":3,
        "adc_noise":4
        }

    item=iodict["quantity"]
    if item in otherkeys:
        col=otherkeys[item]
    elif item in linackeys:
        if "linac" not in iodict:
            print "need a linac key word for this item"
            return False
        if (iodict["linac"]=="last" or iodict["linac"]=="Last"):
            linnum=len(linidx)-1
        else:
            linnum=linidx[iodict["linac"]]

        col=len(otherkeys)+linnum*len(linackeys) \
            +linackeys[item]
    else:
        print "that key does not have a translation"
        return None
    print col
    return col


def versus_plot(datafile,plotcont,connect=None):
    skiprows=int(plotcont.get('skiprows',0))

    y_col=iodict_to_column_num(plotcont["y"],connect=connect)
    scale_y=float(plotcont["y"].get("scale",1.0))
    ylabelis=plotcont.get('ylabel',"{0} times {1}".format(plotcont['y']['quantity'],scale_y))


    if(plotcont.get('x',False)):
        x_col=iodict_to_column_num(plotcont["x"],connect=connect)
        scale_x=float(plotcont["x"].get("scale",1.0))
        xlabelis=plotcont.get('xlabel',"{0} times {1}".format(plotcont['x']['quantity'],scale_x))

        usecols=(x_col,y_col)
        data=np.loadtxt(datafile,dtype=np.complex,usecols=usecols,skiprows=skiprows)
        x=data[:,0]
        y=data[:,1]

    else:
        print ' no x data input plotting y against data number'
        scale_x=1.0;
        xlabelis="data series number"

        usecols=(y_col,)
        y=np.loadtxt(datafile,dtype=np.complex,usecols=usecols,skiprows=skiprows)
        x=np.arange(0,y.size)


    xtoplot=x*scale_x
    ytoplot=y*scale_y

    linetype=plotcont.get('linetype','')
    linewidth=plotcont.get('linewidth',1)
    linelabel=plotcont.get('linelabel',None)

    py.plot(xtoplot,ytoplot,linetype,label=linelabel,linewidth=linewidth)
    py.xlabel(xlabelis)
    py.ylabel(ylabelis)
    py.legend()

    return 0


def TF_plot(indata,outdata,dt,
            scale_in=1.0,scale_out=1.0, OL_suppression=0.0,
            wintype=None,steadyN=0):

    from scipy import fft
    from scipy import signal

    #cut out the beginning transient data if specificied
    trunc_indata=indata[steadyN:]
    trunc_outdata=outdata[steadyN:]

    #use a windowing function in scipy.signal to condition
    #data for fft. No window if not specified in file
    if(wintype):
        windcommand="signal.{0}(len(trunc_indata))".format(wintype)
        window=eval(windcommand)
        trunc_indata=trunc_indata*window
        trunc_outdata=trunc_outdata*window

    #perform the ffts
    fft_in=fft(trunc_indata*scale_in)
    fft_out=fft(trunc_outdata*scale_out)

    fft_rat = fft_out/fft_in;
    #form the magnitude and phase data and scale appropiately
    TF_mag=20*np.log10(np.abs(fft_rat))
    TF_mag=TF_mag-OL_suppression

    TF_pha = np.angle(fft_rat);
    N=fft_in.size
    freq=np.arange(N)/(N*dt)

    N_nyq=np.int(N/2) # only frequencys up to half the sampling frequency are
                      #viable so divide in half and round down to cut out
                      #negative(aliased?) results above 1/(2*dt)

    py.subplot(2,1,1)
    py.semilogx(freq[0:N_nyq],TF_mag[0:N_nyq])
    py.xlabel('frequency [log]')
    py.ylabel('Magnitude')
    py.subplot(2,1,2)
    py.semilogx(freq[0:N_nyq],TF_pha[0:N_nyq]*180/np.pi)
    py.xlabel('frequency [log]')
    py.ylabel('Angle [deg]')

    return 0

# Make the plot autoexecutable
if __name__=="__main__":
    configfile= sys.argv[1]
    plotdata(configfile)
