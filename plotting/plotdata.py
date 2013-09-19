#!/usr/bin/python

# shim for plotting without X server
import matplotlib as mtp
mtp.use('Agg')

import sys
import numpy as np
import pylab as py

def plotdata(plotcfgfile, columnscfgfile):

    import sys
    sys.path.append("../")
    import json
    from readjson.readjson import jsontodict
    from readjson.readjson import readentry
    import numpy as np
    import pylab as py

    plotdict = jsontodict(plotcfgfile) #read json and return python dictionary
    acceldict = jsontodict(plotdict["Accelerator Config"]) #read accelerator config
    columnsdict = jsontodict(columnscfgfile) # read column mapping

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
        local_display=plotcont.get('display',display)

        if plotcont["type"]=='TF':
            in_col=iodict_to_column_num(plotcont["input"],columnsdict,connect)
            out_col=iodict_to_column_num(plotcont["output"],columnsdict,connect)

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

            indict = plotcont['input']
            inlabel = columnsdict[indict['quantity']]['name']
            if not columnsdict[indict['quantity']].get('global',False):
                inlabel = indict['linac'] + ': ' + inlabel
            outdict = plotcont['output']
            outlabel = columnsdict[outdict['quantity']]['name']
            if not columnsdict[outdict['quantity']].get('global',False):
                outlabel = outdict['linac'] + ': ' + outlabel
       
            py.figure(figurenum)
            errcode=TF_plot(indata,outdata,inlabel, outlabel, effective_dt,
                    scale_in=scale_in,scale_out=scale_in,
                    OL_suppression=OL_suppression,wintype=wintype,
                    steadyN=steadyN)

        elif (plotcont["type"]=='versus'):
            py.figure(figurenum)
            errcode=versus_plot(datafile,plotcont,columnsdict,connect)
        else:
            print 'Case for that type {0}  is not present'.format(plotcont["type"])
            continue
        
        if output_fn:
            py.savefig(output_fn, bbox_inches=0)
        if local_display:
            py.show()

    return 0

def iodict_to_column_num(iodict,columndict,connect=None):
    if connect:
        linidx={connect[i]:i for i in range(len(connect))}

    code=iodict["quantity"]
    column = columndict[code]
    if not column:
        print "that key does not have a translation"
        return None
    elif column.get("global",False):
        col=int(column["index"])
    else :
        if "linac" not in iodict:
            print "need a linac key word for this item"
            return False
        if (iodict["linac"]=="last" or iodict["linac"]=="Last"):
            linnum=len(linidx)-1
        else:
            linnum=linidx[iodict["linac"]]

        col=int(columndict["globals"])+linnum*int(columndict["locals"]) \
            +int(column["index"])
    #print col
    return col


def versus_plot(datafile,plotcont,columndict,connect=None):
    skiprows=int(plotcont.get('skiprows',0))

    y_col=iodict_to_column_num(plotcont["y"],columndict,connect)
    y_code = plotcont['y']['quantity']
    y_dict = columndict[y_code]
    y_name = y_dict.get('name',y_code)
    if not y_dict.get('global',False):
        y_name = plotcont['y']['linac'] + ': ' + y_name

    y_units = y_dict.get('units','')
    scale_y=float(plotcont["y"].get("scale",1.0))
    if scale_y != 1.0:
        ylabelis=plotcont.get('ylabel',"{0} times {1} [{2}]".format(y_name,scale_y,y_units))
    else:
        ylabelis=plotcont.get('ylabel','{0} [{1}]'.format(y_name,y_units))

    if plotcont.get('x',False):
        x_col=iodict_to_column_num(plotcont["x"],columndict,connect)
        x_code = plotcont['x']['quantity']
        x_dict = columndict[x_code]
        x_name = x_dict.get("name",x_code)
        if not x_dict.get('global',False):
            x_name = plotcont['x']['linac'] + ': ' + x_name
        x_units = x_dict.get('units','')
        scale_x=float(plotcont["x"].get("scale",1.0))
        if scale_x != 1.0:
            xlabelis=plotcont.get('xlabel',"{0} times {1} [{2}]".format(x_name,scale_x, x_units))
        else:
            xlabelis=plotcont.get('xlabel',"{0} [{1}]".format(x_name, x_units))

        usecols=(x_col,y_col)
        data=np.loadtxt(datafile,dtype=np.complex,usecols=usecols,skiprows=skiprows)
        x=data[:,0]
        y=data[:,1]

    else:
        print ' no x data input plotting y against data number'
        scale_x=1.0
        xlabelis="data series number"

        usecols=(y_col,)
        y=np.loadtxt(datafile,dtype=np.complex,usecols=usecols,skiprows=skiprows)
        x=np.arange(0,y.size)

    xtoplot=x*scale_x
    ytoplot=y*scale_y

    dimension_x = plotcont["x"].get("dimension",None)
    if dimension_x == "amp":
        xtoplot = np.abs(xtoplot)
        xlabelis += " (Amplitude)"
        x_name += " (Amplitude)"
    elif dimension_x == "phase":
        xtoplot = np.angle(xtoplot)        
        xlabelis += " (Phase)"
        x_name += " (Phase)"

    dimension_y = plotcont["y"].get("dimension",None)
    if dimension_y == "amp":
        ytoplot = np.abs(ytoplot)
        ylabelis += " (Amplitude)"
        y_name += " (Amplitude)"
    elif dimension_y == "phase":
        ytoplot = np.angle(ytoplot)        
        ylabelis += " (Phase)"
        y_name += " (Amplitude)"

    linetype=plotcont.get('linetype','')
    linewidth=plotcont.get('linewidth',1)
    linelabel=plotcont.get('linelabel',None)
    if plotcont['x']['quantity'] == 't':
        title=plotcont.get('title',y_name)
    else:
        title=plotcont.get('title','{0} vs. {1}'.format(y_name,x_name))

    py.xlim([xtoplot.min(),xtoplot.max()])
    py.plot(xtoplot,ytoplot,linetype,label=linelabel,linewidth=linewidth)
    py.xlabel(xlabelis)
    py.ylabel(ylabelis)
    py.title(title)
    py.legend()

    return 0


def TF_plot(indata,outdata,indata_label, outdata_label, dt,
            scale_in=1.0,scale_out=1.0, OL_suppression=0.0,
            wintype=None,steadyN=0):

    from scipy import fft
    from scipy import signal

    #steadyN = max(512, steadyN)

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

    # Add an infinitesimal to avoid division by 0
    fft_rat = (np.abs(fft_out) + np.spacing(1))/(np.abs(fft_in) + np.spacing(1))
    #fft_rat = np.abs(fft_out) / np.abs(fft_in)
    #form the magnitude and phase data and scale appropiately
    TF_mag=20*np.log10(fft_rat)
    TF_mag=TF_mag-OL_suppression

    TF_pha = (np.angle(fft_out/fft_in))
    #TF_pha = np.angle(fft_out)-np.angle(fft_in)
    N=fft_in.size
    freq=np.arange(N)/(N*dt)

    N_nyq=np.int(N/2) # only frequencys up to half the sampling frequency are
                      #viable so divide in half and round down to cut out
                      #negative(aliased?) results above 1/(2*dt)

    py.subplot(2,1,1)
    #py.xlim(freq[0], freq[N_nyq-1])
    py.semilogx(freq[0:N_nyq],TF_mag[0:N_nyq])
    #py.xlabel('Frequency [Hz]')
    py.ylabel('Power [dB]')
    py.title('Transfer Function: {0} vs. {1}'.format(indata_label,outdata_label))
    py.subplot(2,1,2)
    #py.xlim(freq[0], freq[N_nyq-1])
    py.semilogx(freq[0:N_nyq],TF_pha[0:N_nyq]*180/np.pi)
    py.xlabel('Frequency [Hz]')
    py.ylabel('Phase [deg]')

    return 0

# Make the plot autoexecutable
if __name__=="__main__":
    configfile= sys.argv[1]
    if len(sys.argv) > 2:
        columns = sys.argv[2]
    else:
        columns = '../columns.json'
    plotdata(configfile,columns)
