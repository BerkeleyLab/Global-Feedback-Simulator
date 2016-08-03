#!/usr/bin/python

# Shim for plotting without X server
# Force matplotlib to not use any X-Windows back-end.
# (use function needs to be called before importing pyplot)
# import matplotlib as mtp
# mtp.use('Agg')

import sys
import numpy as np
# import matplotlib.pyplot as plt

def PlotData(plotcfgfile, columnscfgfile, XWindows=False):
    """ PlotData: Plot data given a configuration file specifying plot properties
    and a file describing the format of the data."""

    from readjson.readjson import jsontodict, readentry
    from readjson.parse_simulation import loadConfig

    import matplotlib as mtp

    # Shim for plotting without X server
    # Force matplotlib to not use any X-Windows back-end.
    # (use function needs to be called before importing pyplot)
    if XWindows == False:
        mtp.use('Agg')

    import matplotlib.pyplot as plt

    plotdict = jsontodict(plotcfgfile)  # Read JSON and return python dictionary
    acceldict = loadConfig(plotdict["Accelerator Config"])  # Read accelerator configuration
    columnsdict = jsontodict(columnscfgfile)  # Read column mapping

    datafile = plotdict["Data"]
    linac_connect = acceldict['Accelerator']['linac_connect']
    Nlinac = len(linac_connect)

    closeflag = plotdict.get('closewindow', 'none')
    display = plotdict.get('display', 'True') in ['True', 'true']
    plt.close(closeflag)

    # Iterate over items of plot configuration
    for (typeplot, plotcont) in plotdict.iteritems():

        plotflag = True

        # Pass if the content of a given entry is not a dictionary (i.e. not a plot description)
        if type(plotcont) != dict or "type" not in plotcont:
            continue

        # Figure number
        figurenum = int(plotcont.get('figure', 1))
        # Output file name
        output_fn = plotcont.get('filename', None)
        local_display = plotcont.get('display', display)

        # If the plot is a Transfer Function (TF)
        if plotcont["type"] == 'TF':

            # Input/Output column index
            in_col = Get_Column_Idx(plotcont["input"], columnsdict, linac_connect)
            out_col = Get_Column_Idx(plotcont["output"], columnsdict, linac_connect)

            # Load data
            data = np.loadtxt(datafile, dtype=np.double, usecols=(in_col, out_col))

            # Data vectors
            indata = data[:, 0]  # Input
            outdata = data[:, 1]  # Output

            # Scaling and suppression
            scale_out = float(plotcont["output"].get("scale", 1.0))
            scale_in = float(plotcont["input"].get("scale", 1.0))
            OL_suppression = float(plotcont.get("OL_suppression", 0.0))

            # Simulation time-step
            Tstep = float(acceldict['Simulation']['Tstep']['value'])
            effective_Tstep = Tstep*float(acceldict['Simulation']['OutputFreq']['value'])

            # Windowing
            wintype = plotcont.get('windowtype', None)
            steadyN = int(readentry(plotcont, plotcont.get('steadyN', 0), acceldict['Simulation']))

            # Get entry describing the input
            indict = plotcont['input']

            inlabel = columnsdict[indict['quantity']]['name']

            # Add Linac prefix if Input quantity is not Global
            if not columnsdict[indict['quantity']].get('global', False):
                inlabel = indict['linac'] + ': ' + inlabel

            # Get entry describing the output
            outdict = plotcont['output']
            outlabel = columnsdict[outdict['quantity']]['name']

            # Add Linac prefix if Output quantity is not Global
            if not columnsdict[outdict['quantity']].get('global', False):
                outlabel = outdict['linac'] + ': ' + outlabel

            # Plot
            plt.figure(figurenum)

            errcode = TF_plot(indata, outdata, inlabel, outlabel, effective_Tstep,
                              scale_in=scale_in, scale_out=scale_in,
                              OL_suppression=OL_suppression, wintype=wintype,
                              steadyN=steadyN, XWindows=XWindows)

            # End of Transfer Function (TF) plot

        # "Versus" plot
        elif (plotcont["type"] == 'versus'):
            plt.figure(figurenum)
            errcode = Versus_Plot(datafile, plotcont, columnsdict, linac_connect, XWindows=XWindows)
        # Only TF or versus types supported
        else:
            print 'Case for that type {0}  is not present'.format(plotcont["type"])
            continue

        # Save figure if output file is specified
        if output_fn:
            plt.savefig(output_fn, bbox_inches=0)

        # Show if local_display is True
        if local_display:
            plt.show()

    return 0

def Get_Column_Idx(iodict, columndict, linac_connect=None):
    """Get_Column_Idx: Get the column index position for a given quantity.
    Inputs:
        iodict: Input/Output dictionary describing the entry,
        columndict: Dictionary containing data file column format,
        linac_connect: List of Linacs (if applicable).
    Output:
        col: Column index for quantity described by iodict."""

    # Associate Linac names with their corresponding index
    if linac_connect:
        linidx = {linac_connect[i]: i for i in range(len(linac_connect))}

    # Quantity name
    code = iodict["quantity"]
    # Find item in dictionary describing data file format
    column = columndict[code]
    # If entry is not found in description return None
    if not column:
        print "that key does not have a translation"
        return None

    # If its a Global quantity then index corresponds to index value
    elif column.get("global", False):
        col = int(column["index"])
    # Else, it needs to be offset by the Linac number
    else:
        if "linac" not in iodict:
            print "need a linac key word for this item"
            return False
        if (iodict["linac"] == "last" or iodict["linac"] == "Last"):
            linnum = len(linidx)-1
        else:
            linnum = linidx[iodict["linac"]]

        col = int(columndict["globals"])+linnum*int(columndict["locals"]) \
            + int(column["index"])

    # Return column number
    return col

def Versus_Plot(datafile, plotcont, columndict, linac_connect=None, XWindows=False):
    """Versus_Plot: Generate plot of Output VS. Input signals.
        Inputs:
            datafile: File containing waveforms from the Simulation,
            plotcont: Plot content information,
            columndict: Data file column format,
            linac_connect: List of Linacs in Accelerator (if applicable)."""

    import matplotlib as mtp

    # Shim for plotting without X server
    # Force matplotlib to not use any X-Windows back-end.
    # (use function needs to be called before importing pyplot)
    if XWindows == False:
        mtp.use('Agg')

    import matplotlib.pyplot as plt

    skiprows = int(plotcont.get('skiprows', 0))

    # Grab information about the Y quantity
    y_col = Get_Column_Idx(plotcont["y"], columndict, linac_connect)
    y_code = plotcont['y']['quantity']
    y_dict = columndict[y_code]
    y_name = y_dict.get('name', y_code)

    # Add Linac prefix to name if Y quantity is not Global
    if not y_dict.get('global', False):
        y_name = plotcont['y']['linac'] + ': ' + y_name

    # Grab more properties
    y_units = y_dict.get('units', '')
    scale_y = float(plotcont["y"].get("scale", 1.0))

    # Determine Y label for plot
    if scale_y != 1.0:
        ylabelis = plotcont.get('ylabel', "{0} times {1} [{2}]".format(y_name, scale_y, y_units))
    else:
        ylabelis = plotcont.get('ylabel', '{0} [{1}]'.format(y_name, y_units))

    # If X quantity is defined then format, otherwise define it as index of Y
    if plotcont.get('x', False):

        # Grab information about the X quantity as it was done for Y
        x_col = Get_Column_Idx(plotcont["x"], columndict, linac_connect)
        x_code = plotcont['x']['quantity']
        x_dict = columndict[x_code]
        x_name = x_dict.get("name", x_code)

        if not x_dict.get('global', False):
            x_name = plotcont['x']['linac'] + ': ' + x_name
        x_units = x_dict.get('units', '')
        scale_x = float(plotcont["x"].get("scale", 1.0))
        if scale_x != 1.0:
            xlabelis = plotcont.get('xlabel', "{0} times {1} [{2}]".format(x_name, scale_x, x_units))
        else:
            xlabelis = plotcont.get('xlabel', "{0} [{1}]".format(x_name, x_units))

        usecols = (x_col, y_col)
        data = np.loadtxt(datafile, dtype=np.double, usecols=usecols, skiprows=skiprows)
        x = data[:, 0]
        y = data[:, 1]

        dimension_x = plotcont["x"].get("dimension", None)
        if dimension_x == "amp":
            xtoplot = np.abs(xtoplot)
            xlabelis += " (Amplitude)"
            x_name += " (Amplitude)"
        elif dimension_x == "phase":
            xtoplot = np.angle(xtoplot)
            xlabelis += " (Phase)"
            x_name += " (Phase)"
    else:
        print ' no x data input plotting y against data number'
        scale_x = 1.0
        xlabelis = "data series number"
        x_name = "Index"

        usecols = (y_col,)
        y = np.loadtxt(datafile, dtype=np.complex, usecols=usecols, skiprows=skiprows)
        x = np.arange(0, y.size)

    xtoplot = x*scale_x
    ytoplot = y*scale_y

    dimension_y = plotcont["y"].get("dimension", None)
    if dimension_y == "amp":
        ytoplot = np.abs(ytoplot)
        ylabelis += " (Amplitude)"
        y_name += " (Amplitude)"
    elif dimension_y == "phase":
        ytoplot = np.angle(ytoplot)
        ylabelis += " (Phase)"
        y_name += " (Amplitude)"

    linetype = plotcont.get('linetype', '')
    linewidth = plotcont.get('linewidth', 1)
    linelabel = plotcont.get('linelabel', None)

    if plotcont.get('x', False):
        if plotcont['x']['quantity'] == 't':
            title = plotcont.get('title', y_name)
        else:
            title = plotcont.get('title', '{0} vs. {1}'.format(y_name, x_name))
    else:
        title = plotcont.get('title', '{0} vs. {1}'.format(y_name, x_name))

    plt.xlim([xtoplot.min(), xtoplot.max()])
    plt.plot(xtoplot, ytoplot, linetype, label=linelabel, linewidth=linewidth)
    plt.xlabel(xlabelis)
    plt.ylabel(ylabelis)
    plt.title(title)

    return 0

def TF_plot(indata, outdata, indata_label, outdata_label, Tstep,
            scale_in=1.0, scale_out=1.0, OL_suppression=0.0,
            wintype=None, steadyN=0, XWindows=False):

    from scipy import fft
    from scipy import signal

    import matplotlib as mtp

    # Shim for plotting without X server
    # Force matplotlib to not use any X-Windows back-end.
    # (use function needs to be called before importing pyplot)
    if XWindows == False:
        mtp.use('Agg')

    import matplotlib.pyplot as plt

    # Cut out transient data at the start  if specified
    trunc_indata = indata[steadyN:]
    trunc_outdata = outdata[steadyN:]

    # Use a windowing function in sciplt.signal to condition
    # data for FFT (no windowing if not specified in file)
    if(wintype):
        windcommand = "signal.{0}(len(trunc_indata))".format(wintype)
        window = eval(windcommand)
        trunc_indata = trunc_indata*window
        trunc_outdata = trunc_outdata*window

    # Perform the FFTs
    fft_in = fft(trunc_indata*scale_in)
    fft_out = fft(trunc_outdata*scale_out)

    # Add an infinitesimal to avoid division by 0
    fft_rat = (np.abs(fft_out) + np.spacing(1))/(np.abs(fft_in) + np.spacing(1))

    # Form the magnitude and phase data and scale appropriately
    TF_mag = 20*np.log10(fft_rat)
    TF_mag = TF_mag-OL_suppression

    TF_pha = (np.angle(fft_out/fft_in))

    N = fft_in.size
    freq = np.arange(N)/(N*Tstep)

    # Only frequencies up to half the sampling frequency are
    # present, so divide in half and round down to cut out
    # aliased results above 1/(2*Tstep)
    N_nyq = np.int(N/2)

    plt.subplot(2, 1, 1)
    plt.semilogx(freq[0:N_nyq], TF_mag[0:N_nyq])
    plt.ylabel('Power [dB]')
    plt.title('Transfer Function: {0} vs. {1}'.format(indata_label, outdata_label))
    plt.subplot(2, 1, 2)
    plt.semilogx(freq[0:N_nyq], TF_pha[0:N_nyq]*180/np.pi)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Phase [deg]')

    return 0

# Make the plot auto-executable
if __name__ == "__main__":
    configfile = sys.argv[1]
    if len(sys.argv) > 2:
        columns = sys.argv[2]
    else:
        columns = 'source/plotting/columns.json'
    PlotData(configfile, columns)
