import numpy as np
import oct2py

fftout = oct2py.octave.call('load', "results_TF_charge_dE3.dat")

fftx = fftout["ffts"]["fft_x"][0]
ffty = fftout["ffts"]["fft_y"][0]

rawout = oct2py.octave.call('load', "SSToutput.dat")
rawout = rawout[-6000:, :]

for n in range(rawout.shape[1]):
    rawfft = oct2py.octave.call('fft', rawout[:, n])
    print n
    print max(abs(fftx.real-rawfft.real[0]))


freq_domain = fftout["ffts"]["freq_domain"][0]
