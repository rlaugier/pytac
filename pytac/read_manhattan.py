

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
#
from pathlib import Path
import scipy.signal as sig
#import scipy.interpolate import interp1d
import scipy.interpolate as interp
#import MNII_functions as mn2
#import pandas as pd

def plot_waterfall(data, fmin, fmax,
                 data_names=["Data"],
                 fs=None, colorbar=False,
                 allnps=None, show=True,
                 vmin=0, vmax=None,
                 target_freq=None,
                 **kwargs):
    """
    Plot a series of aterfall plots for different sampling intervals

    **Arguments** :
    * data : Series of data (Time at the first index)
    * fmin : minimum frequency displayd [Hz]
    * fmax : maximum frequency displayd [Hz]
    * data_names : Labels for the data channels
    * fs : Sampling frequency of the data [Hz]
    * colorbar : (default False)
    * allnps : a list of nps values to plot for
    * show : (default True)
    * vmin : (default 0)
    * vmax : (defautl None)
    * kwargs : to be passed to subplots 
    """
    print(data[:3,0])
    fig, axarr = plt.subplots(allnps.shape[0], data.shape[1] ,
                              sharex=False, sharey=True,
                              **kwargs)

    for channel_index, data in enumerate(data.T): 
        for nps_index, anps in enumerate(allnps):
            plt.sca(axarr[nps_index, channel_index])
            stf, stt, st_spectra  = sig.stft(data[:], fs=fs, nperseg=anps)
            for astt, ast_spec in zip(stt[1:-2], st_spectra.T[1:-2]):
                plt.plot(stf, np.abs(ast_spec), label=f"{astt:.2f}")
            stt = stt[1:-2]
            st_spectra = st_spectra[:,1:-2]
            fmask = (stf>=fmin) * (stf<=fmax)

            extent = (fmin, fmax, stt[-1], stt[0])
            plt.imshow(np.abs(st_spectra[fmask,:].T),extent=extent,
                       vmin=vmin, vmax=vmax,
                       interpolation="nearest")
            if target_freq is not None:
                plt.axhvline(target_freq, linewidth=0.5, color="k")
            if colorbar:
                plt.colorbar()
            plt.title(f"{data_names[channel_index]}")
            if (channel_index == 0):
                plt.ylabel(f"Time [s](nps={anps:.0f})")
            if nps_index == allnps.shape[0]:
                plt.xlabel("Frequency [Hz]")
            #plt.colorbar()
            plt.gca().set_aspect("auto")
    plt.tight_layout()
    if show:
        plt.show()
    return fig


def double_integrate_psd(f,acc_psd):
    """
        Double integration in the of a PSD in the Fourier domain: 

                **Arguments** :

                * f          : Frequencies [Hz]
                * acc_psd    : Acceleration values PSD [(m/s^2)^2 / Hz]

                **returns** : 
                * `f`    : The frequencies [Hz]
                * `psd_pos` : The PSD in position    [m^2 / Hz]

    """
    psd_pos = 1/(2*np.pi)**4 * acc_psd[1:] / f[1:]**4
    
    return( f[1:], psd_pos)

base_faulty = [[0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [2, 2, 2, 2],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0]]

def get_mosaic(file, plot_path="./",
                            plot_suffix=".png",
                            showall=True,
                            saveall=True,
                            faulty_matrix=None):
    """
            get_mosaic: 
            **Arguments** :
    * `file`       : Path to fits file
      OR hdu list
    * `plot_path`  : ('./') Path to output some plots
    * `showall`    : (True) Whether to save the plots
    * `saveall`    : (True) Whether to 

    **returns** : 
    * figs :  a list of figures
    """
    print("#################")
    print(file, flush=True)
    V2acc_gain = 0.01
    coeffs = [V2acc_gain*np.sqrt(2.0), #Acc1
              V2acc_gain*np.sqrt(2.0), #Acc2
              V2acc_gain*2.0,          #Acc3
              V2acc_gain*1.0,          #Acc4 # What is angle for M8?
              -V2acc_gain*2.0,         #Acc5
                  -V2acc_gain*2.0,         #Acc6
              -V2acc_gain*2.0,         #Acc7
              -V2acc_gain*2.0,         #Acc8
              V2acc_gain*np.sqrt(2.0), #Acc9
              V2acc_gain*1.9941,       #Acc10
              V2acc_gain*1.8083,       #Acc11
              V2acc_gain*1.9822]       #Acc12

    # 0: OK (white), 1: Faulty (Orange), 2: not-yet-installed (Grey)
    if faulty_matrix is None:
        faulty = [[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [2, 2, 2, 2],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]]
    else:
        faulty = faulty_matrix

    # %% Accelerometers [m/s^2]

    fig1, axarr1 = plt.subplots(12, 4, sharex=True, sharey=True, figsize=(12, 16))

    if isinstance(file, str):
        print("Loading from the path")
        hdul = fits.open(file)
    elif isinstance(file, list):
        print("Using file as hdul")
        hdul = file
    for tel in range(4):
        accel = hdul[f"MAH-UT{np.abs(tel+1)}"].data
        # accel = fits.getdata(filename, f'MAH-UT{np.abs(tel+1)}')
        for chan in range(12):
            axarr1[chan, tel].plot(accel['RAW'][:1000, chan]*coeffs[chan])
    hdul.close()
    axarr1[0, 0].set_ylim(-0.1, 0.1)
    for tel in range(4):
        axarr1[0, tel].set_title(f'UT{tel+1}')
        for chan in range(12):
            if faulty[chan][tel]==1:
                axarr1[chan, tel].set_facecolor('orange')
            if faulty[chan][tel]==2:
                axarr1[chan, tel].set_facecolor('grey')

    for chan in range(12):
        axarr1[chan, 0].set_ylabel(f'ACC{chan+1}')

    fig1.tight_layout()
    if saveall:
        fig1.savefig(f"{plot_path}AccelWaveform{plot_suffix}")
    if showall:
        plt.show()

    # %% Accelerometers PSD [m^2/s^4/Hz]

    fig2, axarr2 = plt.subplots(12, 4, sharex=True, sharey=True, figsize=(12, 16))

    for tel in range(4):
        axarr2[0, tel].set_title(f'UT{tel+1}')
        accel = hdul[f"MAH-UT{np.abs(tel+1)}"].data
        # accel = fits.getdata(file, f'MAH-UT{tel+1}')
        for chan in range(12):
            freq, psd = sig.welch(accel['RAW'][:, chan]*coeffs[chan], fs=4000, nperseg=2**11)
            if faulty[chan][tel]==1:
                axarr2[chan, tel].set_facecolor('orange')
            if faulty[chan][tel]==2:
                axarr2[chan, tel].set_facecolor('grey')
            axarr2[chan, tel].loglog(freq, psd)
            axarr2[chan, tel].grid()

    for chan in range(12):
        axarr2[chan, 0].set_ylabel(f'ACC{chan+1}')

    fig2.tight_layout()
    if saveall:
        fig2.savefig(f"{plot_path}AccelSpectrum{plot_suffix}")
    if showall:
        plt.show()

    # %% Positions PSD [m^2/Hz]

    fig3, axarr3 = plt.subplots(12, 4, sharex=True, sharey=True, figsize=(12, 16))

    for tel in range(4):
        axarr3[0, tel].set_title(f'UT{tel+1}')
        accel = hdul[f"MAH-UT{np.abs(tel+1)}"].data
        # accel = fits.getdata(file, f'MAH-UT{tel+1}')
        for chan in range(12):
            freq, psd = sig.welch(accel['RAW'][:, chan]*coeffs[chan], fs=4000, nperseg=2**11)
            freq, psd_pos = double_integrate_psd(freq, psd)
            if faulty[chan][tel]==1:
                axarr3[chan, tel].set_facecolor('orange')
            if faulty[chan][tel]==2:
                axarr3[chan, tel].set_facecolor('grey')
            axarr3[chan, tel].loglog(freq, psd_pos)
            axarr3[chan, tel].grid()

    for chan in range(12):
        axarr3[chan, 0].set_ylabel(f'POS{chan+1}')

    fig3.tight_layout()
    if saveall:
        fig3.savefig(f"{plot_path}PosSpectrum{plot_suffix}")
    if showall:
        plt.show()

    # %% Rev. cumulative [m^2]

    fig4, axarr4 = plt.subplots(12, 4, sharex=True, sharey=False, figsize=(12, 16))

    for tel in range(4):
        axarr4[0, tel].set_title(f'UT{tel+1}')
        accel = hdul[f"MAH-UT{np.abs(tel+1)}"].data
        # accel = fits.getdata(file, f'MAH-UT{tel+1}')
        for chan in range(12):
            freq, psd = sig.welch(accel['RAW'][:, chan]*coeffs[chan], fs=4000, nperseg=2**11)
            freq, psd_pos = double_integrate_psd(freq, psd)
            rc_pos = np.cumsum(psd_pos[::-1])[::-1] * np.diff(freq[1:])[0]
            if faulty[chan][tel]==1:
                axarr4[chan, tel].set_facecolor('orange')
            if faulty[chan][tel]==2:
                axarr4[chan, tel].set_facecolor('grey')
            axarr4[chan, tel].plot(freq, rc_pos)
            axarr4[chan, tel].set_xscale('log')
            axarr4[chan, tel].set_yscale('linear')
            axarr4[chan, tel].grid()

    for chan in range(12):
        axarr4[chan, 0].set_ylabel(f'RC{chan+1}')

    fig4.tight_layout()
    if saveall:
        fig4.savefig(f"{plot_path}RevCumSpectrum{plot_suffix}")
    if showall:
        plt.show()

    return [fig1, fig2, fig3, fig4]