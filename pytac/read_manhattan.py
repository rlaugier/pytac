

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
from pytac import get_TF

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



def GTF_plot(master_time_s, sig_a, sig_b, nperseg=2e4, force_indices=False,
            threshold=0.9, mask_fmin=10., mask_fmax=500., 
            showfig=True, saveit=False,
            save_prefix="plots/f_corr_", save_suffix=f".pdf",
            amplitude_lims="default",
            labels=None,
            DL_TF_ph_coeff=np.nan,
            DL_TF_amp_power=np.nan,
            dl_tf_amp=np.nan,
            ref_freqs=[],
            amplims=None,
            flims=None,
            phlims=None,
            display_gain=False):
    if flims is None:
        flims = (1.0, 300.0)
    if phlims is None:
        phlims = (-4.0,4.0)
    if force_indices is not False:
        indices = force_indices
    else:
        indices = np.arange(sig_a.shape[1])
    if labels is None:
        labels = indices
    if flims is None:
        flims = (10., 500.)
    if amplims is None:
        amplims = [None, None]
    elif amplims == "default":
        amplims = (1e-2, 1e2)
    master_dt = np.gradient(master_time_s).mean()
    coh_FT_vib = []
    for i in indices:
        coh = sig.coherence(sig.detrend(sig_a[:,i]), sig.detrend(sig_b[:,i]),
                            fs=1/master_dt, nperseg=nperseg)
        coh_FT_vib.append(coh)
    coh_FT_vib = np.array(coh_FT_vib).T

    alltfs = np.array([get_TF(master_time_s, sig_a[:,u], sig_b[:,u]*1e-6,
                                           get_coh=True, nps=nperseg) for u in indices])

    f1, TFsig, pos_coherence = alltfs[0,0,:], alltfs[:,1,:].T, alltfs[:,2,:].T
    coh_mask = (pos_coherence>=threshold) * (f1<=mask_fmax)[:,None] * (f1>=mask_fmin)[:,None]
    TFsig2 = np.where(coh_mask, TFsig, np.nan)
    coh_2 = np.where(coh_mask, pos_coherence, np.nan)
    dewrap = np.pi* np.cumsum(np.gradient(np.angle(TFsig), axis=0)>=3., axis=0)

    unwrap_phases = np.angle(TFsig) - dewrap

    #f_coherences, pos_coherence = sig.coherence(DL_positions[:,u],total_DL_commands[:,u], fs=1/master_dt, nperseg=1e3)

    target_freq = 150
    #mask_ampli = np.abs(mfft_DL_commands).mean(axis=1) >= 5e-5 
    res = np.abs(f1 - target_freq)
    
    freqmask = (f1>mask_fmin) * (f1<mask_fmax)
    meangain_unknown = np.median(np.abs(TFsig2[freqmask[:,None] * coh_mask]))
    meanphase_unknown = np.median(np.angle(TFsig2[freqmask[:,None] * coh_mask]))
    print(f"Mean gain unaccounted for: {meangain_unknown:.2e}")

    phase_on_target = np.angle(TFsig[np.argmin(res),:])
    amp_on_target = np.abs(TFsig[np.argmin(res),:])

    for i in indices:
        print(f"For {labels[i]} $A =$ {amp_on_target[i]:.2e}, $\\varphi = $ {phase_on_target[i]:.1f} rad")
        print(f"Delay ata {target_freq} Hz = {(phase_on_target[i]/(2*np.pi)*1/target_freq)}")
    delay_func = (unwrap_phases)/(2*np.pi)*1/f1[:,None]
    delay_func_2 = delay_func[coh_mask]
    #meandelay_unknown = np.abs(np.median(delay_func_2))
    meandelay_unknown = 0.78e-3 
    print(f"Mean phase unaccounted for: {meanphase_unknown:.2e}")
    print(f"Mean delay unaccounted for: {meandelay_unknown*1e3:.2f} ms, coming from FT integration")
    print(np.abs(TFsig))
    #print(TFsig2)
    print(indices)
    print(master_dt)
    print(master_time_s)


    fig = plt.figure(figsize=(8,8),dpi=150)
    
    # Plotting amplitude
    plt.subplot(311)
    plt.title(save_prefix)
    for i in indices:
        plt.plot(f1, np.abs(TFsig2[:,i]), color=f"C{i}", label=f"{labels[i]}", marker="+",
                 markersize=2., markeredgewidth=0.4, linewidth=0.2)
        plt.plot(f1, np.abs(TFsig[:,i]), color=f"C{i}", alpha=0.2, linewidth=0.5)
    for afreq in ref_freqs:
        plt.axvline(afreq, linewidth=0.2, color="k")
    plt.axhline(1., linewidth=0.5, color="k")
    if display_gain:
        plt.text(2*flims[0], meangain_unknown, f"Extra gain: {meangain_unknown:.3f}")
        plt.axhline(meangain_unknown, linestyle="--", linewidth=0.5, color="k")
    # plt.plot(f1, dl_tf_amp**DL_TF_amp_power, color="k", alpha=0.5, label="Typical DL amp")
    plt.xlim(*flims)
    plt.ylabel("Ampllitude")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim(*amplims)
    plt.legend(fontsize="x-small", loc="upper right")
    
    # Plotting phase
    plt.subplot(312)
    plt.plot(f1, -2*np.pi*f1 * meandelay_unknown, "k--",
             label=f"Expected from {meandelay_unknown*1e3:.2f} ms\n of delay")
    for i in indices:
        plt.plot(f1, np.angle(TFsig2[:,i]), color=f"C{i}", label=f"{labels[i]}", marker="+",
                 markersize=2., markeredgewidth=0.4, linewidth=0.2)
        plt.plot(f1, np.angle(TFsig[:,i]), color=f"C{i}", alpha=0.2, linewidth=0.5)
    plt.axhline(0, linewidth=0.5, color="k")
    for afreq in ref_freqs:
        plt.axvline(afreq, linewidth=0.2, color="k")
    # plt.plot(f1, DL_TF_ph_coeff * np.angle(dl_tf), color="k", alpha=0.5, label="Typical DL phase")
    plt.xlim(*flims)
    plt.ylim(*phlims)
    plt.ylabel("Phase [rad]")
    plt.xscale("log")
    plt.legend(fontsize="xx-small", loc="upper right")
    
    # Plotting coherence
    plt.subplot(313)
    for i in indices:
        plt.plot(f1, coh_2[:,i], color=f"C{i}", label=f"{labels[i]}", marker="+",
                 markersize=2., markeredgewidth=0.4, linewidth=0.2)
        plt.plot(f1, pos_coherence[:,i], color=f"C{i}", alpha=0.1)
    for afreq in ref_freqs:
        plt.axvline(afreq, linewidth=0.2, color="k")
    plt.xlim(*flims)
    plt.xlabel("Freq [Hz]")
    plt.ylabel("Coherence")
    plt.xscale("log")
    plt.legend(fontsize="x-small", loc="lower right")
    #plt.tight_layout()
    if saveit:
        plt.savefig(f"{save_prefix}_{save_suffix}", bbox_inches="tight")
    if showfig:
        plt.show()
    return fig
