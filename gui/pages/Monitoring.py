import streamlit as st

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import pytac
from pathlib import Path
from sys import exit

import scipy.signal as sig
from scipy.interpolate import interp1d
import io
# from pytac_vlti_context import *
import pytac_vlti_context as pt
import control
import pickle
tf_123_file = st.file_uploader("upload tf_123.p", type=["p"])
if tf_123_file is not None:
    tf_123 = pickle.load(tf_123_file)
else:
    exit(0)

def resample(x, y, master, be=True):
    """Just a macro for interp1d"""
    values = interp1d(x, y, fill_value=np.nan, bounds_error=be, )(master)
    return values
def all_valid(data):
    """Data is a list of arrays,
    time dimension is first index"""
    isbad = np.zeros(data[0].shape[0], dtype=bool)
    for adata in data:
        print(adata.shape)
        np.logical_or(isbad, adata)
    return isbad

def show_and_save(figure, name):
    st.pyplot(figure)
    # if st.button("Save the plot", key=f"button_save{i}"):
    format = st.selectbox("Format to save", options=["pdf", "png"], key=f"format_{name.strip()}")
    if format is not None:
        with io.BytesIO() as buffer:
            figure.savefig(buffer, bbox_inches="tight",format="pdf", dpi=200)
            st.download_button("Download pdf", data=buffer,
                        file_name=f"{name.strip()}.pdf",
                        key=f"dl_button{name.strip()}")

pfreqs = np.array([83.0,
        84.5,
        85.0,
        87.0,
        90.0,
        95.,
        100.,
        110.,
        120.,
        166.0,
        175.0,
        200.])
amps_1 = np.array([9.9e-8,
                9.25e-8,
                8.6e-8,
                4.5e-8,
                1.62e-8,
                5.3e-9,
                2.5e-9,
                9.0e-10,
                4.3e-10,
                5.65e-11,
                4.3e-11,
                 2.19e-11])

all_PLLs = [[{'FINIT': 150.0,
               'FMIN': 145.0,
               'FMAX': 155.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.13333333333333333,
               'name': 'PLL_1',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 173.0,
               'FMIN': 168.0,
               'FMAX': 178.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.11560693641618497,
               'name': 'PLL_2',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 83.5,
               'FMIN': 78.5,
               'FMAX': 88.5,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.23952095808383234,
               'name': 'PLL_3',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 47.8,
               'FMIN': 44.8,
               'FMAX': 50.8,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.4184100418410042,
               'name': 'PLL_4',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1}],
             [{'FINIT': 127.5,
               'FMIN': 125,
               'FMAX': 145,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 10.0,
               'PDTAU': 0.1568627450980392,
               'name': 'PLL_1',
               'input_block': 'Fork_TF1',
               'gain': 0.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 100.0,
               'FMIN': 96.0,
               'FMAX': 104.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 0.1,
               'PDTAU': 0.2,
               'name': 'PLL_2',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 83.5,
               'FMIN': 80.5,
               'FMAX': 86.5,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.23952095808383234,
               'name': 'PLL_3',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 42.0,
               'FMIN': 37.0,
               'FMAX': 47.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.47619047619047616,
               'name': 'PLL_4',
               'input_block': 'Fork_TF1',
               'gain': 0.4,
               'phase': -2.0,
               'sum_input': 1}],
             [{'FINIT': 150.0,
               'FMIN': 145.0,
               'FMAX': 155.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.13333333333333333,
               'name': 'PLL_1',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 210.0,
               'FMIN': 202.0,
               'FMAX': 218.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.09523809523809523,
               'name': 'PLL_2',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 81.0,
               'FMIN': 73.0,
               'FMAX': 89.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.24691358024691357,
               'name': 'PLL_3',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 100.0,
               'FMIN': 96.0,
               'FMAX': 104.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.2,
               'name': 'PLL_4',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1}],
             [{'FINIT': 150.0,
               'FMIN': 140,
               'FMAX': 155,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 0.05,
               'PDTAU': 0.13333333333333333,
               'name': 'PLL_1',
               'input_block': 'Fork_TF1',
               'gain': 0.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 200.0,
               'FMIN': 195.0,
               'FMAX': 205.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.1,
               'name': 'PLL_2',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 100.0,
               'FMIN': 95.0,
               'FMAX': 105.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.2,
               'name': 'PLL_3',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1},
              {'FINIT': 75.0,
               'FMIN': 71.0,
               'FMAX': 79.0,
               'KP': 2.0,
               'KI': 0.5,
               'KLP': 1.0,
               'PDTAU': 0.26666666666666666,
               'name': 'PLL_4',
               'input_block': 'Fork_TF1',
               'gain': 1.0,
               'phase': 0.0,
               'sum_input': 1}]]

x_filter, y_filter = pfreqs/83., amps_1/np.max(amps_1)
band_filter = interp1d(x_filter, y_filter, kind="linear", bounds_error=False, fill_value=(1., 1e-4))
def filter_peak(f, amp, f0, band_filter):
    # print(np.abs(f/f0))
    amp_out = amp * band_filter(1. + np.abs((f-f0)/f0))
    return amp_out


amps_2 = np.array([1.98e-7,
                    1.85e-7,
                    1.72e-7,
                    9.0e-8,
                    3.24e-8,
                    1.065e-8,
                    5.05e-9,
                    1.8e-9,
                    8.65e-10,
                    1.13e-10,
                    8.55e-11,
                    4.375e-11,])

def test_pll_peaks(peaks, PLL_list, rel_amp_margin = 0.05, verbose=False):
    for apll in PLL_list:
        pll_info = f"-> **{apll['name']}**\n\n{apll['FMIN']} to {apll['FMAX']} Hz\n\n"
        pll_active = True
        if apll["gain"] == 0.:
            pll_active = False
        peak_distances = np.abs(peaks[:,0] - apll["FINIT"])
        peak_rel_distances = peak_distances / (apll["FMAX"] - apll["FMIN"])
        peak_centered = peak_rel_distances <= 0.1
        peak_present = peak_rel_distances <= 1.
        if np.count_nonzero(peak_present) == 0:
            if pll_active:
                st.warning(f"{pll_info}warning: no peak found in the interval")
            else:
                st.info(f"{pll_info}Deactivated\n\nInfo: no peak found in the interval")
        elif (np.count_nonzero(peak_present) == 1) and (np.count_nonzero(peak_centered) == 0):
            st.warning(f"{pll_info}Warning: peak at {peaks[peak_present][0]} poorly centered")
        if np.count_nonzero(peak_present) > 1:
            multi_peak=True
            st.error(f"{pll_info}Warning: extra peak in the filter's band:")
            st.write("( f [Hz], a [µm] )")
            st.write(peaks[peak_present])
        else:
            multi_peak = False

        for i, apeak in enumerate(peaks):
            peak_freq = apeak[0]
            peak_amp = apeak[1]
            amp_out = filter_peak(peak_freq, peak_amp, apll["FINIT"], band_filter)
            relative_amp = amp_out/peak_amp
            if peak_present[i]:
                if multi_peak:
                    st.error(f"{pll_info} amp={amp_out:.2f} for {peak_amp:.2f}")
                    st.write(f"f = {peak_freq:.2f} This is assumed to peak targeted at {apll['FINIT']:.2f}")
                else:
                    st.success(f"{pll_info} amp={amp_out:.2f} for {peak_amp:.2f}")
                    st.write(f"f = {peak_freq:.2f} This is assumed to peak targeted at {apll['FINIT']:.2f}")
                
            else:
                if verbose:
                    st.write(f"Relative amplitude = {relative_amp:.2f}")
                if relative_amp >= rel_amp_margin: 
                    st.error(f"Warning: relative amplitude margin exceeded")


with st.sidebar:
    global_nps = st.number_input("nps (number of samples per segment)", value=10000, step=100)
st.title("Vibration input-output exploration")

include_ft = st.checkbox("include FT", value=True)

st.header("Uploading a file")
myfile = st.file_uploader("Choose a file", accept_multiple_files=False, type=["fits"])
myfilename = None

tab_settings, tab_global,\
tab_mirror, tab_gen_tf = st.tabs(["Settings", "Global view",
                                                    "Mirror view",
                                                    "Generalized TF"])
with tab_settings:

    st.header("Telescope - DL pairing")
    columns = st.columns(len(pt.UT_names))
    selected_mah = []
    selected_dl = []
    input_mask = []
    for i_col, acol in enumerate(columns):
        with acol:
            st.write(f"input {i_col}")
            maskit = st.checkbox("Activate", value=True, key=f"mask_{i_col}")
            input_mask.append(maskit)
            my_selected_tel = st.selectbox(label="MAH :", options=pt.UT_names, key=f"tel_{i_col}",
                                            index=i_col, disabled=(not maskit))
            if my_selected_tel is None:
                my_selected_tel = pt.UT_names[i_col]
            if maskit is False:
                my_selected_tel = None
            else:
                selected_mah.append(my_selected_tel)
            if maskit:
                my_selected_dl = st.selectbox(label="DL :", options=pt.all_dl,
                    index=pt.ut2dl[pt.ut_names2indices[my_selected_tel]] - 1,
                    key=f"dl_{i_col}", disabled=(not maskit))
                selected_dl.append(my_selected_dl)
        
    ut2dl_new = {a:b for a,b in zip(selected_mah, selected_dl)}
    list_indices = [pt.ut2ind[pt.ut_names2indices[amah]] for amah in selected_mah]
    selected_dl_indices = [pt.all_dl_names2indices[f"DL{anindex}"] for anindex in selected_dl]
    list_dl_names = [f"DL{anindex}" for anindex in selected_dl]
    base_indices = np.arange(len(selected_dl_indices))
    # selected_dl_indices = [pt.dl_number2indices[adl] for adl in selected_dl]

    with st.expander("Show debug lists of selected ports", ):
        st.text("selected_mah")
        st.write(selected_mah)
        st.text("selected_dl")
        st.write(selected_dl)
        st.text("list_indices")
        st.write(list_indices)
        st.text("ut2dl_new")
        st.write(ut2dl_new)
        st.text("selected_dl_indices")
        st.write(selected_dl_indices)
        # st.write(selected_dl_indices)
        st.write("DL_names")
        st.write(pt.DL_names)
        st.write(list_dl_names)

    if myfile is not None:
        myfilename = myfile.name.split(".")[0]
        hdul = fits.open(myfile)
        anhdu = hdul
        with st.expander("More informations"):
            st.write(f"Date of acquisition = {hdul[0].header['DATE-OBS']}")
            st.text([anhdu[0].header])

            crop_samples = st.number_input(label="Number of samples discarted from \
                                                each ends (to enable synchronization)",
                            min_value=0,
                            value=10,
                            step=1,)
    
            # Creating a master time to resample everything
            timeref = st.selectbox("Reference for the time", options=pt.UT_names)
    else:
        st.write("No file uploaded")
        exit()

    with st.spinner("Processing the data"):
        onetime = hdul[timeref].data["TIME"]
        master_time = np.arange(onetime[crop_samples], np.max(onetime), 250.)[:-crop_samples]
        master_time_s = master_time * 1e-6
        master_dt = np.gradient(master_time*1e-6).mean()
        del onetime


        anhdu = hdul
        sensors = np.arange(12)
        if include_ft:
            #Resampling DL_OFFSET
            FT_DL_commands = np.array([resample(anhdu["OPDC"].data["TIME"], anhdu["OPDC"].data["VLTI_DL_OFFSET"][:,i],
                                            master_time) for i in list_indices]).T
        # Resampling from other source TTR-110.0016
        MAN_DL_commands =  np.array([resample(anhdu[f"MAH-{aname}"].data["TIME"], anhdu[f"MAH-{aname}"].data["DPL"],
                                            master_time) for aname in list_dl_names]).T
        #Resampling DL positions
        DL_positions =  np.array([resample(anhdu[aname].data["TIME"], anhdu[aname].data["POS"],
                                            master_time) for aname in list_dl_names]).T
        total_DL_commands = MAN_DL_commands # - FT_DL_commands + MAN_DL_commands 

        if include_ft:
            #Resampling DL_OFFSET
            FT_OPD_raw = np.array([resample(anhdu["OPDC"].data["TIME"], anhdu["OPDC"].data["OPD"][:,i],
                                                master_time) for i in list_indices]).T
        if include_ft:
            opdc = anhdu["OPDC"].data
            t = opdc['TIME']
            opd = opdc['OPD']
            kopd = opdc['KALMAN_OPD']
            kpiezo = opdc['KALMAN_PIEZO']
            vlti = opdc['VLTI_DL_OFFSET']
            unwrapped_opd = ((opd - kopd + np.pi) % (2*np.pi) - np.pi + kopd)
            FT_UW_OPD = np.array([resample(anhdu["OPDC"].data["TIME"], unwrapped_opd[:,i],
                                                master_time) for i in range(6)]).T

            raw_POL = 2.2 / (2*np.pi) * (unwrapped_opd + (pt.T2B @ (kpiezo - vlti*2*np.pi/2.2e-6).T).T)
            FT_POL = np.array([resample(anhdu["OPDC"].data["TIME"], raw_POL[:,i],
                                                master_time) for i in range(6)]).T

            raw_UT_POL = pt.Ap.dot(raw_POL.T).T 
            UT_POL = np.array([resample(anhdu["OPDC"].data["TIME"], raw_UT_POL[:,i],
                                                master_time) for i in range(4)]).T

        all_sensors = []
        all_mirrors = []
        for aname in pt.MAN_names:
            raw_sensors = np.array([resample(anhdu[aname].data["TIME"], anhdu[aname].data["RAW"][:,i],
                                            master_time) for i in sensors])
            all_sensors.append(raw_sensors)

            raw_mirrors = []
            raw_mirrors.append(-2*pt.volt2acc*raw_sensors[pt.mirror_indices[1]].sum(axis=0)/4)
            raw_mirrors.append(2*pt.volt2acc*raw_sensors[pt.mirror_indices[2]].sum(axis=0))
            raw_mirrors.append(np.sqrt(2)*pt.volt2acc*raw_sensors[pt.mirror_indices[3]].sum(axis=0)/2)
            raw_mirrors.append(np.sqrt(2)*pt.volt2acc*raw_sensors[pt.mirror_indices[4]].sum(axis=0))
            raw_mirrors.append(1.9941*pt.volt2acc*raw_sensors[pt.mirror_indices[5]].sum(axis=0))
            raw_mirrors.append(1.8083*pt.volt2acc*raw_sensors[pt.mirror_indices[6]].sum(axis=0))
            raw_mirrors.append(1.9820*pt.volt2acc*raw_sensors[pt.mirror_indices[7]].sum(axis=0))
            raw_mirrors.append(2*pt.volt2acc*raw_sensors[pt.mirror_indices[8]].sum(axis=0))
            all_mirrors.append(raw_mirrors)
        all_sensors = np.array(all_sensors).T
        all_mirrors = np.array(all_mirrors).T

        sos = sig.butter(4, 5., btype="high", fs= 1/master_dt, output="sos")
        better_mirrors = sig.sosfiltfilt(sos, all_mirrors, axis=0)
            

    
with tab_global:
    if st.checkbox("Show power spectrum", value=True):
        my_nps = st.slider("nperseg", value=int(1e3), step=1,
                        min_value=10, max_value=int(1.0e4))
        fig = plt.figure()
        for i in base_indices:
            plt.plot(*sig.welch(DL_positions[:,i], fs=4000, nperseg=my_nps))
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Power spectrum [m^2 /Hz]")
        plt.show()

        st.pyplot(fig)

    st.header("Time series")
    with st.expander("Plot tweaks"):
        st.write("For the display only:")

        time_start = st.slider("Start value [s]", min_value=0., max_value=np.max(master_time_s),
                        value=0., key="time_start_slider")

        time_end = time_start + st.number_input("Length [s]", min_value=0., max_value=np.max(master_time_s),
                        value=1.0, key="time_end_slider")
        enable_local_extrema = st.checkbox("Local y min max")
        amax, amin = np.nan, np.nan
        time_mask = (master_time_s>time_start) * (master_time_s<time_end)
        mask_fig = plt.figure(figsize=(8,0.2))
        plt.plot(master_time_s, time_mask)
        st.pyplot(mask_fig)
    fig_temporal = plt.figure(dpi=200)
    for i, (i_tel, i_dl) in enumerate(zip(base_indices, base_indices)):
        detrend_pos = sig.detrend(DL_positions[:,i_tel], axis=0)
        plt.plot(master_time_s, detrend_pos,
                color=f"C{i}", linewidth=1, alpha=0.5,
                label=f"Positions DL{selected_dl[i]}")
        detrend_commands = sig.detrend(total_DL_commands[:,i_dl], axis=0)
        plt.plot(master_time_s, detrend_commands,
                color=f"C{i}", linewidth=0.5, linestyle="--",
                label=f"Commands {selected_mah[i]}")
        amax = np.max([np.nan_to_num(amax),
            np.max(detrend_pos[time_mask]),\
            np.max(detrend_commands[time_mask])])

        amin = np.min([np.nan_to_num(amin),
            np.min(detrend_pos[time_mask]),
            np.min(detrend_commands[time_mask])])
    plt.xlabel("Time [s]")
    plt.ylabel("Position [m]")
    plt.xlim(time_start, time_end)
    if enable_local_extrema:
        plt.ylim(amin, amax)
    plt.legend(fontsize="x-small")
    plt.show()
    show_and_save(fig_temporal, "fig_temporal")
    # st.pyplot(fig_temporal)

    # ######################################################
with tab_mirror:
    st.header("Monitoring")
    verbose = st.checkbox("Verbose")
    peaks = np.array([[47.5, 5.],
                      [82.,10.],
                      [90.0, 8.0],
                      [100.,8.]])
    tel_columns = st.columns(4)
    for i_col, (acol, pll_list) in enumerate(zip(tel_columns, all_PLLs)):
        with acol:
            st.write(f"### UT{i_col+1}")
            test_pll_peaks(peaks, pll_list, verbose=verbose)
        
    
    mirrors = np.arange(1, 8+1)
    # mirror_mask = []
    # for i, acol in enumerate(st.columns(len(mirrors))):
    #     with acol:
    #         mirror_mask.append(st.checkbox(label=f"M{mirrors[i]}", key=f"global_check_M_{i+1}"))
    # mirror_mask = np.array(mirror_mask)
    # st.write(mirror_mask)
    st.header("Manhattan positions")
    with st.sidebar:
        with st.expander("Plot options"):
            flims = []
            ylims = [None,None]
            cols_f = st.columns(2)
            with cols_f[0]:
                flims.append(st.number_input("f start [Hz]", value=40.0, step=10.))
                y_mode_log = st.checkbox("Log scale", value=False, key="man_pos_log" )
                if y_mode_log:
                     ylims[0] = st.number_input("y min [m]", value=1.0e-9, step=0.5e-7, format="%.e")
                else:
                    ylims[0] = 0
            with cols_f[1]:
                flims.append(st.number_input("F end [Hz]", value=300., step=10.))
                ylims[1] = st.number_input("y max [m]", value=1.0e-7, step=0.5e-7, format="%.e")
            # freq_cols = st.columns(4)
            ref_freqs = []
            for i, acol in enumerate(range(4)):
                # with acol:
                    ref_freqs.append(st.number_input("Frequency input", min_value=0., key=f"frequency_ref_1_{i}"))
            
    with st.sidebar:
        mirror_names_selection = st.multiselect("Which mirrors to include", mirrors, [4,5,6,7])
        mirror_mask = np.array([(amirror in mirror_names_selection) for amirror in mirrors])

        ut_names_selection = st.multiselect("Which mirrors to include", selected_mah, selected_mah)
        telescope_mask = np.array([(aut_name in ut_names_selection) for aut_name in selected_mah])

#############################
# Important computations here
    mirror_filtered_pos_4567 = []
    for i in base_indices:
        aresp = control.forced_response(tf_123, T=master_time_s, U=all_mirrors[:,mirror_mask>0,i].sum(axis=1))
        mirror_filtered_pos_4567.append(aresp.y[0,:])
    mirror_filtered_pos_4567 = np.array(mirror_filtered_pos_4567).T

    f0, psd_mirrors_opl = sig.welch(mirror_filtered_pos_4567, fs=1/master_dt, nperseg=global_nps, axis=0)
    df = np.mean(np.gradient(f0))
    # psd_rcum = np.cumsum(psd_mirrors_opl[::-1], axis=0)[::-1]*df
    psd_rcum = np.sqrt(np.cumsum((psd_mirrors_opl[::-1] * df), axis=0))[::-1]
    vib_on_BL_full_4567 = pt.UT2B.dot(mirror_filtered_pos_4567.T).T
    #vib_on_BL_full_4567 = vib_on_BL_full_3
    MAN_DL_commands_on_BL = pt.UT2B.dot(MAN_DL_commands.T).T
    mirror_filtered_pos_3 = []
    for i in range(4):
        aresp = control.forced_response(tf_123, T=master_time_s, U=all_mirrors[:,:3,i].sum(axis=1))
        mirror_filtered_pos_3.append(aresp.y[0,:])
    mirror_filtered_pos_3 = np.array(mirror_filtered_pos_3).T



    mirror_filtered_pos = []
    for i in range(4):
        aresp = control.forced_response(tf_123, T=master_time_s, U=all_mirrors[:,:-1,i].sum(axis=1))
        mirror_filtered_pos.append(aresp.y[0,:])
    mirror_filtered_pos = np.array(mirror_filtered_pos).T


    vib_on_BL_full = pt.UT2B.dot(mirror_filtered_pos.T).T
    vib_on_BL_full_3 = pt.UT2B.dot(mirror_filtered_pos_3.T).T


# Most computation are done
    
    fig_all_manhattan = plt.figure(dpi=200)
    for i, mah_name in enumerate(selected_mah):
        if telescope_mask[i]:
            plt.plot(f0, np.sqrt(psd_mirrors_opl[:,i]),
                     linewidth=1., color=f"C{i}", label=f"{mah_name}")
            plt.plot(f0, psd_rcum[:,i], color=f"C{i}", label=f"{mah_name}")
        plt.xlim(*flims)
    for afreq in ref_freqs:
        plt.axvline(afreq, linewidth=0.5, color="k")
    plt.ylabel("PSD [m/Hz^1/2]\n rcumsum [m]")
    plt.grid(visible=True)
    plt.yticks(np.arange(0., 2.0e-6, 0.2e-7))
    #plt.yscale("log")
    man_fig_title = f"{myfilename}_{mirror_mask.astype(int)}"
    plt.title(man_fig_title)
    plt.xscale("log")
    if y_mode_log:
        plt.yscale("log")
    plt.ylim(*ylims)
    plt.legend(fontsize="x-small", loc="upper right")
    # plt.savefig(f"plots/levels_{filepath.stem}_{mirror_mask}.pdf", bbox_inches="tight")
    show_and_save(fig_all_manhattan, f"ut{telescope_mask.astype(int)}{myfilename}_{mirror_mask.astype(int)}")

    st.header("Waterfall")
    if st.checkbox("Compute waterfall"):
        n_nps = st.selectbox("N_nps", [1,2,3,4,5,6,7,8,9,10], index=4)
        # wf_data = st.selectbox("Dataset", [mirror_filtered_pos, mirror_filtered_pos_3, mirror_filtered_pos_4567,
        #                                 vib_on_BL_full, vib_on_BL_full_3, vib_on_BL_full_4567], index=2)
        # wf_labels = st.selectbox("Dataset", [selected_mah, selected_dl], index=0)

        wf_data = mirror_filtered_pos_4567
        wf_labels = selected_mah
        allnps = np.linspace(500, global_nps, n_nps)
        fig_waterfall = pytac.plot_waterfall(wf_data, *flims,
                                data_names=selected_mah, fs=1/master_dt,
                                allnps=allnps, figsize=(20, 20), dpi=100, colorbar=True,
                                vmin=ylims[0], vmax=ylims[1])
        show_and_save(fig_waterfall, f"Waterfall_{myfilename}")

    st.header("Signal presence")
    
    tel_index = st.selectbox("Select telescope ", options=base_indices,)
    f1, sig_vib_pos = sig.welch(mirror_filtered_pos_4567[:,tel_index], fs=1/master_dt, nperseg=global_nps)
    f1, sig_vib_commands = sig.welch(MAN_DL_commands[:,tel_index], fs=1/master_dt, nperseg=global_nps)
    f1, sig_dl_pos = sig.welch(DL_positions[:,tel_index], fs=1/master_dt, nperseg=global_nps)
    
    fig_presence = plt.figure(dpi=200)
    plt.subplot(211)
    plt.title(f"Detecting the presence of output signal to DL {selected_mah[tel_index]}")
    plt.plot(f1, sig_vib_pos, label="A. Integrated acc.")
    plt.plot(f1, sig_vib_commands, label="B. Commands")
    plt.plot(f1, sig_dl_pos, label="C. DL position" )
    for afreq in ref_freqs:
        plt.axvline(afreq, linewidth=0.5, color="k")
    plt.legend(fontsize="x-small")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("PSD")
    plt.ylim(1e-22, 1e-14)
    plt.xlim(*flims)
    plt.subplot(212)
    plt.plot(*sig.coherence(mirror_filtered_pos_4567[:,tel_index], MAN_DL_commands[:,tel_index], fs=1/master_dt, nperseg=global_nps), label="A<->B")
    plt.plot(*sig.coherence(mirror_filtered_pos_4567[:,tel_index], DL_positions[:,tel_index], fs=1/master_dt, nperseg=global_nps), label="A<->C")
    for afreq in ref_freqs:
        plt.axvline(afreq, linewidth=0.5, color="k")
    plt.legend(fontsize="x-small")
    plt.xscale("log")
    plt.xlim(*flims)
    plt.ylabel("Coherence")
    plt.xlabel("Frequency [Hz]")
    plt.ylim(0., 1.)
    # plt.savefig(f"plots/UT{tel_index+1}_presence_{myfilename}.pdf", dpi=120, bbox_inches="tight")
    plt.show()
    show_and_save(fig_presence, name=f"Presence_{myfilename}")
    
    st.header("Gain computation (intersting for ON files)")
    if st.checkbox("Compute gains", value=False):
        gain_pos_command = sig_vib_commands / sig_vib_pos
        gain_pos_dl_pos = sig_dl_pos / sig_vib_pos
        fig_gain = plt.figure(dpi=200)
        plt.axhline(1., color="k", linewidth=0.5)
        plt.plot(f1, gain_pos_command, label=f"Gain A-> B")
        plt.plot(f1, gain_pos_dl_pos, label=f"Gain A-> C")
        for afreq in ref_freqs:
            plt.axvline(afreq, linewidth=0.5, color="k")
        plt.legend(fontsize="x-small")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(*flims)
        plt.ylabel("Coherence")
        plt.xlabel("Frequency [Hz]")
        plt.ylim(1e-1, 1e1)
        plt.show()
        show_and_save(fig_gain, name=f"Gain_{myfilename}")

        for afreq in ref_freqs:
            myindex = np.argmin(np.abs(f1-afreq))
            st.write(f"f = {afreq:.1f}, gain A-> B = {gain_pos_command[myindex]:.2f}, gain B-> A = {gain_pos_dl_pos[myindex]:.2f}")

