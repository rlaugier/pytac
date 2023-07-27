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
tf_man_file = st.file_uploader("upload tf_manhattan.p", type=["p"])
if tf_123_file is not None:
    tf_123 = pickle.load(tf_123_file)
else:
    exit(0)
if tf_man_file is not None:
    tf_man = pickle.load(tf_man_file)
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

def show_info(anhdu):
    with st.expander("More informations"):
        st.write(f"Date of acquisition = {hdul[0].header['DATE-OBS']}")
        st.text([anhdu[0].header])
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


with st.sidebar:
    global_nps = st.number_input("nps (number of samples per segment)", value=10000, step=100)
st.title("Vibration input-output exploration")

include_ft = st.checkbox("include FT", value=False)


tab_settings, tab_global,\
tab_mirror, tab_combinations = st.tabs(["Settings", "Global view",
                                                    "Mirror view",
                                                    "Combination"])
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

    st.header("Uploading a file")

    st.write("Upload a file")

    myfile = st.file_uploader("Choose a file", accept_multiple_files=False, type=["fits"])
    myfilename = None
    if myfile is not None:
        myfilename = myfile.name.split(".")[0]
        hdul = fits.open(myfile)
        show_info(hdul)
    else:
        st.write("No file uploaded")
        exit()


    crop_samples = st.number_input(label="Number of samples discarted from \
                                        each ends (to enable synchronization)",
                    min_value=0,
                    value=10,
                    step=1,)
    
    # Creating a master time to resample everything
    timeref = st.selectbox("Reference for the time", options=pt.UT_names)
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
    show_info(anhdu)

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
                label=f"Positions {selected_mah[i]}")
        detrend_commands = sig.detrend(total_DL_commands[:,i_dl], axis=0)
        plt.plot(master_time_s, detrend_commands,
                color=f"C{i}", linewidth=0.5, linestyle="--",
                label=f"Commands {selected_dl[i]}")
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
    mirrors = np.arange(1, 8+1)
    # mirror_mask = []
    # for i, acol in enumerate(st.columns(len(mirrors))):
    #     with acol:
    #         mirror_mask.append(st.checkbox(label=f"M{mirrors[i]}", key=f"global_check_M_{i+1}"))
    # mirror_mask = np.array(mirror_mask)
    # st.write(mirror_mask)
    st.header("Manhattan positions")
    with st.expander("Plot options"):
        flims = []
        ylims = [None,None]
        cols_f = st.columns(2)
        with cols_f[0]:
            flims.append(st.number_input("f start [Hz]", value=1., step=10.))
            y_mode_log = st.checkbox("Log scale", value=False, key="man_pos_log" )
            if y_mode_log:
                 ylims[0] = st.number_input("y min [m]", value=1.0e-9, step=0.5e-7, format="%.e")
            else:
                ylims[0] = 0
        with cols_f[1]:
            flims.append(st.number_input("F end [Hz]", value=300., step=10.))
            ylims[1] = st.number_input("y max [m]", value=2.0e-7, step=0.5e-7, format="%.e")
        freq_cols = st.columns(4)
        ref_freqs = []
        for i, acol in enumerate(freq_cols):
            with acol:
                ref_freqs.append(st.number_input("Frequency input", min_value=0., key=f"frequency_ref_1_{i}"))
            
    with st.sidebar:
        mirror_names_selection = st.multiselect("Which mirrors to include", mirrors, [4])
        mirror_mask = np.array([(amirror in mirror_names_selection) for amirror in mirrors])

        ut_names_selection = st.multiselect("Which mirrors to include", selected_mah, selected_mah)
        telescope_mask = np.array([(aut_name in ut_names_selection) for aut_name in selected_mah])

    mirror_filtered_pos_4567 = []
    for i in base_indices:
        aresp = control.forced_response(tf_123, T=master_time_s, U=all_mirrors[:,mirror_mask>0,i].sum(axis=1))
        mirror_filtered_pos_4567.append(aresp.y[0,:])
    mirror_filtered_pos_4567 = np.array(mirror_filtered_pos_4567).T

    f0, psd_mirrors_opl = sig.welch(mirror_filtered_pos_4567, fs=1/master_dt, nperseg=global_nps, axis=0)
    df = np.mean(np.gradient(f0))
    # psd_rcum = np.cumsum(psd_mirrors_opl[::-1], axis=0)[::-1]*df
    psd_rcum = np.sqrt(np.cumsum((psd_mirrors_opl[::-1] * df), axis=0))[::-1]
    
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

    st.header("Signal presence")
    

st.header("Transfer function")

indices = np.arange(4)
st.write("list_indices")
st.write(list_indices)
alltfs = np.array([pytac.get_TF(master_time_s, total_DL_commands[:,i_tel], DL_positions[:,i_dl],
                                       get_coh=True) for i_tel, i_dl in zip(base_indices, base_indices)])

f1, TFsig, pos_coherence = alltfs[0,0,:], alltfs[:,1,:].T, alltfs[:,2,:].T
dewrap = np.pi* np.cumsum(np.gradient(np.angle(TFsig), axis=0)>=3., axis=0)

unwrap_phases = np.angle(TFsig) - dewrap

#f_coherences, pos_coherence = sig.coherence(DL_positions[:,u],total_DL_commands[:,u], fs=1/master_dt, nperseg=1e3)

st.write("Pick a frequency of interest for the cursor")
target_freq = st.slider("Target frequency", min_value=10., max_value=500.,value=150.,
                        step=0.1)
res = np.abs(f1 - target_freq)

phase_on_target = unwrap_phases[np.argmin(res),:]
amp_on_target = np.abs(TFsig[np.argmin(res),:])

my_columns = st.columns(len(list_indices))
for i, (u, col, k) in enumerate(zip(list_indices, my_columns, base_indices)):
    with col:
        st.write(f"## {pt.all_dl_names[selected_dl_indices[i]]}")
        st.write(f"From **{pt.UT_names[u]}**")
        st.write(f"at {target_freq:.1f}Hz")
        st.write(f"$A = $ {amp_on_target[i]:.2e},")
        st.write(f"$\\varphi = $ {phase_on_target[i]:.1f} rad")
        st.write(f"Delay at {target_freq} Hz = {(phase_on_target[i]/(2*np.pi)*1/target_freq):.4f}s")
delay_func = (unwrap_phases)/(2*np.pi)*1/f1[:,None]


flims = (1., 2000)
fig_tf_1 = plt.figure(figsize=(8,8))
plt.subplot(311)
plt.title("Transfer function")
for i, u in enumerate(list_indices):
    plt.plot(f1, np.abs(TFsig[:,i]), label=f"DL{selected_dl[i]}({selected_mah[i]})")
    st.write()
# plt.plot(lf, ltf_a, label=f"Used last", linestyle="--")
#plt.plot(shift(master_freqs), shift(np.abs(fft_DL_pos)), label="DL_pos")
#plt.fill_between(shift(master_mfft_freqs), 1e-2, 1e2, where=mask_ampli, alpha=0.1)
plt.axvline(target_freq, linewidth=0.5, color="k")
plt.axhline(1., linewidth=0.5, color="k")
plt.xlim(flims[0], flims[1])
plt.ylim(1e-3, 2.)
#plt.xlabel("Freq [Hz]")
plt.ylabel("Ampllitude")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.subplot(312)
for i, u in enumerate(list_indices):
    plt.plot(f1, unwrap_phases[:,i], label=f"DL{selected_dl[i]}({selected_mah[i]})")
# plt.plot(lf, ltf_ph, linestyle="--", label=f"Used last")
#plt.plot(shift(master_freqs), shift(np.abs(fft_DL_pos)), label="DL_pos")
#plt.fill_between(shift(master_mfft_freqs), -4., 4., where=mask_ampli, alpha=0.1)
plt.axhline(0, linewidth=0.5, color="k")
plt.axvline(target_freq, linewidth=0.5, color="k")
for i, u in enumerate(list_indices):
    plt.axhline(phase_on_target[i], linewidth=0.5, linestyle="--", color=f"C{i}")
plt.xlim(flims[0], flims[1])
#plt.ylim(-12, 1)
#plt.xlabel("Freq [Hz]")
plt.ylabel("Phase [rad]")
#plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.subplot(313)
for i, u in enumerate(list_indices):
    plt.plot(f1, pos_coherence[:,i], label=f"DL{selected_dl[i]}({selected_mah[i]})")
plt.axvline(target_freq, linewidth=0.5, color="k")
plt.xlim(flims[0], flims[1])
plt.xlabel("Freq [Hz]")
plt.ylabel("Coherence")
#plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.tight_layout()
plt.show()



fig_tf_2 = plt.figure(figsize=(8,3))
for i, u in enumerate(list_indices):
    plt.plot(f1, delay_func[:,i], label=f"DL{selected_dl[i]}({selected_mah[i]})")
marks = (-2.2e-3, -2.6e-3)
markstyles = ("--", ":")
for amark, astyle in zip(marks, markstyles):
    plt.axhline(amark, linewidth=1, linestyle=astyle, color="k")
    plt.text(2., amark, f"{amark:.2e} s")
plt.legend(fontsize="x-small")
plt.xscale("log", )
plt.xlim(flims[0], flims[1])
plt.ylim(-5e-3, 0)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Delay [s]")
plt.show()

st.pyplot(fig_tf_1)
st.pyplot(fig_tf_2)


output_columns = st.columns(len(selected_dl))
for i, acol in enumerate(output_columns):
    with acol:
        st.write(f"DL{selected_dl[i]}")
        good_co = pos_coherence[:,i] >= 0.4
        clean_f = np.abs(f1[good_co])
        clean_TF = TFsig[good_co]
        clean_TF_amp = (np.abs(TFsig))[good_co]
        clean_TF_ph = (np.angle(TFsig) - dewrap)[good_co]

        with io.BytesIO() as buffer:
            np.savetxt(fname=buffer, 
                       header="#Measured DL_transfer function, f[Hz]",
                      X=clean_f,
                      delimiter=";")
            st.download_button(label="Download the frequency file, f[Hz]", data=buffer,
                                file_name="clean_f.dat")

        with io.BytesIO() as buffer:
            np.savetxt(fname=buffer, 
                       header="#Measured DL_transfer function, Amp, cpx",
                      X=clean_TF,
                      delimiter=";")
            st.download_button(label="Download function, Amp, cpx", data=buffer,
                                file_name="clean_TF.dat")
            
        with io.BytesIO() as buffer:
            np.savetxt(fname=buffer, 
                       header="#Measured DL_transfer function, Amp, real",
                      X=clean_TF_amp,
                      delimiter=";")
            st.download_button(label="Download function, Amp, real", data=buffer,
                                file_name="clean_amp.dat")

        with io.BytesIO() as buffer:
            np.savetxt(fname=buffer, 
                       header="#Measured DL_transfer function, Phase, [rad]",
                      X=clean_TF_ph,
                      delimiter=";")
            st.download_button(label="Download function, Phase[rad], real", data=buffer,
                                file_name="clean_ph.dat")

