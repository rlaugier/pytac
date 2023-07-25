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



st.title("Transfer function computation")


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

if st.checkbox("Show debug lists of selected ports", value=False, ):
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
if myfile is not None:
    hdul = fits.open(myfile)
    st.text([hdul[0].header])
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
onetime = hdul[timeref].data["TIME"]
master_time = np.arange(onetime[crop_samples], np.max(onetime), 250.)[:-crop_samples]
master_time_s = master_time * 1e-6
master_dt = np.gradient(master_time*1e-6).mean()
del onetime


anhdu = hdul
sensors = np.arange(12)
include_ft = st.checkbox("include FT", value=False)
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

st.write(str(hdul.info()))


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

st.header("Temporal signal")

st.write("For the display only:")

time_start = st.slider("Start value", min_value=0., max_value=np.max(master_time_s),
                value=0., key="time_start_slider")

time_end = st.slider("End value", min_value=0., max_value=np.max(master_time_s),
                value=1.0, key="time_end_slider")


st.header("Time series")

fig_temporal = plt.figure(dpi=200)
for i, (i_tel, i_dl) in enumerate(zip(base_indices, base_indices)):
    st.write(i)
    plt.plot(master_time_s, sig.detrend(DL_positions[:,i_tel], axis=0),
            color=f"C{i}", linewidth=1, alpha=0.5, label=f"Positions {pt.UT_names[i_tel]}")
    st.write(pt.all_dl_names[selected_dl_indices[i_dl]])
    plt.plot(master_time_s, sig.detrend(total_DL_commands[:,i_dl], axis=0),
             color=f"C{i}", linewidth=0.5, linestyle="--", label=f"Commands {pt.DL_names[i_dl]}")
plt.xlabel("Time [s]")
plt.ylabel("Position [m]")
plt.xlim(time_start, time_end)
plt.legend()
plt.show()
st.pyplot(fig_temporal)


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

target_freq = st.slider("Target frequency", min_value=10., max_value=500.,value=150.,
                        step=0.1)
res = np.abs(f1 - target_freq)

phase_on_target = unwrap_phases[np.argmin(res),:]
amp_on_target = np.abs(TFsig[np.argmin(res),:])

st.write(f"{amp_on_target[0]:.2e}")

my_columns = st.columns(len(list_indices))
for i, (u, col, k) in enumerate(zip(list_indices, my_columns, base_indices)):
    with col:
        st.write(pt.all_dl_names)
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
        clean_f = f1[good_co]
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

