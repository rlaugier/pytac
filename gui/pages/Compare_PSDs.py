import streamlit as st

import scipy.signal as sig
import astropy.io.fits as fits
import numpy as np
from pytac.read_manhattan import get_mosaic
import io
from pytac_vlti_context import rev_search, mirror_indices 
import pytac_vlti_context as ptc
import matplotlib.pyplot as plt

st.title("Creating mosaic plots")
tab_names = ["Setting",
    "Colors",
    "Comparison",
    "Comparison MC"]
tabs = st.tabs(tab_names)

def get_std(x, sig_x, y, sig_y):
    x_mc = np.random.normal(loc=x, scale=sig_x, size=(1000, x.shape[0]))
    y_mc = np.random.normal(loc=y, scale=sig_y, size=(1000, y.shape[0]))
    x_over_y = x_mc/y_mc
    amean, botq, topq = np.mean(x_over_y, axis=0),\
                        np.quantile(x_over_y,0.32, axis=0),\
                        np.quantile(x_over_y,0.68, axis=0)
    return amean, botq, topq
def get_mc(x, sig_x, y, sig_y, n_traces):
    x_mc = np.random.normal(loc=x, scale=sig_x, size=(n_traces, x.shape[0]))
    y_mc = np.random.normal(loc=y, scale=sig_y, size=(n_traces, y.shape[0]))
    x_over_y = x_mc/y_mc
    x_over_y = np.where(x_over_y>0, x_over_y,  np.nan)

    return x_over_y


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
    #          "C0",
    #          "C0",
    #          "C0",
    #          "C1",
    #          "C1",
    #          "C1",
    #          "C1"]

def plot_on_off_times(file_list, headers, colors, group_labels):
    from matplotlib.dates import DateFormatter
    time_list = []
    datetime_list = []
    for aheader in headers:
        atime = aheader["DATE"][:-5]
        time_list.append(atime)
        datetime_list.append(datetime.fromisoformat(atime))
    fig_timeline = plt.figure(figsize=(10,1))
    for i, afile in enumerate(file_list):
        mylabel = f"{afile.name} ({group_labels[i]})"
        plt.plot_date(datetime_list[i], 0., fmt="+", tz="UTC",
                    color=colors[i], markersize=4., label=mylabel)
        plt.gca().xaxis.set_major_formatter(DateFormatter('%Y-%m-%dT%H:%M:%S'))
        plt.gca().xaxis.set_tick_params(rotation=40)
    plt.legend(fontsize="x-small")
    plt.xlabel("Acquisition time (UTC)")
    return fig_timeline, time_list
    
    
with tabs[0]:
    file_list = st.file_uploader("Upload a file for mosaic",
                    accept_multiple_files=True, type=["fits"])
    
    typ = 'VIBMAH'
    # filenames = [f'./data/TTR-111.0009/on_10.fits',
    #              f'./data/TTR-111.0009/on_11.fits',
    #              f'./data/TTR-111.0009/on_12.fits',
    #              f'./data/TTR-111.0009/on_13.fits',
    #              f'./data/TTR-111.0009/all_PLLs_1.fits',
    #              f'./data/TTR-111.0009/all_PLLs_2.fits',
    #              f'./data/TTR-111.0009/all_PLLs_3.fits',
    #              f'./data/TTR-111.0009/all_PLLs_4.fits',]

    st.header("Identifying files")
    st.write("Identifiers for the data")
    on_label = st.text_input("Ticked are ", value="on")
    off_label = st.text_input("Unticked are: ", value="off")
    st.write("Pick the color for plotting")
    on_color = st.text_input("Ticked are ", value="C1")
    off_color = st.text_input("Unticked are: ", value="C0")
    colors = []
    labels = []
    group_on = []
    if len(file_list) == 0:
        st.write("Ending here")
        exit(0)
    st.write("### Tick the files ON")
    for i, afile in enumerate(file_list):
        thecheck = st.checkbox(str(afile.name), value=False,
                    key=f"grouping_{i}")
        if thecheck:
            colors.append(on_color)
            labels.append(on_label)
            group_on.append(thecheck)
        else:
            colors.append(off_color)
            labels.append(off_label)
            group_on.append(thecheck)
    group_on = np.array(group_on)
    group_labels = np.where(group_on, on_label, off_label)
    # colors = ["C0",
    n_on = np.count_nonzero(group_on)
    n_off = np.count_nonzero(np.logical_not(group_on))
    
    # mean_off = np.concatenate((1/n_off*np.ones(n_off), np.zeros(n_off)))
    mean_off = 1/n_off * np.logical_not(group_on)
    # mean_on = np.concatenate((np.zeros(n_on), 1/n_on*np.ones(n_on),))
    mean_on = 1/n_on * group_on
    mean_matrix = np.concatenate((mean_off[None,:], mean_on[None,:]))
    # mean_matrx: o e (out, experiment)
    # all_pol_ps: e b f (frequency, baseline, experiment)
    
    st.header("Computation preferences")
    sampling_period = st.number_input("Sampling period [s]", value=1.1e-3, step=1., format="%.2e")
    end_discard = st.number_input("Samples to discard at end of file", value=0, step=10, format="%d")
    nps = st.number_input(label="nperseg", min_value=10, step=1, value=3000)
    target_freq = st.number_input(label="Target frequency", step=0.1, value=200.)
    st.header("Plotting preferences")
    fstart = st.number_input(label="Start frequency", step=0.1, value=80.)
    fend = st.number_input(label="End frequency", step=0.1, value=300.)
    # 1.0e-2, 1e1
    ylims = st.number_input(label="y min", value=1.0e-2, format="%.1e"),\
            st.number_input(label="y max", value=2.e1, format="%.1e")


    all_pol_fs = []
    all_pol_ps = []
    all_headers = []

    if file_list is []:
        exit(0)
    for i, afile in enumerate(file_list):
        #_, ft_state, vib_state = filename.split('/')[-1].split('.')[0].split('_')[:3]
    
        hdul = fits.open(afile)
        opdc = hdul["OPDC"].data
        # opdc = fits.getdata(afile, 'OPDC')
        all_headers.append(hdul[0].header)
        hdul.close()
        t = opdc['TIME'][:-end_discard]
        opd = opdc['OPD'][:-end_discard]
        kopd = opdc['KALMAN_OPD'][:-end_discard]
        kpiezo = opdc['KALMAN_PIEZO'][:-end_discard]
        vlti = opdc['VLTI_DL_OFFSET'][:-end_discard]

        unwrapped_opd = 2.2 / (2*np.pi) * ((opd - kopd + np.pi) % (2*np.pi) - np.pi + kopd)
        pol = unwrapped_opd + (ptc.T2B @ (kpiezo - vlti*2*np.pi/2.2e-6).T).T
        one_pol_vector = []
        for apol in pol.T:
            af, aps = sig.welch(apol, fs=1/sampling_period, nperseg=nps, detrend='linear')
            one_pol_vector.append(aps)
        one_pol_vector = np.array(one_pol_vector)
        all_pol_fs.append(af)
        all_pol_ps.append(one_pol_vector)
    all_pol_ps = np.array(all_pol_ps).T
    all_pol_fs = np.array(all_pol_fs).T

    st.write(f"{on_label} = {n_on}, {off_label} = {n_off}")
    mean_ps = np.einsum("o e , f b e -> f b o",mean_matrix, all_pol_ps)

    normalization_choice = st.selectbox("Normalization mode", options=["mean","sample"], index=0)

    # Showing the timeline
    st.header("Timeline")
    check_plot_times = st.checkbox("Plot the time of acquisitions?")
    if check_plot_times:
        from datetime import datetime
        fig_times, time_list = plot_on_off_times(file_list, headers=all_headers,
                                colors=colors, group_labels=group_labels)
        # show_and_save(fig_times, "Timeline")
        show_and_save(fig_times, "timeline")

    offdata = all_pol_ps[:,:,mean_off>0]
    ondata = all_pol_ps[:,:,mean_on>0]
    if normalization_choice == "mean":
        mean_std = np.array([get_std(np.mean(ondata[:,i], axis=-1), np.std(ondata[:,i], axis=-1)/np.sqrt(ondata.shape[-1]),
                                 np.mean(offdata[:,i], axis=-1), np.std(offdata[:,i], axis=-1)/np.sqrt(offdata.shape[-1])) for i in range(6)]).T
    elif normalization_choice == "sample":
        mean_std = np.array([get_std(np.mean(ondata[:,i], axis=-1), np.std(ondata[:,i], axis=-1),
                                     np.mean(offdata[:,i], axis=-1), np.std(offdata[:,i], axis=-1)) for i in range(6)]).T

with tabs[1]:
    st.header("Colors")
    if check_plot_times:
        st.pyplot(fig_times)
    ylims_color = st.number_input(label="y min", value=1.0e-5, format="%.1e"),\
            st.number_input(label="y max", value=1.e0, format="%.1e")
    thealpha = st.number_input(label="alpha", value=0.5, step=0.05)
    thelinewidth = st.number_input(label="line width", value=1., step=0.05)
    plot_rcumsum = st.checkbox("Plot reverse cumulative", value=True,)

    with st.sidebar: # st.expander("Plot options_1"):
        cols_f = st.columns(2)
        flims = []
        with cols_f[0]:
            flims.append(st.number_input("f start [Hz]", value=40., step=10.))
            y_mode_log = st.checkbox("Log scale", value=True, key="man_pos_log" )
        with cols_f[1]:
            flims.append(st.number_input("F end [Hz]", value=300., step=10.))
        freq_cols = st.columns(4)
        ref_freqs = []
        for i, acol in enumerate(freq_cols):
            # with acol:
            ref_freqs.append(st.number_input("Frequency input", min_value=0., key=f"frequency_ref_1_{i}"))
            
    
    
    fig_color, axarr = plt.subplots(6, 1, sharex=True, sharey=True, figsize=(10, 16), dpi=200)

    for i, afile in enumerate(file_list):
        #_, ft_state, vib_state = filename.split('/')[-1].split('.')[0].split('_')[:3]
    
        for iBase in range(6):
            axarr[iBase].plot(all_pol_fs[:,0], all_pol_ps[:,iBase,i], color=colors[i],
                                alpha=thealpha, linewidth=thelinewidth,
                                label=afile.name)
            for afreq in ref_freqs:
                axarr[iBase].axvline(afreq, linewidth=0.5, color="k")
            # axarr[iBase].plot(x_pol, y_pol, label=f'{filepath.name}',color=colors[i], alpha=0.5)
            if plot_rcumsum:
                df = np.mean(np.gradient(all_pol_fs[:,0]))
                thercum = np.sqrt(np.cumsum(all_pol_ps[::-1, iBase, i]*df, axis=0))[::-1]
                axarr[iBase].plot(all_pol_fs[:,0], thercum, color=colors[i], alpha=0.5)
            plt.xscale("log")
            if y_mode_log:
                plt.yscale("log")
            #axarr[iBase].axvline(target_freq-1, color="k", linewidth=0.5)
            #axarr[iBase].axvline(target_freq+1, color="k", linewidth=0.5)
    for iBase in range(6):
        if plot_rcumsum:
            axarr[iBase].set_ylabel(f"{ptc.base2name[iBase]} [µm²/Hz]\nRev. cum. [µm]")
        else:
            axarr[iBase].set_ylabel(f"{ptc.base2name[iBase]} [µm²/Hz]")
    axarr[iBase].legend(loc="upper left", fontsize="xx-small")

    axarr[0].set_xlim(*flims)
    axarr[0].set_ylim(*ylims_color)
    show_and_save(fig_color, "fig_colors")



with tabs[2]:
    # with st.expander("Plot options"):
    #     cols_f = st.columns(2)
    #     flims = []
    #     with cols_f[0]:
    #         flims.append(st.number_input("f start [Hz]", value=1., step=10.,key="flim_1"))
    #         y_mode_log = st.checkbox("Log scale", value=False, key="man_pos_log_2" )
    #     with cols_f[1]:
    #         flims.append(st.number_input("F end [Hz]", value=300., step=10., key="flims_2"))
    #     freq_cols = st.columns(4)
    #     ref_freqs = []
    #     for i, acol in enumerate(freq_cols):
    #         with acol:
    #             ref_freqs.append(st.number_input("Frequency input", min_value=0., key=f"frequency_ref_2_{i}"))
            
    fig_comp_1, axarr = plt.subplots(6, 1, sharex=True, sharey=True, figsize=(10, 14), dpi=200)
    for bl_index in range(6):
        plt.sca(axarr[bl_index])
        plt.plot(all_pol_fs[:,0], mean_ps[:,bl_index,1]/mean_ps[:,bl_index,0], color="k")
        #plt.plot(all_pol_fs[:,0], mean_std[:,0, bl_index])
        #plt.plot(all_pol_fs[:,0], mean_std[:,1, bl_index])
        #plt.plot(all_pol_fs[:,0], mean_std[:,2, bl_index])
        plt.fill_between(all_pol_fs[:,0], mean_std[:,1, bl_index], mean_std[:,2, bl_index],
                         color="k", alpha=0.5, linewidth=0.)
        # plt.plt(all_pol_fs[:,0], )
        #plt.plot(all_pol_fs[:,0], mean_ps[:,bl_index,1])
        plt.axhline(1., color="k", linewidth=0.5, linestyle=":")
        plt.axvline(75.0, color="k", linewidth=0.2)
        for afreq in ref_freqs:
            plt.axvline(afreq, linewidth=0.5, color="k")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlim(*flims)
        plt.ylim(*ylims)
        plt.ylabel(f"Amplification {ptc.base2name[bl_index]}")
    plt.xlabel("Frequency [Hz]")
    # plt.savefig("plots/on_on_improvement_TTR-111-0009.pdf", dpi=200, bbox_inches="tight")

    show_and_save(fig_comp_1, name=tab_names[0])

with tabs[3]:
    ntraces = int(1e5/nps)
    fig_comp_2, axarr = plt.subplots(6, 1, sharex=True, sharey=True, figsize=(10, 14), dpi=200)
    for bl_index in range(6):
        plt.sca(axarr[bl_index])
        all_traces = get_mc(np.mean(ondata[:,bl_index], axis=-1), np.std(ondata[:,bl_index], axis=-1)/np.sqrt(ondata.shape[-1]),
                             np.mean(offdata[:,bl_index], axis=-1), np.mean(offdata[:,bl_index], axis=-1)/np.sqrt(offdata.shape[-1]),
                           n_traces=ntraces)
        plt.plot(all_pol_fs[:,0], mean_ps[:,bl_index,1]/mean_ps[:,bl_index,0], color="k")
        for atrace in all_traces:
            plt.plot(all_pol_fs[:,0], atrace, color="k", alpha=0.05,
                    linewidth=2.5)
        #plt.plot(all_pol_fs[:,0], mean_std[:,0, bl_index])
        #plt.plot(all_pol_fs[:,0], mean_std[:,1, bl_index])
        #plt.plot(all_pol_fs[:,0], mean_std[:,2, bl_index])
        # plt.plt(all_pol_fs[:,0], )
        #plt.plot(all_pol_fs[:,0], mean_ps[:,bl_index,1])
        plt.axhline(1., color="k", linewidth=0.5, linestyle=":")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlim(fstart, fend)
        plt.ylim(*ylims)
        plt.ylabel(f"Amplification {ptc.base2name[bl_index]}")
    plt.xlabel("Frequency [Hz]")
    # plt.savefig("plots/improvement_TTR-111-0009.png",
    #             dpi=200, bbox_inches="tight")

    show_and_save(fig_comp_2, name=tab_names[1])


