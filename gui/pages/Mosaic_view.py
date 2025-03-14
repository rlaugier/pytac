import streamlit as st

import astropy.io.fits as fits
import numpy as np
from pytac.read_manhattan import get_mosaic
import io
from pytac_vlti_context import rev_search, mirror_indices 

st.title("Creating mosaic plots")
file_target = st.file_uploader("Upload a file for mosaic",
                accept_multiple_files=False, type=["fits"])
matrix_colors = None
matrix_faulty = np.zeros((12,4), dtype=bool)
st.header("Coloring of plots")
if st.checkbox(label="Orange", value=False):
    all_faulty = []
    cols = st.columns(4+1)
    for i, acol in enumerate(cols[1:], start=1):
        column_faulty = []
        with acol:
            st.write(f"## UT{i}")
            for j in range(1, 12+1):
                mirror_name = rev_search(mirror_indices, j-1)
                column_faulty.append(st.checkbox(label=f"Acc{j} M{mirror_name}", key=f"orange_{i}_{j}"))
        all_faulty.append(column_faulty)
    all_faulty = np.array(all_faulty).T
    matrix_faulty = np.where(all_faulty, 1, matrix_faulty)
    matrix_colors = matrix_faulty
if st.checkbox(label="Grey", value=False):
    all_faulty = []
    cols = st.columns(4+1)
    for i, acol in enumerate(cols[1:], start=1):
        column_faulty = []
        with acol:
            st.write(f"## UT{i}")
            for j in range(1, 12+1):
                mirror_name = rev_search(mirror_indices, j-1)
                column_faulty.append(st.checkbox(label=f"Acc{j} M{mirror_name}", key=f"grey_{i}_{j}"))
        all_faulty.append(column_faulty)
    all_faulty = np.array(all_faulty).T
    matrix_faulty = np.where(all_faulty, 2, matrix_colors)
    matrix_colors = matrix_faulty

myfs = st.number_input("Sampling frequency [Hz]", value=1000, step=10, format="%f")
    
st.write("### The colors requested:")
st.write(matrix_colors)
st.header("The resutling plots")
if file_target is not None:
    hdul = fits.open(file_target)
    figs = get_mosaic(hdul, showall=False,
                saveall=False, faulty_matrix=matrix_colors,
                fs=myfs)

    tab_names = ["Accelerations",
        "Acc. power spectrum",
        "Position power spectrum",
        "Revers cumulative position"]
    tabs = st.tabs(tab_names)
    for i, atab in enumerate(tabs):
        with atab:
            st.pyplot(figs[i])
            # if st.button("Save the plot", key=f"button_save{i}"):
            format = st.selectbox("Format to save", options=["pdf", "png"], key=f"format_{i}")
            if format is not None:
                with io.BytesIO() as buffer:
                    figs[i].savefig(buffer, bbox_inches="tight",format="pdf", dpi=200)
                    st.download_button("Download pdf", data=buffer,
                                file_name=f"mosaic_{tab_names[i]}.pdf",
                                key=f"button_download{i}")

