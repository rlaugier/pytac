
import streamlit as st
from pytac import tac_obj
from io import BytesIO
from io import StringIO


example_source = """
# This is example source code
TAC_BLOCK Noise_Gain Gain 1e-7                     #//color=forestgreen
TAC_BLOCK TF1_Noise DigitalTF 9.8472955380e-02 2.0000041742e+00 1.0000041743e+00 -4.1903685209e-01 3.5560304945e-01  #//color=forestgreen
TAC_BLOCK TF2_Noise DigitalTF 1.0000000000e+00 5.6044466135e-07 -9.9999109108e-01 -1.1582810619e+00 1.5854546526e-01  #//color=forestgreen
TAC_BLOCK TF3_Noise DigitalTF 1.0000000000e+00 -2.0000047347e+00 1.0000047347e+00 -1.9996858539e+00 9.9968595265e-01   #//color=forestgreen

TAC_BLOCK Out_Switch ManualSwitch 2 2  #//color=forestgreen



TAC_LINK Raw_noise Noise 1 Noise_Gain 1  #//color=forestgreen
TAC_LINK ngained Noise_Gain 1 TF1_Noise 1  #//color=forestgreen
TAC_LINK nfilt1 TF1_Noise 1 TF2_Noise 1  #//color=forestgreen
TAC_LINK nfilt2 TF2_Noise 1 TF3_Noise 1  #//color=forestgreen


TAC_LINK  F_Noise_to_switch TF3_Noise 1 TST_Switch 4  #//color=forestgreen

# Now we feed the output signal through a switch where one
# can select to use noise instead

TAC_LINK  Opl_to_Switch  Opl_Sign_Gain  1  Out_Switch 1    #//color=forestgreen
TAC_LINK  Filtered_Noise TF3_Noise      1  Out_Switch 2    #//color=forestgreen
TAC_LINK  Opl_or_Noise   Out_Switch     1  UT_Vib 2     #//color=forestgreen



# Extra switch for 
TAC_BLOCK Out_Info_Switch ManualSwitch 2,2     #//color=forestgreen
TAC_LINK Vib_Opl Opl_tot 1 Out_Info_Switch 1     #//color=forestgreen
TAC_LINK Filtered_Noise_Info TF3_Noise 1 Out_Info_Switch 2     #//color=forestgreen
TAC_LINK  Opl_or_Noise_Info   Out_Info_Switch  1   UT_Raw2 1     #//color=forestgreen"""

st.header("TAC schematics")
with st.expander("Help"):
    st.write("Allows to display and save TAC schematics to pdf")
    st.write("This renderer supports formating. write your formating\
        commands at the end of a `TAC_BLOCK` or `TAC_LINK` command\
        with after escaping it with `#//`.")
    st.write("Attributes can be found in [the graphviz documentation](https://graphviz.org/doc/info/attrs.html).")
with st.expander("Options"):
    encoding = st.text_input("Encoding", value="utf-8")
myfile = st.file_uploader("Upload a TAC source code:", type=["tacsrc"], )
if myfile is None:
    st.write("Or paste some source code in the text field.")
    mysrc = st.text_area("Source code:",
                value=example_source)
    mylines = mysrc.splitlines()
else:
    mysrc = StringIO(myfile.getvalue().decode(encoding))
    mylines = mysrc.readlines()
with st.expander("Informations"):
    if st.checkbox("Display source", value=False):
        st.write(mylines)

    st.write(f"File length: {len(mylines)} lines")

    myobj = tac_obj(src=mylines)
    st.write(f"Parsed lines: {len(myobj.lines)}")
    st.write(f"Links: {len(myobj.links)}")

# proj_name = st.text_input("Project name", value="tac_graph")
if st.checkbox("Show (ok on small graphs)"):
    myobj.graph

if st.checkbox("Open (preferable for large graphs)"):
    myobj.graph.render(view=True)
if st.checkbox("Save"):
    with BytesIO() as buffer:
        myobj.graph.render(buffer)
        st.write(buffer.name)
        st.download_button(label="Download graph",
                                data=buffer)