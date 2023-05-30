
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy.signal as sig
import control

import graphviz

z_control = None

s, z = sp.symbols("s, z")

import graphviz
import re

    

class tac_obj(object):
    """
    A TAC application object used to graph the application from a source code
    """
    def __init__(self, filepath=None,src=None,
                build=True):
        """
        A TAC application object used to graph the application from a source code
        
        There is an additional interface to add instructions for this in the source code comments:
        In a comment after a #, add a '//' then you can add attrs for blocks or linkgs, to the exception of:
        * `name`   : already used for block identifier
        * `shape`  : already provided depending on nature of block
        * `comment`: already grabbed from # comments
        
        
        **Arguments: **
        * filepath : (None) Path to a TAC sourcecode `.tacsrc`
        * src      : Source code in a list of string (or a multiline string)
        * build    : If False: does not build.
        
        **After creating: **
        * `self.graph` in IPython.display: directly show the graph
        * `self.graph.view()` Open the .pdf file.
        """
        if src is not None:
            if not isinstance(src, list):
                ansrc = "\n".split(src)
            else:
                ansrc = src
            self.src = src 
        elif filepath is not None:
            with open(filepath, "r") as afile:
                self.src = afile.readlines()
        else:
            raise ArgumentError("Need to provide either")
        self.blocks = []
        
        if build:
            self.parse_source()
            self.build_graph()
        
    def parse_line(self, aline):
        splitted = aline.split("#")
        code_line = splitted[0]
        if len(splitted) >= 2 : # if there are xomments
            comment = splitted[1:]
            if len(splitted[-1].split("//")) >= 2:
                comment_last = splitted[-1].split("//")[-1]
            else :
                comment_last = None
                kwargs = {}
        else:
            comment_last = None
            comment = ""
            kwargs = {}
        
        if comment_last is not None:
            lskwargs = comment_last.replace(" ", "").split(",")
            kwargs = {}
            for akwarg in lskwargs:
                elements = akwarg.split("=")
                kwargs[elements[0]] = elements[1]
            kwargs["comment"] = ("".join(comment)).split("//")[0]

        if code_line is not None:
            args = re.split(r"\s+", code_line)
            
        else:
            return None
        
        return args, kwargs
    def parse_source(self):
        """
        Parse the source code of the TAC application
        """
        
        self.lines = []
        self.links = []
        for aline in self.src:
            if aline == "":
                pass
            elif aline[0] == "#":
                pass
            else:
                self.lines.append(self.parse_line(aline))
                
        for aline in self.lines:
            if aline[0][0] == "TAC_BLOCK":
                block_line = {"type":aline[0][0],
                            "name":aline[0][1],
                            "nature":aline[0][2],
                            "args":aline[0][3:],
                            "attrs":aline[1]}
                self.blocks.append(block_line)
            elif aline[0][0] == "TAC_LINK":
                link_line = {"type":aline[0][0],
                            "name":aline[0][1],
                            "block_a":aline[0][2],
                            "out_a":aline[0][3],
                            "block_b":aline[0][4],
                            "in_b":aline[0][5],
                            "attrs":aline[1]}
                self.links.append(link_line)
    def build_graph(self):
        """
        Build the graph
        
        The graph is accessible as `self.graph`
        """
        self.graph = graphviz.Digraph("TAC Source")
        self.graph.graph_attr['rankdir'] = 'LR' 
        for ablock in self.blocks:
            if ablock["nature"] == "Sum":
                ashape = "circle"
            elif ablock["nature"] == "Gain":
                ashape = "cds"
            else:
                ashape = "rectangle"
            argstr = str(ablock["args"])[:20] + "..."
            self.graph.node(name=ablock["name"],
                            label=f"{ablock['name']}\n{ablock['nature']}\n{argstr}",
                            shape=ashape,
                           fontsize="8",
                           **ablock["attrs"])
        for alink in self.links:
            self.graph.edge(tail_name=alink["block_a"],
                           head_name=alink["block_b"],
                           label=alink["name"],
                           fontsize="7",
                           **alink["attrs"])

def PLL_summary(the_PLLs):
    """
    Output a table of all PLLs involved
    
    * the_PLLs : a list of lists of dictionaries holding the PLL parameters
    
    Returns:
    * A table of the parameters, markdown format
    """
    mytab = []
    mytab.append(f"| UT | Identifier | Frequency [Hz] | Fmin [Hz] | Fmax[Hz] |tau [s] | Amp measured | Phase Measured [rad] |\n")
    mytab.append("|----|------------|----------------|-----------|----------|--------|--------------|----------------------|\n")
    for i, aUT in enumerate(the_PLLs, start=1):
        for apll in aUT:
            mytab.append(f"|{i}\
            |{apll['name']}\
            |{apll['FINIT']:.1f}\
            |{apll['FMIN']:.1f}\
            |{apll['FMAX']:.1f}|\
            {apll['PDTAU']:.2f}\
            |{apll['gain']:.2f}\
            |{apll['phase']:.2f}|\n")
    mytab = "".join(mytab)
    return mytab

def stripsrc(src="", block_names=[],match_names=[],
            attr_string=" # Remove this //color=red, style=dashed",
            verbose=True):
    """
    **Arguments**:
                
        * src    : ('') 
        * block_names    : ([]) A list of block names. Will remove all concerned lines
        * match_names    : ([]) A list of strings. Will remove lines that contain that string
        * attr_string    : (' # Remove this //color=red,style=dashed') 
        * verbose    : (True) ')
    
    **Returns**:
    * 
    """
    if isinstance(src, list):
        inputas = "list"
        mysrc = src
    elif isinstance(src, str):
        inputas = "str"
        mysrc = src.splitlines()
    newlines = []
    untouched = 0
    modified = 0
    removed = 0
    for aline in mysrc:
        code_line = aline.split("#")[0]
        allnames = re.split(r"\s+", code_line)
        tripped = False
        for ablock in block_names:
            if ablock in allnames:
                tripped = True
        for amatch in match_names:
            if amatch in aline:
                tripped = True
        if tripped:
            if attr_string is not None:
                print("Dashing : ", aline)
                newlines.append(aline+attr_string)
                modified += 1
            else:
                print("Not adding : ", aline)
                removed += 1
                pass
        else:
            untouched += 1
            newlines.append(aline)
    print(f"Untouched = {untouched}\nModified = {modified}\nRemoved = {removed}")
    print(f"Total = {len(newlines)}")
    if inputas == "list":
        out = newlines
    elif inputas == "str":
        out = "\n".join(newlines)
    return out

def replacelines(src, added=[], removed=[], match_names=[],
                 attr_string=" # Remove this //color=red, style=dashed",
                add_vibid=False,
                showdiff=True,
                verbose=False):
    """
    replacelines
        **Arguments** :
        
        * src    : 
        * added    : ([]) 
        * removed    : ([]) 
        * attr_string    : ('# Remove this //color=red,style=dashed') 
        * showdiff    : (True) 
        **returns** :
        
    **Returns:**
    * The modified source code
    * The diff source code
    """
    print(f"Initial length {len(src.splitlines())}")
    diffsrc = stripsrc(src=src, block_names=removed,
                       match_names=match_names,
                      attr_string=attr_string,
                      verbose=verbose)
    newsrc = stripsrc(src=src, block_names=removed,
                      match_names=match_names,
                      attr_string=None,
                      verbose=verbose)
    
    if add_vibid is not False:
        if verbose:
            print(f"Adding {add_vibid}")
        newsrc = re.sub(r'(\nTAC_VERSION.*\n)', f"\\1\n{add_vibid}", newsrc)
    if isinstance(added, list):
        addlist = added
    elif isinstance(added, str):
        addlist = added.splitlines()
    for anewline in addlist:
        if verbose:
            print("Appending ",anewline)
        diffsrc = "\n".join((diffsrc,anewline))
        newsrc = "\n".join((newsrc,anewline))
    
    if showdiff:
        atac = tac_obj(src=diffsrc.splitlines())
        atac.graph.view()
    
    print(f"New length {len(newsrc.splitlines())}")
    print(f"Diff length {len(diffsrc.splitlines())}")
    return newsrc, diffsrc


def assemble_PLL(pll_params, basename=None, u=None,
             gain=None, phase=None,
             input_block="In", output_block="Out",
            attr_string=" # //color=navyblue",
            mode="string",
            input_filter=None,
            sim_input=False,
            inverse_output=True,
            enable=1):
    """
    make_PLL
    
    **Arguments** :
    
    * pll_params    : [dict]  
        - "FINIT"
        - "FMIN"
        - "FMAX"
        - "KP"
        - "KI"
        - "KLP"
        - "PDTAU"
        - "name"
        - "phase"
        - "sum_input" : The index of input to which the output
        - "input_block" The block to use as input (supersedes)
          is connected
    * basename : (None) name used fo the 
    * u    : (None) Complex number for output gain (replaces `ug` and `phase`)
    * ug    : (None) The amplitude of gain
    * phase    : (None) [rad]
    * input_block    : ('In') 
    * output_block    : ('Out') 
    * attr_string    : ('# Remove this //color=red,style=dashed') 
    * mode           : ("string") otherwise: will return a list of lines
    * inverse_output : (True) caters for unexpected sign of the PLL
    * enable         : (int 0 or 1) the value to the enable line
    
    **returns** : Either the  
    * """
    if "sum_input" in pll_params.keys():
        output_index = pll_params["sum_input"]
    else:
        output_index = 1
    if "input_block" in pll_params.keys():
        input_block = pll_params["input_block"]
    if basename is None:
        basename = pll_params["name"]
    if u is not None:
        theta = np.angle(u)
        ug = np.abs(u)
    else:
        # Caution: sign of ug is changed later by 
        # `inverse_output`
        if gain is not None:
            ug = gain
        else:
            ug = 1.
        if phase is not None:
            theta = phase
        else:
            theta = 0
    if inverse_output:
        ug = -1 * ug
    else:
        ug = ug
    # ua and ub constitute first row of rotaion matrix
    ua = np.cos(theta)
    ub = - np.sin(theta)
    
    # The PLL arArithmeticErrorguments
    argnames = ["FINIT", "FMIN", "FMAX", "KP", "KI", "KLP", "PDTAU"]
    list_args = [format(pll_params[aname], ".2e") for aname in argnames]
    param_string = ",".join(list_args)
    
    src = []
    src.append(f"#\n")
    src.append(f"# Line of filter for {basename} \n")
    src.append(f"#\n")
    
    print(len(src))
    if input_filter is not None:
        src.append(input_filter)
        print(len(src))
    src.append(f"TAC_BLOCK {basename} PLL {param_string} {attr_string}\n")
    
    src.append(f"TAC_BLOCK {basename}_En Constant {enable} {attr_string}\n")
    src.append(f"TAC_BLOCK {basename}_Reset Constant 0 {attr_string}\n")
    src.append(f"TAC_BLOCK ua_{basename} Gain {ua:.3e} {attr_string}\n")
    src.append(f"TAC_BLOCK ub_{basename} Gain {ub:.3e} {attr_string}\n")
    src.append(f"TAC_BLOCK ug_{basename} Gain {ug:.3e} {attr_string}\n")
    # sum_string = "\""
    src.append(f"TAC_BLOCK {basename}_Sum Sum \\\"++\\\" {attr_string}\n")
    
    src.append(f"\n")
    
    if input_block is not None: # If we have added input filter, then there is no need to add an input link
        src.append(f"TAC_LINK {basename}_Input    {input_block}    1    {basename}    1 {attr_string}\n" )
    src.append(f"TAC_LINK {basename}_En_Link    {basename}_En    1    {basename}    2 {attr_string}\n" )
    src.append(f"TAC_LINK {basename}_Reset_Link    {basename}_Reset    1 {basename}    3 {attr_string}\n" )
    src.append(f"TAC_LINK {basename}_cos    {basename}    3    ua_{basename}    1 {attr_string}\n")
    src.append(f"TAC_LINK {basename}_sin    {basename}    2    ub_{basename}    1 {attr_string}\n")
    src.append(f"TAC_LINK {basename}_ua_Link    ua_{basename}    1    {basename}_Sum    1 {attr_string}\n")
    src.append(f"TAC_LINK {basename}_ub_Link    ub_{basename}    1    {basename}_Sum    2 {attr_string}\n")
    src.append(f"TAC_LINK {basename}_ug_Link    {basename}_Sum    1    ug_{basename}    1 {attr_string}\n")
    src.append(f"TAC_LINK {basename}_Output    ug_{basename}    1    {output_block}    {output_index} {attr_string}\n")
    
    # Adding simulated input
    
    if mode == "string":
        src = "".join(src)
    else:
        pass
    return src


def get_order(atf):
    numorder = atf.num[0][0].shape[0]-1
    denorder = atf.den[0][0].shape[0]-1
    return(np.max(numorder, denorder))
def get_degree_sp(expr):
    if expr == 1:
        return 0
    else:
        return sp.Poly(sp.expand(expr)).degree()
    
    
def tf2tac_sos(atf, zc , verbose=False, tacsource=True):
    """
    Output params in the way they are read by TAC: series of 
    gains and coefficients of order 2 TFs with coefficients
    of 1/z polynomials.
    
    *Arguments:*
    
    * atf   : The transfer function control object
    * zc    : The z symbol for the TF (used for sampling period)
    
    """
    lz = sp.symbols("z")
    pi = sp.symbols("p:6")
    sp.Add()
    sos = sig.tf2sos(atf.num[0][0], atf.den[0][0])
    params = []
    blocs = []
    for abloc in sos:
        anum, aden = abloc[:3], abloc[3:]
        atf = control.tf(anum, aden, dt=zc.dt)
        #print(abloc)
        if verbose:
            display(atf)
            #display((anum*zc**-2)/(aden*zc**-2))
        blocs.append(atf)
        gain = anum[0] / aden[0]
        a2a1 = anum/anum[0]
        b2b1 = aden/aden[0]
        a2 = a2a1[-1]
        a1 = a2a1[-2]
        b2 = b2b1[-1]
        b1 = b2b1[-2]
        if tacsource:
            print(f"{gain:.10e} {a1:.10e} {a2:.10e} {b1:.10e} {b2:.10e}")
        else:
            print(f"gain = {gain:.10e}")
            print(f"a1 = {a1}, a2 = {a2}, ")
            print(f"b1 = {b1}, b2 = {b2}, ")
        tfs = gain * (a2 + a1*z + z**2)/(b2 + b1*z +z**2)
        #tfs2 = gain * (1 + a1*z**-1 + a2*z**-2)/(1 + b1*z**-1 + b2*z**-2)
        params.append([gain, a1, a2, b1, b2])
        if verbose: 
            print("tfs")
            display(tfs)
            #print("tfs2")
            #display(tfs2)
    params = np.array(params)
    return params

def params2source(params, basename="Filter_M0",
                 input_block="<in_block>",
                 output_block="<out_block>",
                 comments="Filter transfer function",
                 attr_string=" # //color=darkgreen",
                 printit=True):
    
    # Creating the list of block names
    block_names = [f"{basename}_TF{i}" for i, dump in enumerate(params)]
    # Adding the output and input blocks
    block_all = block_names.copy()
    block_all.insert(0,input_block)
    if output_block is not None:
        block_all.append(output_block)
    
    src = []
    src.append(f"#\n")
    src.append(f"# Line of filter for {basename}\n")
    src.append(f"#\n")
    for i, (ablockname, ablockparam) in enumerate(zip(block_names, params)):
        param_string = " ".join([format(f"{apar:.15e}") for apar in ablockparam])
        if attr_string is not None:
            param_string += attr_string
        src.append(f"TAC_BLOCK {ablockname} DigitalTF {param_string}\n")
    src.append(f"\n")
    
    for i, (block_a, block_b) in enumerate(zip(block_all[:-1], block_all[1:])):
        src.append(f"TAC_LINK {basename}_L{i}    {block_a} 1    {block_b}    1 {attr_string}\n")
    
    
    src.append(f"\n")
    src.append(f"# End of filter {basename}\n")
    src.append(f"#\n")
    src = "".join(src)
    if printit:
        print(src)
    return src
    
        
        

def read_params(params, tacsource=True):
    """
    Prints the parameters with symbols
    (and a given float format)
    """
    if len(params.shape)==2:
        for aparam in params:
            gain, a1, a2, b1, b2 = aparam
            if tacsource:
                print(f"{gain:.10e} {a1:.10e} {a2:.10e} {b1:.10e} {b2:.10e}")
            else:
                print(f"gain = {gain:.10e}")
                print(f"a1 = {a1:.10e}, a2 = {a2:.10e}, ")
                print(f"b1 = {b1:.10e}, b2 = {b2:.10e}, ")

def save_filter(filename, params):
    """
    Save the parameters into a csv file
    """
    np.savetxt(filename, params, header="gain, a1, a2, b1, b2",
           delimiter="  ,   ", fmt="%.10f")
                    
def tac2tf(tac_params, dt):
    """
    Translates TAC parameters into a control TF
    """
    if len(tac_params.shape)==2:
        tfs = []
        for aparam in tac_params:
            gain, a1, a2, b1, b2 = aparam
            tfnum = gain*np.array([1, a1, a2])
            tfden = np.array([1, b1, b2])
            atf = control.tf(tfnum, tfden, dt)
            tfs.append(atf)
        return tfs
def tf_list_combine(tf_list):
    """
    Combine the list of control TF into a single TF
    """
    combined_tf = 1
    for atf in tf_list:
        combined_tf = combined_tf * atf
    return combined_tf

def check_tac(tac_array, original_tf, zc=z_control):
    """
    Comparison (with bode plot) of a TAC parameter set
    with the transfer function it is meant to represent.
    """
    tflist = tac2tf(tac_array, zc.dt)
    tf_combined = tf_list_combine(tflist)
    plt.figure()
    orig_amp, orig_phases, omega = control.bode(original_tf,
                                Hz=True, label="Original")
    compiled_amp, compiled_phases, omega = control.bode(tf_combined,
                                Hz=True, label="Compiled tf")
    plt.legend()
    plt.show()
    print("Max error on amp:", np.max(np.abs(orig_amp - compiled_amp)[:-2]))
    print("Max error on phase:", np.max(np.abs(orig_phases - compiled_phases)[:-2]))
    print("All close (1e-6)")
    print(np.allclose(orig_amp[:-2], compiled_amp[:-2], rtol=1e-6))
    print(np.allclose(orig_phases[:-2], compiled_phases[:-2], rtol=1e-6))



def get_stable_approx_inverse(atf, verbose=False, z=z,
                             z_control=z_control,
                              regularize=True):
    """
    Produces an approximate inverse following the method in Maggio2020
    
    *Parameters:*
    * atf   : transfer function of the actuator.
    * verbose: if `True`, will plot some information
    * z     : The sympy symbol to use for manipulating the transfer function
      internally
    * regulariz: If `True`, will add a delay to try to make the function causal
    
    """
    B = control.tf(atf.num, 1, z_control.dt)
    A = control.tf(atf.den, 1, z_control.dt)
    plt.figure()
    poles, zeros = control.pzmap(atf, grid=True)
    if not verbose:
        plt.close()
    else :
        plt.show()
        print(f"zeros : {np.abs(zeros)}")
    isoutside = np.abs(zeros)> 1
    outsidezeros = zeros[isoutside]
    insidezeros = zeros[np.logical_not(isoutside)]
    if verbose:
        print(zeros)
    Asp = control2sp(A, z)
    Bsp = control2sp(B, z)
    Bpoly = sp.Poly(Bsp)
    Bgain = Bpoly.all_coeffs()[0]
    Again = sp.Poly(Asp).all_coeffs()[0]
    #Bgain = sp.factor(Bsp).args[0]
    #Again = sp.factor(Asp).args[0]
    Bplus = sp.Mul(*[(z-azero) for azero in insidezeros])
    Bminus = Bgain*sp.Mul(*[(z-azero) for azero in outsidezeros])
    degA = sp.Poly(Asp).degree()
    degB = sp.Poly(Bsp).degree()
    if regularize:
        d_condition = 1*((degA - degB)>=1 )
    else:
        d_condition = 0
    degBminus = sp.Poly(sp.expand(Bminus)).degree()
    degBplus = get_degree_sp(Bplus)
    #Asp/(z**d*)
    #print("poles:",B.zeros())

    Bstarminus = sp.expand(z**degBminus * Bminus.subs([(z, 1/z)]))

    Bstarminus_tf = sp2control(Bstarminus, T_symbol=z_control, dt=z_control.dt)
    Bminus_tf = sp2control(Bminus, T_symbol=z_control, dt=z_control.dt)
    if verbose:
        plt.figure()
        control.pzmap(Bminus_tf, grid=True)
        plt.show()
        plt.figure()
        control.pzmap(Bstarminus_tf, grid=True)
        plt.show()

    Hdagger_sp = Asp/sp.expand(z**d_condition * Bplus * Bstarminus)
    Hdagger = sp2control(Hdagger_sp, z_control, dt=z_control.dt)
    #control.pzmap(Hdagger, grid=True)
    if verbose:
        plt.figure()
        dump, dump, omega_real = control.bode(
            1/Hdagger,
                    wrap_phase=True, Hz=True,  omega_limits=(10e-1, fmax),
                    label=f"$1/H^{{\dagger}}$")
        dump, dump, omega_real = control.bode(
            atf,
                    wrap_phase=True, Hz=True,  omega_limits=(10e-1, fmax),
                    label=f"$H_{{approx}}$")
        plt.legend()
    return Hdagger


def get_gain(reference_tf, freq):
    """
    Returns the complex gain of the the reference TF at a given
    frequency 
    
    *Arguments:*
    * reference_tf : `control.tf` object to use as a reference
    * freq         : The frequency to look up [Hz]
    
    *Returns:* The *Complex amplitude* of the TF at this freq
    """
    plt.figure()
    mag_ref, phase_ref, freq_ref = control.bode_plot(reference_tf, Hz=True, dB=True,
                                      plot=True,
                                      omega_limits=(1., 10000),
                                     wrap_phase=False)
    
    plt.close()
    res = np.abs(freq_ref - freq*2*np.pi)
    target_mag = mag_ref[np.argmin(res)]
    target_phase = mag_ref[np.argmin(res)]
    return target_mag*np.exp(1j*target_phase)

def build_notchfilter(f0_notch, reference_tf, order, 
                      notch_width, dt,
                      width_mode="relative",
                      conpensate_skewness=0.9):
    """
    Returns the TF of a notch MANHATTAN (feedfoward) to attenuate a given
    frequency. The filter is based on a butterworth bandpass filter.
    
    *Arguments:*
    * f0_notch: The frequency to cut
    * reference_tf : `control.tf` object to use as a reference.
      This should contain the idealized manhattan filter (even if it is unstable)
    * order : The order of the notch to use (the TF will be order*2)
    * dt    : The sampling period of the system
    * width_mode : a string default: `relative`
        - `relative` The width is defined as relative
        - `absolute` The width is defined in Hz
    * conpensate_skewness : Adjust the central frequency to account for the
      skewness of the double integration. Typically near ~1.
        - 1.2 for order 1.
        - 0.8 for order 2.
    """
    ref_gain = np.abs(get_gain(reference_tf, f0_notch))
    if width_mode == "relative":
        band_edges = np.array([(1-notch_width)*f0_notch, (1+notch_width)*f0_notch])
        if conpensate_skewness is not None: # Compensating the skewness (empirical rule of thumb)
            band_edges = band_edges + conpensate_skewness * f0_notch*notch_width/2
    elif width_mode == "absolute":
        band_edges = np.array([f0_notch - notch_width/2, f0_notch + notch_width/2])
        if conpensate_skewness: # Compensating the skewness (empirical rule of thumb)
            band_edges = band_edges + conpensate_skewness * notch_width/2
    else:
        raise KeyError("Mode unknown")
    numden = sig.butter(order, band_edges, fs=1/dt, btype="bandpass")
    the_tf_notch = -ref_gain * control.tf(*numden, dt=dt, )
    return the_tf_notch

from scipy.interpolate import interp1d
def resample(x, y, master, be=True):
    """Just a macro for interp1d"""
    values = interp1d(x, y, fill_value=np.nan, bounds_error=be, )(master)
    return values

def get_TF(t, x, y, nps=None, axis=0, get_coh=False):
    """      
    Compute the TF from measurements
            **Arguments** :
            
            * t    : Timestamps (used to get sample freq)
            * x    : Input time series
            * y    : Output time series
            * nps    : (None) nperseg for scipy.signal functions
            * axis    : (0) axis for scip.signal functions
            
            **returns** : 
            * """
    if nps is None:
        nps = 1e3
    fs = 1/np.mean(np.gradient(t))
    f1, csd = sig.csd(x, y,
                      nperseg=nps, fs=fs, axis=axis )
    f2, psd = sig.welch(x, nperseg=nps, fs=fs,
                       axis=axis)
    TFsig = csd/psd
    print("Same sampling=", np.allclose(f1, f2))
    if not get_coh:
        return f1, TFsig
    else:
        f3, coh = sig.coherence(x, y, fs=fs, nperseg=nps)
        print("Same sampling=", np.allclose(f1, f3))
        return f1, TFsig, coh


#import sympy as sp
#s, z = sp.symbols("s, z")
def control2sp(atf, asymbol):
    alist = [a * asymbol**i  for i, a in enumerate(np.flip(atf.num[0][0]))]
    num = sp.Add(*alist)
    alist = [a * asymbol**i  for i, a in enumerate(np.flip(atf.den[0][0]))]
    den = sp.Add(*alist)
    return num/den
def sp2control(atf, T_symbol, dt=1):
    
    num, den = atf.subs([(T_symbol, dt)]).as_numer_denom()
    if num.as_poly() is not None:
        numfloat = np.array(num.as_poly().all_coeffs(), dtype=np.float64)
    else:
        numfloat = 1
    if den.as_poly() is not None:
        denfloat = np.array(den.as_poly().all_coeffs(), dtype=np.float64)
    else:
        denfloat = 1
    H = control.tf(numfloat, denfloat , dt=dt)
    return H

