#!/bin/env python-qt5.9

import argparse
import os
import shutil
import time
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

import ccs
import ddl
import tac
import model


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('ddlIndex', type=int)
    parser.add_argument('--axis', type=int, choices=[1, 2, 3])
    args = parser.parse_args()

    env = f'lddl{(args.ddlIndex-1)%2+1}'
    ddlx = f'ddl{args.ddlIndex}'

    ccs.CcsInit()

    ddlServer = ddl.DdlServer(env)
    tacServer = tac.TacServer(env)

    # Send DDL to standby
    ddlServer.standby()

    # Enable disturbance
    tacServer.modblck('OffsetDis', [ddl.DISTURBANCE_ENABLE, 'Noise_white_1mu.txt'])
    time.sleep(1)

    # Disable stage controller
    tacServer.modblck(f'{ddlx}_StageControl', [ddl.PID_OPEN], [0])

    # Set piezo controller to through
    tacServer.modblck(f'{ddlx}_PiezoControl', [ddl.PID_THROUGH], [0])

    # If requested, do not command some of the axes
    if args.axis:
        for axis in [1, 2, 3]:
            if axis != args.axis:
                tacServer.modblck(f'{ddlx}_Piezo{axis}_m2V', [0,5])

    # Configure mixer for disturbance input only
    tacServer.modblck(f'{ddlx}_Mixer', [ddl.MIXER_DIST_FBACK])

    # Record telemetry
    timestamp, filepath = tacServer.record(f'{ddlx}_Probe', 40000)

    # Configure mixer for no input
    tacServer.modblck(f'{ddlx}_Mixer', [ddl.MIXER_NONE])

    # Disable disturbance
    tacServer.modblck('OffsetDis', [ddl.DISTURBANCE_DISABLE, 'Noise_3mu.txt'])

    # Re-enable all axes
    if args.axis:
        for axis in [1, 2, 3]:
            if axis != args.axis:
                tacServer.modblck(f'{ddlx}_Piezo{axis}_m2V', [-166666.666,5])

    # Move record to current directory
    destination = os.path.join(os.getcwd(), os.path.basename(filepath))
    shutil.move(filepath, destination)
    filepath = destination

    # Read data
    data = np.loadtxt(filepath, dtype=[('timestamp', 'M8[us]'), ('setpoint', float), ('fforward', float), ('metrology', float), ('error', float), ('stage', float), ('piezo', float)])
    setpoint = data['setpoint']
    metrology = np.roll(data['metrology'], -1)
    dt = np.round(np.mean(np.diff(data['timestamp']).astype(float)))/1e6

    # Compute transfer function
    freq, Pxy = signal.csd(setpoint, metrology, fs=1/dt)
    freq, Pxx = signal.welch(setpoint, fs=1/dt)
    freq, coh = signal.coherence(setpoint, metrology, fs=1/dt)
    tf = Pxy/Pxx

    # Estimate model
    mask = (coh>0.98) & (freq<1000)
    H = model.TFEst(freq[mask], tf[mask], dt).estimate([1, 1, 1], [1, 1])
    print(model.digitalTF(H))
    _, tfEst = H.freqresp(freq*2*np.pi*dt)

    # Custom display if axis specific
    if args.axis:
        title = f'DDL{args.ddlIndex} AXIS{args.axis} Piezo transfer function'
        filename = f"Fig-TransferFunctionPiezo-DDL{args.ddlIndex}-AXIS{args.axis}-{timestamp}.png"
    else:
        title = f'DDL{args.ddlIndex} Piezo transfer function'
        filename = f"Fig-TransferFunctionPiezo-DDL{args.ddlIndex}-{timestamp}.png"
    
    # Plot result
    fig, axarr = plt.subplots(3, 1, sharex=True)
    axarr[0].set_title(title)
    axarr[0].axvline(300, ls='--', color='k', lw=1)
    axarr[0].axhline(0.98, ls='--', color='k', lw=1)
    axarr[0].semilogx(freq, coh, lw=2)
    axarr[1].axvline(300, ls='--', color='k', lw=1)
    axarr[1].loglog(freq, np.abs(tf), lw=2, label='measurement')
    axarr[1].loglog(freq, np.abs(tfEst), label='model')
    axarr[2].axvline(300, ls='--', color='k', lw=1)
    axarr[2].axhline(-90, ls='--', color='k', lw=1)
    axarr[2].semilogx(freq, np.rad2deg(np.angle(tf)), lw=2, label='measurement')
    axarr[2].semilogx(freq, np.rad2deg(np.angle(tfEst)), label='model')
    axarr[0].set_xlim(1, 2e3)
    axarr[0].set_ylim(0, 1)
    axarr[1].set_ylim(1e-4, 1e1)
    axarr[2].set_ylim(-180, +180)
    axarr[0].grid()
    axarr[1].grid()
    axarr[2].grid()
    axarr[0].set_ylabel('Coherence')
    axarr[1].set_ylabel('Amplitude')
    axarr[2].set_ylabel('Phase')
    axarr[1].legend(loc='lower left')
    fig.tight_layout()
    fig.savefig(filename)
    plt.show()
