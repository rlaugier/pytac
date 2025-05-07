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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('ddlIndex', type=int)
    args = parser.parse_args()

    env = f'lddl{(args.ddlIndex-1)%2+1}'
    ddlx = f'ddl{args.ddlIndex}'

    ccs.CcsInit()

    ddlServer = ddl.DdlServer(env)
    tacServer = tac.TacServer(env)

    # Send DDL to standby
    ddlServer.standby()

    # Reset piezo controller, stage controller, and metrology
    tacServer.modblck(f'{ddlx}_StageControl', [ddl.PID_OPEN], [0])
    tacServer.modblck(f'{ddlx}_PiezoControl', [ddl.PID_OPEN], [0])
    tacServer.reset(f'{ddlx}_StageControl')
    tacServer.reset(f'{ddlx}_PiezoControl')
    time.sleep(1)
    tacServer.modblck(f'{ddlx}_Metrology', [ddl.METROLOGY_RESET], [2])

    # Enable disturbance
    tacServer.modblck('OffsetDis', [ddl.DISTURBANCE_ENABLE, 'Noise_white_1mu.txt'])
    time.sleep(1)

    # Configure mixer for disturbance input only
    tacServer.modblck(f'{ddlx}_Mixer', [ddl.MIXER_DIST_FFORWARD])

    # Enable piezo controller
    tacServer.modblck(f'{ddlx}_PiezoControl', [ddl.PID_CLOSE], [0])

    # Record telemetry
    timestamp, filepath = tacServer.record(f'{ddlx}_Probe', 40000)

    # Disable and reset piezo controller
    tacServer.modblck(f'{ddlx}_PiezoControl', [ddl.PID_OPEN], [0])
    tacServer.reset(f'{ddlx}_PiezoControl')

    # Configure mixer for no input
    tacServer.modblck(f'{ddlx}_Mixer', [ddl.MIXER_NONE])

    # Disable disturbance
    tacServer.modblck('OffsetDis', [ddl.DISTURBANCE_DISABLE, 'Noise_3mu.txt'])

    # Move record to current directory
    destination = os.path.join(os.getcwd(), os.path.basename(filepath))
    shutil.move(filepath, destination)
    filepath = destination

    # Read data
    data = np.loadtxt(filepath, dtype=[('timestamp', 'M8[us]'), ('setpoint', float), ('fforward', float), ('metrology', float), ('error', float), ('stage', float), ('piezo', float)])
    fforward = data['fforward']
    metrology = np.roll(data['metrology'], -1)
    dt = np.round(np.mean(np.diff(data['timestamp']).astype(float)))/1e6

    # Compute transfer function
    freq, Pxy = signal.csd(fforward, metrology, fs=1/dt)
    freq, Pxx = signal.welch(fforward, fs=1/dt)
    freq, coh = signal.coherence(fforward, metrology, fs=1/dt)
    tf = Pxy/Pxx

    # Estimate model
    #H = model.TFEst(freq[coh>0.98], tf[coh>0.98], dt).estimate([1, 1, 1], [1, 1])
    #print(model.digitalTF(H))
    #_, tfEst = H.freqresp(freq*2*np.pi*dt)

    # Plot result
    fig, axarr = plt.subplots(3, 1, sharex=True)
    axarr[0].set_title(f'DDL{args.ddlIndex} Feed Forward transfer function')
    axarr[0].axvline(300, ls='--', color='k', lw=1)
    axarr[0].axhline(0.98, ls='--', color='k', lw=1)
    axarr[0].semilogx(freq, coh, lw=2)
    axarr[1].axvline(300, ls='--', color='k', lw=1)
    axarr[1].loglog(freq, np.abs(tf), lw=2, label='measurement')
    #axarr[1].loglog(freq, np.abs(tfEst), label='model')
    axarr[2].axvline(300, ls='--', color='k', lw=1)
    axarr[2].axhline(-90, ls='--', color='k', lw=1)
    axarr[2].semilogx(freq, np.rad2deg(np.unwrap(np.angle(tf))), lw=2, label='measurement')
    #axarr[2].semilogx(freq, np.rad2deg(np.unwrap(np.angle(tfEst))), label='model')
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
    fig.savefig(f"Fig-TransferFunctionFeedForward-DDL{args.ddlIndex}-{timestamp}.png")
    plt.show()
