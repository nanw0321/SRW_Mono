import numpy as np
import os, copy, h5py

''' bandwidth '''
def get_spectra_from_file(fname):
    with h5py.File(fname,'r') as f:
        aw = f["spectrum/index"][:]
        axis_ev = f["spectrum/axis"][:]
        int_ev = f["spectrum/spectra"][:]
    return aw, axis_ev, int_ev

def calc_bandwidth(aw, axis_ev):
	return np.abs(axis_ev[aw.max()] - axis_ev[aw.min()])

''' throughput '''
def get_power_from_file(fname):
    axis, lineout, ori = get_lineout_from_file(fname)
    daxis = axis[1]-axis[0]
    power = lineout.sum()*daxis
    return power

def get_throughput_from_file(fname_in, fname_OE):
    power_in = get_power_from_file(fname_in)
    power_oe = get_power_from_file(fname_OE)
    return power_oe/power_in