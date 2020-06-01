import numpy as np
import os
import pylab as plt

# import base wavefront and beamline class
from wpg import Wavefront, Beamline

# import base OE wrappers
from wpg.optical_elements import Aperture, Drift, CRL, Empty, Use_PP

# Gaussian beam generator
from wpg.generators import build_gauss_wavefront

#import SRW core functions
from wpg.srwlib import srwl, srwl_opt_setup_CRL, SRWLOptD, SRWLOptA, SRWLOptC, SRWLOptT, SRWLOptCryst

#import some helpers functions
from wpg.wpg_uti_wf import propagate_wavefront, plot_t_wf, get_intensity_on_axis
from wpg.wpg_uti_oe import show_transmission

def E2L(e):
    hbar = 6.582119569e-16
    omega = e/hbar
    frequency = omega /2/np.pi
    wavelength = 3e8/frequency
    return wavelength

def calc_stretching(thetaB, ang_as, range_xy):
    L = range_xy/np.tan(thetaB-ang_as)
    delta = np.cos(thetaB+ang_as)/np.sin(thetaB-ang_as)*range_xy - L
    trange = np.abs(delta)/3e8
    return trange

def mkdir_p(path):
    """
    Create directory with subfolders (like Linux mkdir -p)

    :param path: Path to be created
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def plot_spatial(wf):
    srwl.SetRepresElecField(wf._srwl_wf, 't')
    [xmin, xmax, ymin, ymax] = wf.get_limits()
    img = wf.get_intensity().sum(axis=-1)
    plt.figure()
    plt.imshow(img,cmap='jet',
        extent = [xmin*1e6, xmax*1e6, ymin*1e6, ymax*1e6])
    plt.colorbar()
    plt.xlabel(r'x ($\mu$m)',fontsize=18)
    plt.ylabel(r'y ($\mu$m)',fontsize=18)

def plot_lineout(wf, label=None, fov=1e30, ori='V', if_log=1):
    srwl.SetRepresElecField(wf._srwl_wf, 't')
    mesh = wf.params.Mesh
    [xmin, xmax, ymin, ymax] = wf.get_limits()
    img = wf.get_intensity().sum(axis=-1)
    if ori == 'V':
        # vertical plane
        axis = np.linspace(ymin, ymax, mesh.ny) * 1e6
        lineout = img[:,int(mesh.nx/2)]
        aname = 'y'
    else:
        # horizontal plane
        axis = np.linspace(xmin, xmax, mesh.nx) * 1e6
        lineout = img[int(mesh.ny/2),:]
        aname = 'x'
    index = np.argwhere(np.abs(axis)<fov/2)
    axis = axis[index.min():index.max()]
    lineout = lineout[index.min():index.max()]
    if if_log == 1:
        plt.plot(axis, np.log(lineout), label=label)
        ylabel = 'log beam intensity (a.u.)'
    else:
        plt.plot(axis, lineout, label=label)
        ylabel = 'beam intensity (a.u.)'
    plt.legend(fontsize=18)
    plt.title('lineout', fontsize=18)
    plt.xlabel(aname+r' ($\mu$m}', fontsize=18)
    plt.ylabel(ylabel, fontsize=18)

def get_spectra(wf):
    # change the wavefront into frequency domain, where each slice in z represent a photon energy.
    srwl.SetRepresElecField(wf._srwl_wf, 'f')
    mesh = wf.params.Mesh
    dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)        # spatial sampling resolution
    dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0)   # intensity per slice
    #int0_00 = wf.get_intensity()[int(mesh.ny/2), int(mesh.nx/2),:]  # central intensity
    dSlice = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices-1)   # photon energy sampling resolution
    axis_ev = np.arange(mesh.nSlices) * dSlice + mesh.sliceMin  # photon energy axis

    # get meaningful slices (>1% maximum intensity)
    int0max = max(int0)                                 # maximum intensity slice
    aw = [a[0] for a in np.argwhere(int0>int0max*0.01)] # meaningful indices
    #int0_mean = int0[min(aw):max(aw)]                   # meaningful intensity slices
    #axis_mean = np.arange(min(aw),max(aw)) * dSlice + mesh.sliceMin # meaningful photon energy axis
    return aw, axis_ev, int0

def get_temporal(wf):
    # change the wavefront into time domain, where each slice in z represent a time
    srwl.SetRepresElecField(wf._srwl_wf, 't')
    mesh = wf.params.Mesh
    dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)        # spatial sampling resolution
    dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0)   # intensity per slice
    #int0_00 = wf.get_intensity()[int(mesh.ny/2), int(mesh.nx/2),:]  # central intensity
    dSlice = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices-1)   # time sampling resolution
    axis_t = np.arange(mesh.nSlices) * dSlice + mesh.sliceMin  # time axis

    # get meaningful slices (>1% maximum intensity)
    int0max = max(int0)                                 # maximum intensity slice
    aw = [a[0] for a in np.argwhere(int0>int0max*0.01)] # meaningful indices
    #int0_mean = int0[min(aw):max(aw)]                   # meaningful intensity slices
    #axis_t = np.arange(min(aw),max(aw)) * dSlice + mesh.sliceMin # meaningful time axis
    return aw, axis_t, int0

def plot_spectra(wf, label=None):
    aw, axis_ev, int0 = get_spectra(wf)
    plt.plot(axis_ev, int0/int0.max(), label=label)
    plt.ylim([-0.1,1.1])
    plt.legend(fontsize=18)
    plt.title('spectral energy, meaningful range: {}-{} eV'.format(
        round(axis_ev[aw[0]],2),round(axis_ev[aw[-1]],2)),fontsize=18)
    plt.xlabel('eV',fontsize=18)
    plt.ylabel('normalized spectral energy', fontsize=18)
    return aw, axis_ev, int0

def plot_temporal(wf, label=None, fov=1e30, pulse_duration = None):
    aw, axis_t, int0 = get_temporal(wf)
    index = np.argwhere(np.abs(axis_t)<fov/2)
    norm_t = int0/int0.max()

    if norm_t.min() > 0.01:
        print(str(round(norm_t.min(),3))+', not enough sampling!!!')
    plt.plot(axis_t[index.min():index.max()]*1e15,
        norm_t[index.min():index.max()], label=label)
    plt.ylim([-0.1,1.1])
    plt.legend(fontsize=18)
    plt.title('temporal structure ({} fs rms pulse)'.format(
        pulse_duration*1e15),fontsize=18)
    plt.xlabel('t (fs)',fontsize=18)
    plt.ylabel('normalized temporal energy', fontsize=18)
    return aw, axis_t, int0

def set_crystal_orient(cryst, e, ang_dif_pl):
    nvx, nvy, nvz = cryst.find_orient(
        e,ang_dif_pl)[0][2]             # outward normal vector to crystal surface
    tvx, tvy, _ = cryst.find_orient(
        e,ang_dif_pl)[0][0]      # central tangential vector
    cryst.set_orient(nvx,nvy,nvz,tvx,tvy)