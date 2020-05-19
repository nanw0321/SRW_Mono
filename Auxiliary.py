import numpy as np
import os
import pylab

def calculate_theta_fwhm_cdr(ekev,qnC):
    """
    Calculate angular divergence using formula from XFEL CDR2011

    :param ekev: Energy in keV
    :param qnC: e-bunch charge, [nC]
    :return: theta_fwhm [units?]
    """
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC))*1e-6/ekev**0.85
    return theta_fwhm

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

def integral_intensity(wf, dirname, threshold=0.01, bPlot=True):
    """
    plot the slice-to-slice integral intensity averaged over a meaningful range

    :param wf: wavefront structure
    :param threshold: defined the threshold for slices, integrated_slice_intensity_max*threshold
    :param bPlot: if True plot temporary structure or spectrum in the meaningful range
    :return: intensity averaged over 'meaningful' slices, i.e. above 1% threshold, mainly needed for processing spiky FEL source

    """
    J2eV = 6.24150934e18
    # total0=wf.get_intensity().sum();
    mesh = wf.params.Mesh
    dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0)  # I(slice_num)
    int0 = int0 * (dx * dy * 1.e6)  # wf amplitude units sqrt(W/mm^2)

    # Get center pixel numbers.
    center_nx = int(mesh.nx / 2)
    center_ny = int(mesh.ny / 2)
    int0_00 = wf.get_intensity()[center_ny, center_nx, :]
    int0max = max(int0)

    # Get meaningful slices.
    aw = [a[0] for a in np.argwhere(int0 > int0max * threshold)]
    int0_mean = int0[min(aw):max(aw)]  # meaningful range of pulse
    if bPlot:
        if mesh.nSlices > 1:
            dSlice = (mesh.sliceMax - mesh.sliceMin) / (mesh.nSlices - 1)
        else:
            dSlice = 0
        pylab.figure()
        pylab.plot(np.arange(mesh.nSlices) * dSlice + mesh.sliceMin, int0)
        pylab.plot(np.arange(min(aw), max(aw)) *
                   dSlice + mesh.sliceMin, int0_mean, 'ro')
        if(wf.params.wDomain == 'time'):
            pylab.title('Power')
            pylab.xlabel('s')
            pylab.ylabel('W')
        else:  # frequency domain
            pylab.title('Spectral Energy')
            pylab.xlabel('eV')
            pylab.ylabel('J/eV')
        pylab.savefig(dirname+'spectrum.png')
        pylab.figure()
        pylab.plot(np.arange(mesh.nSlices) *
                   dSlice + mesh.sliceMin, int0_00)
        pylab.plot(np.arange(min(aw), max(aw)) * dSlice +
                   mesh.sliceMin, int0_00[min(aw):max(aw)], 'ro')
        if(wf.params.wDomain == 'time'):
            pylab.title('On-Axis Power Density')
            pylab.xlabel('s')
            pylab.ylabel('W/mm^2')
        else:  # frequency domain
            pylab.title('On-Axis Spectral Fluence')
            pylab.xlabel('eV')
            pylab.ylabel('J/eV/mm^2')
        pylab.savefig(dirname+'spectral_fluence.png')
    averaged = int0_mean.sum() / len(int0_mean)
    print('number of meaningful slices:', len(int0_mean))
    if(wf.params.wDomain == 'time'):
        dt = (mesh.sliceMax - mesh.sliceMin) / (mesh.nSlices - 1)
        print('Pulse energy {:1.2g} J'.format(int0_mean.sum() * dt))
    return averaged
