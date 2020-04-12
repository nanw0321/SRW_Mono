#!/usr/bin/env python
import os
try:
    __IPYTHON__
    import sys
    del sys.argv[1:]
except:
    pass


import srwl_bl
import srwlib
import srwlpy
import srwl_uti_smp


def set_optics(v=None):
    el = []
    pp = []
    names = ['Crystal', 'Crystal_s2_on_2nd_crystal', 's2_on_2nd_crystal_Crystal2', 'Crystal2', 'Crystal2_s3_on_1st_mirror', 's3_on_1st_mirror_Elliptical_Cylinder', 'Elliptical_Cylinder', 'CRL_s4_on_slit', 's4_on_slit', 's4_on_slit_Aperture', 'Aperture', 'Aperture_s5_on_2nd_mirror', 's5_on_2nd_mirror', 'Elliptical_Cylinder2', 'Elliptical_Cylinder2_CRL', 'CRL_s6_on_3rd_crystal', 's6_on_3rd_crystal', 's6_on_3rd_crystal_Crystal3', 'Crystal3', 'Crystal3_s7_on_4th_crystal', 's7_on_4th_crystal', 'Crystal4', 'Crystal4_on_detector', 'on_detector']
    for el_name in names:
        if el_name == 'Crystal':
            # Crystal: crystal 140.1m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_Crystal_d_sp,
                _psi0r=v.op_Crystal_psi0r,
                _psi0i=v.op_Crystal_psi0i,
                _psi_hr=v.op_Crystal_psiHr,
                _psi_hi=v.op_Crystal_psiHi,
                _psi_hbr=v.op_Crystal_psiHBr,
                _psi_hbi=v.op_Crystal_psiHBi,
                _tc=v.op_Crystal_tc,
                _ang_as=v.op_Crystal_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_Crystal_nvx,
                _nvy=v.op_Crystal_nvy,
                _nvz=v.op_Crystal_nvz,
                _tvx=v.op_Crystal_tvx,
                _tvy=v.op_Crystal_tvy,
            )
            el.append(crystal)
            pp.append(v.op_Crystal_pp)

        elif el_name == 'Crystal_s2_on_2nd_crystal':
            # Crystal_s2_on_2nd_crystal: drift 140.1m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Crystal_s2_on_2nd_crystal_L,
            ))
            pp.append(v.op_Crystal_s2_on_2nd_crystal_pp)
        elif el_name == 's2_on_2nd_crystal_Crystal2':
            # s2_on_2nd_crystal_Crystal2: drift 140.11788537m
            el.append(srwlib.SRWLOptD(
                _L=v.op_s2_on_2nd_crystal_Crystal2_L,
            ))
            pp.append(v.op_s2_on_2nd_crystal_Crystal2_pp)
        elif el_name == 'Crystal2':
            # Crystal2: crystal 140.11788538m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_Crystal2_d_sp,
                _psi0r=v.op_Crystal2_psi0r,
                _psi0i=v.op_Crystal2_psi0i,
                _psi_hr=v.op_Crystal2_psiHr,
                _psi_hi=v.op_Crystal2_psiHi,
                _psi_hbr=v.op_Crystal2_psiHBr,
                _psi_hbi=v.op_Crystal2_psiHBi,
                _tc=v.op_Crystal2_tc,
                _ang_as=v.op_Crystal2_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_Crystal2_nvx,
                _nvy=v.op_Crystal2_nvy,
                _nvz=v.op_Crystal2_nvz,
                _tvx=v.op_Crystal2_tvx,
                _tvy=v.op_Crystal2_tvy,
            )
            el.append(crystal)
            pp.append(v.op_Crystal2_pp)

        elif el_name == 'Crystal2_s3_on_1st_mirror':
            # Crystal2_s3_on_1st_mirror: drift 140.11788538m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Crystal2_s3_on_1st_mirror_L,
            ))
            pp.append(v.op_Crystal2_s3_on_1st_mirror_pp)
        elif el_name == 's3_on_1st_mirror_Elliptical_Cylinder':
            # s3_on_1st_mirror_Elliptical_Cylinder: drift 140.13788537m
            el.append(srwlib.SRWLOptD(
                _L=v.op_s3_on_1st_mirror_Elliptical_Cylinder_L,
            ))
            pp.append(v.op_s3_on_1st_mirror_Elliptical_Cylinder_pp)
        elif el_name == 'Elliptical_Cylinder':
            # Elliptical_Cylinder: ellipsoidMirror 140.13788538m
            el.append(srwlib.SRWLOptMirEl(
                _p=v.op_Elliptical_Cylinder_p,
                _q=v.op_Elliptical_Cylinder_q,
                _ang_graz=v.op_Elliptical_Cylinder_ang,
                _size_tang=v.op_Elliptical_Cylinder_size_tang,
                _size_sag=v.op_Elliptical_Cylinder_size_sag,
                _nvx=v.op_Elliptical_Cylinder_nvx,
                _nvy=v.op_Elliptical_Cylinder_nvy,
                _nvz=v.op_Elliptical_Cylinder_nvz,
                _tvx=v.op_Elliptical_Cylinder_tvx,
                _tvy=v.op_Elliptical_Cylinder_tvy,
                _x=v.op_Elliptical_Cylinder_x,
                _y=v.op_Elliptical_Cylinder_y,
            ))
            pp.append(v.op_Elliptical_Cylinder_pp)

        elif el_name == 'CRL_s4_on_slit':
            # CRL_s4_on_slit: drift 140.13788538m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL_s4_on_slit_L,
            ))
            pp.append(v.op_CRL_s4_on_slit_pp)
        elif el_name == 's4_on_slit':
            # s4_on_slit: watch 140.437856138m
            pass
        elif el_name == 's4_on_slit_Aperture':
            # s4_on_slit_Aperture: drift 140.437856138m
            el.append(srwlib.SRWLOptD(
                _L=v.op_s4_on_slit_Aperture_L,
            ))
            pp.append(v.op_s4_on_slit_Aperture_pp)
        elif el_name == 'Aperture':
            # Aperture: aperture 140.437856138m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Aperture_shape,
                _ap_or_ob='a',
                _Dx=v.op_Aperture_Dx,
                _Dy=v.op_Aperture_Dy,
                _x=v.op_Aperture_x,
                _y=v.op_Aperture_y,
            ))
            pp.append(v.op_Aperture_pp)
        elif el_name == 'Aperture_s5_on_2nd_mirror':
            # Aperture_s5_on_2nd_mirror: drift 140.437856138m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Aperture_s5_on_2nd_mirror_L,
            ))
            pp.append(v.op_Aperture_s5_on_2nd_mirror_pp)
        elif el_name == 's5_on_2nd_mirror':
            # s5_on_2nd_mirror: watch 140.737826895m
            pass
        elif el_name == 'Elliptical_Cylinder2':
            # Elliptical_Cylinder2: ellipsoidMirror 140.737826895m
            el.append(srwlib.SRWLOptMirEl(
                _p=v.op_Elliptical_Cylinder2_p,
                _q=v.op_Elliptical_Cylinder2_q,
                _ang_graz=v.op_Elliptical_Cylinder2_ang,
                _size_tang=v.op_Elliptical_Cylinder2_size_tang,
                _size_sag=v.op_Elliptical_Cylinder2_size_sag,
                _nvx=v.op_Elliptical_Cylinder2_nvx,
                _nvy=v.op_Elliptical_Cylinder2_nvy,
                _nvz=v.op_Elliptical_Cylinder2_nvz,
                _tvx=v.op_Elliptical_Cylinder2_tvx,
                _tvy=v.op_Elliptical_Cylinder2_tvy,
                _x=v.op_Elliptical_Cylinder2_x,
                _y=v.op_Elliptical_Cylinder2_y,
            ))
            pp.append(v.op_Elliptical_Cylinder2_pp)

        elif el_name == 'Elliptical_Cylinder2_CRL':
            # Elliptical_Cylinder2_CRL: drift 140.737826895m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Elliptical_Cylinder2_CRL_L,
            ))
            pp.append(v.op_Elliptical_Cylinder2_CRL_pp)
        elif el_name == 'CRL_s6_on_3rd_crystal':
            # CRL_s6_on_3rd_crystal: drift 140.75m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL_s6_on_3rd_crystal_L,
            ))
            pp.append(v.op_CRL_s6_on_3rd_crystal_pp)
        elif el_name == 's6_on_3rd_crystal':
            # s6_on_3rd_crystal: watch 140.757826895m
            pass
        elif el_name == 's6_on_3rd_crystal_Crystal3':
            # s6_on_3rd_crystal_Crystal3: drift 140.757826895m
            el.append(srwlib.SRWLOptD(
                _L=v.op_s6_on_3rd_crystal_Crystal3_L,
            ))
            pp.append(v.op_s6_on_3rd_crystal_Crystal3_pp)
        elif el_name == 'Crystal3':
            # Crystal3: crystal 140.757826895m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_Crystal3_d_sp,
                _psi0r=v.op_Crystal3_psi0r,
                _psi0i=v.op_Crystal3_psi0i,
                _psi_hr=v.op_Crystal3_psiHr,
                _psi_hi=v.op_Crystal3_psiHi,
                _psi_hbr=v.op_Crystal3_psiHBr,
                _psi_hbi=v.op_Crystal3_psiHBi,
                _tc=v.op_Crystal3_tc,
                _ang_as=v.op_Crystal3_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_Crystal3_nvx,
                _nvy=v.op_Crystal3_nvy,
                _nvz=v.op_Crystal3_nvz,
                _tvx=v.op_Crystal3_tvx,
                _tvy=v.op_Crystal3_tvy,
            )
            el.append(crystal)
            pp.append(v.op_Crystal3_pp)

        elif el_name == 'Crystal3_s7_on_4th_crystal':
            # Crystal3_s7_on_4th_crystal: drift 140.757826895m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Crystal3_s7_on_4th_crystal_L,
            ))
            pp.append(v.op_Crystal3_s7_on_4th_crystal_pp)
        elif el_name == 's7_on_4th_crystal':
            # s7_on_4th_crystal: watch 140.775712275m
            pass
        elif el_name == 'Crystal4':
            # Crystal4: crystal 140.775712275m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_Crystal4_d_sp,
                _psi0r=v.op_Crystal4_psi0r,
                _psi0i=v.op_Crystal4_psi0i,
                _psi_hr=v.op_Crystal4_psiHr,
                _psi_hi=v.op_Crystal4_psiHi,
                _psi_hbr=v.op_Crystal4_psiHBr,
                _psi_hbi=v.op_Crystal4_psiHBi,
                _tc=v.op_Crystal4_tc,
                _ang_as=v.op_Crystal4_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_Crystal4_nvx,
                _nvy=v.op_Crystal4_nvy,
                _nvz=v.op_Crystal4_nvz,
                _tvx=v.op_Crystal4_tvx,
                _tvy=v.op_Crystal4_tvy,
            )
            el.append(crystal)
            pp.append(v.op_Crystal4_pp)

        elif el_name == 'Crystal4_on_detector':
            # Crystal4_on_detector: drift 140.775712275m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Crystal4_on_detector_L,
            ))
            pp.append(v.op_Crystal4_on_detector_pp)
        elif el_name == 'on_detector':
            # on_detector: watch 150.0m
            pass
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)


varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'Hasan_Mono', 'simulation name'],

#---Data Folder
    ['fdir', 's', '', 'folder (directory) name for reading-in input and saving output data files'],


    ['gbm_x', 'f', 0.0, 'average horizontal coordinates of waist [m]'],
    ['gbm_y', 'f', 0.0, 'average vertical coordinates of waist [m]'],
    ['gbm_z', 'f', 0.0, 'average longitudinal coordinate of waist [m]'],
    ['gbm_xp', 'f', 0.0, 'average horizontal angle at waist [rad]'],
    ['gbm_yp', 'f', 0.0, 'average verical angle at waist [rad]'],
    ['gbm_ave', 'f', 4401.0, 'average photon energy [eV]'],
    ['gbm_pen', 'f', 0.001, 'energy per pulse [J]'],
    ['gbm_rep', 'f', 1, 'rep. rate [Hz]'],
    ['gbm_pol', 'f', 1, 'polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left'],
    ['gbm_sx', 'f', 4e-05, 'rms beam size vs horizontal position [m] at waist (for intensity)'],
    ['gbm_sy', 'f', 4e-05, 'rms beam size vs vertical position [m] at waist (for intensity)'],
    ['gbm_st', 'f', 1e-16, 'rms pulse duration [s] (for intensity)'],
    ['gbm_mx', 'f', 0, 'transverse Gauss-Hermite mode order in horizontal direction'],
    ['gbm_my', 'f', 0, 'transverse Gauss-Hermite mode order in vertical direction'],
    ['gbm_ca', 's', 'c', 'treat _sigX, _sigY as sizes in [m] in coordinate representation (_presCA="c") or as angular divergences in [rad] in angular representation (_presCA="a")'],
    ['gbm_ft', 's', 't', 'treat _sigT as pulse duration in [s] in time domain/representation (_presFT="t") or as bandwidth in [eV] in frequency domain/representation (_presFT="f")'],

#---Calculation Types
    # Electron Trajectory
    ['tr', '', '', 'calculate electron trajectory', 'store_true'],
    ['tr_cti', 'f', 0.0, 'initial time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_ctf', 'f', 0.0, 'final time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_np', 'f', 10000, 'number of points for trajectory calculation'],
    ['tr_mag', 'i', 1, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
    ['tr_fn', 's', 'res_trj.dat', 'file name for saving calculated trajectory data'],
    ['tr_pl', 's', '', 'plot the resulting trajectiry in graph(s): ""- dont plot, otherwise the string should list the trajectory components to plot'],

    #Single-Electron Spectrum vs Photon Energy
    ['ss', '', '', 'calculate single-e spectrum vs photon energy', 'store_true'],
    ['ss_ei', 'f', 100.0, 'initial photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ef', 'f', 20000.0, 'final photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ne', 'i', 10000, 'number of points vs photon energy for single-e spectrum vs photon energy calculation'],
    ['ss_x', 'f', 0.0, 'horizontal position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_y', 'f', 0.0, 'vertical position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_meth', 'i', 1, 'method to use for single-e spectrum vs photon energy calculation: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['ss_prec', 'f', 0.01, 'relative precision for single-e spectrum vs photon energy calculation (nominal value is 0.01)'],
    ['ss_pol', 'i', 6, 'polarization component to extract after spectrum vs photon energy calculation: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['ss_mag', 'i', 1, 'magnetic field to be used for single-e spectrum vs photon energy calculation: 1- approximate, 2- accurate'],
    ['ss_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
    ['ss_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['ss_fn', 's', 'res_spec_se.dat', 'file name for saving calculated single-e spectrum vs photon energy'],
    ['ss_pl', 's', '', 'plot the resulting single-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

    #Multi-Electron Spectrum vs Photon Energy (taking into account e-beam emittance, energy spread and collection aperture size)
    ['sm', '', '', 'calculate multi-e spectrum vs photon energy', 'store_true'],
    ['sm_ei', 'f', 100.0, 'initial photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ef', 'f', 20000.0, 'final photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ne', 'i', 10000, 'number of points vs photon energy for multi-e spectrum vs photon energy calculation'],
    ['sm_x', 'f', 0.0, 'horizontal center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_rx', 'f', 0.001, 'range of horizontal position / horizontal aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_nx', 'i', 1, 'number of points vs horizontal position for multi-e spectrum vs photon energy calculation'],
    ['sm_y', 'f', 0.0, 'vertical center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ry', 'f', 0.001, 'range of vertical position / vertical aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ny', 'i', 1, 'number of points vs vertical position for multi-e spectrum vs photon energy calculation'],
    ['sm_mag', 'i', 1, 'magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate'],
    ['sm_hi', 'i', 1, 'initial UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_hf', 'i', 15, 'final UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_prl', 'f', 1.0, 'longitudinal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_pra', 'f', 1.0, 'azimuthal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_meth', 'i', -1, 'method to use for spectrum vs photon energy calculation in case of arbitrary input magnetic field: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler", -1- dont use this accurate integration method (rather use approximate if possible)'],
    ['sm_prec', 'f', 0.01, 'relative precision for spectrum vs photon energy calculation in case of arbitrary input magnetic field (nominal value is 0.01)'],
    ['sm_nm', 'i', 1, 'number of macro-electrons for calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_na', 'i', 5, 'number of macro-electrons to average on each node at parallel (MPI-based) calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons) for intermediate intensity at calculation of multi-electron spectrum in case of arbitrary input magnetic field'],
    ['sm_type', 'i', 1, 'calculate flux (=1) or flux per unit surface (=2)'],
    ['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['sm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['sm_fn', 's', 'res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
    ['sm_pl', 's', '', 'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],
    #to add options for the multi-e calculation from "accurate" magnetic field

    #Power Density Distribution vs horizontal and vertical position
    ['pw', '', '', 'calculate SR power density distribution', 'store_true'],
    ['pw_x', 'f', 0.0, 'central horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_rx', 'f', 0.015, 'range of horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_nx', 'i', 100, 'number of points vs horizontal position for calculation of power density distribution'],
    ['pw_y', 'f', 0.0, 'central vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ry', 'f', 0.015, 'range of vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ny', 'i', 100, 'number of points vs vertical position for calculation of power density distribution'],
    ['pw_pr', 'f', 1.0, 'precision factor for calculation of power density distribution'],
    ['pw_meth', 'i', 1, 'power density computation method (1- "near field", 2- "far field")'],
    ['pw_zst', 'f', 0., 'initial longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_zfi', 'f', 0., 'final longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
    ['pw_fn', 's', 'res_pow.dat', 'file name for saving calculated power density distribution'],
    ['pw_pl', 's', '', 'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #Single-Electron Intensity distribution vs horizontal and vertical position
    ['si', '', '', 'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position', 'store_true'],
    #Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],
    #Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    ['w_e', 'f', 4401.0, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', -1.0, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', 0.001, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0.0, 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', 0.001, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],
    ['w_smpf', 'f', 0, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 2, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],

    ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '', 'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    ['wm_nm', 'i', 1000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    ['wm_na', 'i', 5, 'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    ['wm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
    ['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],
    ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
    ['wm_x0', 'f', 0, 'horizontal center position for mutual intensity cut calculation'],
    ['wm_y0', 'f', 0, 'vertical center position for mutual intensity cut calculation'],
    ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
    ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['wm_am', 'i', 0, 'multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance'],
    ['wm_fni', 's', 'res_int_pr_me.dat', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],

    #to add options
    ['op_r', 'f', 140.1, 'longitudinal position of the first optical element [m]'],

    # Former appParam:
    ['rs_type', 's', 'g', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:
    # Crystal: crystal
    ['op_Crystal_hfn', 's', '', 'heightProfileFile'],
    ['op_Crystal_dim', 's', 'x', 'orientation'],
    ['op_Crystal_d_sp', 'f', 3.13557135638, 'dSpacing'],
    ['op_Crystal_psi0r', 'f', -5.11321939903e-05, 'psi0r'],
    ['op_Crystal_psi0i', 'f', 3.58887560736e-06, 'psi0i'],
    ['op_Crystal_psiHr', 'f', -2.71341902391e-05, 'psiHr'],
    ['op_Crystal_psiHi', 'f', 2.50565837652e-06, 'psiHi'],
    ['op_Crystal_psiHBr', 'f', -2.71341902391e-05, 'psiHBr'],
    ['op_Crystal_psiHBi', 'f', 2.50565837652e-06, 'psiHBi'],
    ['op_Crystal_tc', 'f', 0.01, 'crystalThickness'],
    ['op_Crystal_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_Crystal_nvx', 'f', 0.0, 'nvx'],
    ['op_Crystal_nvy', 'f', 0.893336362055, 'nvy'],
    ['op_Crystal_nvz', 'f', -0.449388633846, 'nvz'],
    ['op_Crystal_tvx', 'f', 0.0, 'tvx'],
    ['op_Crystal_tvy', 'f', 0.449388633846, 'tvy'],
    ['op_Crystal_ang', 'f', 0.466080858374, 'grazingAngle'],
    ['op_Crystal_amp_coef', 'f', 1.0, 'heightAmplification'],

    # Crystal_s2_on_2nd_crystal: drift
    ['op_Crystal_s2_on_2nd_crystal_L', 'f', 0.01788537, 'length'],

    # s2_on_2nd_crystal_Crystal2: drift
    ['op_s2_on_2nd_crystal_Crystal2_L', 'f', 9.99997951112e-09, 'length'],

    # Crystal2: crystal
    ['op_Crystal2_hfn', 's', '', 'heightProfileFile'],
    ['op_Crystal2_dim', 's', 'x', 'orientation'],
    ['op_Crystal2_d_sp', 'f', 3.13557135638, 'dSpacing'],
    ['op_Crystal2_psi0r', 'f', -5.11321939903e-05, 'psi0r'],
    ['op_Crystal2_psi0i', 'f', 3.58887560736e-06, 'psi0i'],
    ['op_Crystal2_psiHr', 'f', -2.71341902391e-05, 'psiHr'],
    ['op_Crystal2_psiHi', 'f', 2.50565837652e-06, 'psiHi'],
    ['op_Crystal2_psiHBr', 'f', -2.71341902391e-05, 'psiHBr'],
    ['op_Crystal2_psiHBi', 'f', 2.50565837652e-06, 'psiHBi'],
    ['op_Crystal2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_Crystal2_ang_as', 'f', -0.366519142919, 'asymmetryAngle'],
    ['op_Crystal2_nvx', 'f', -2.4158420418e-09, 'nvx'],
    ['op_Crystal2_nvy', 'f', -0.672975300169, 'nvy'],
    ['op_Crystal2_nvz', 'f', -0.739664954802, 'nvz'],
    ['op_Crystal2_tvx', 'f', -2.65524409917e-09, 'tvx'],
    ['op_Crystal2_tvy', 'f', -0.739664954802, 'tvy'],
    ['op_Crystal2_ang', 'f', 0.83257236543, 'grazingAngle'],
    ['op_Crystal2_amp_coef', 'f', 1.0, 'heightAmplification'],

    # Crystal2_s3_on_1st_mirror: drift
    ['op_Crystal2_s3_on_1st_mirror_L', 'f', 0.01999999, 'length'],

    # s3_on_1st_mirror_Elliptical_Cylinder: drift
    ['op_s3_on_1st_mirror_Elliptical_Cylinder_L', 'f', 1.00000079328e-08, 'length'],

    # Elliptical_Cylinder: ellipsoidMirror
    ['op_Elliptical_Cylinder_hfn', 's', '', 'heightProfileFile'],
    ['op_Elliptical_Cylinder_dim', 's', 'x', 'orientation'],
    ['op_Elliptical_Cylinder_p', 'f', 10.0, 'firstFocusLength'],
    ['op_Elliptical_Cylinder_q', 'f', 0.3, 'focalLength'],
    ['op_Elliptical_Cylinder_ang', 'f', 0.00698131700798, 'grazingAngle'],
    ['op_Elliptical_Cylinder_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_Elliptical_Cylinder_size_tang', 'f', 0.1, 'tangentialSize'],
    ['op_Elliptical_Cylinder_size_sag', 'f', 0.001, 'sagittalSize'],
    ['op_Elliptical_Cylinder_nvx', 'f', 0.0, 'normalVectorX'],
    ['op_Elliptical_Cylinder_nvy', 'f', 0.999975630705, 'normalVectorY'],
    ['op_Elliptical_Cylinder_nvz', 'f', -0.00698126029796, 'normalVectorZ'],
    ['op_Elliptical_Cylinder_tvx', 'f', 0.0, 'tangentialVectorX'],
    ['op_Elliptical_Cylinder_tvy', 'f', 0.00698126029796, 'tangentialVectorY'],
    ['op_Elliptical_Cylinder_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Elliptical_Cylinder_y', 'f', 0.0, 'verticalOffset'],

    # CRL_s4_on_slit: drift
    ['op_CRL_s4_on_slit_L', 'f', 0.29997075755, 'length'],

    # s4_on_slit_Aperture: drift
    ['op_s4_on_slit_Aperture_L', 'f', 1.00044417195e-11, 'length'],

    # Aperture: aperture
    ['op_Aperture_shape', 's', 'r', 'shape'],
    ['op_Aperture_Dx', 'f', 0.001, 'horizontalSize'],
    ['op_Aperture_Dy', 'f', 5e-05, 'verticalSize'],
    ['op_Aperture_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Aperture_y', 'f', 0.0, 'verticalOffset'],

    # Aperture_s5_on_2nd_mirror: drift
    ['op_Aperture_s5_on_2nd_mirror_L', 'f', 0.29997075719, 'length'],

    # Elliptical_Cylinder2: ellipsoidMirror
    ['op_Elliptical_Cylinder2_hfn', 's', '', 'heightProfileFile'],
    ['op_Elliptical_Cylinder2_dim', 's', 'x', 'orientation'],
    ['op_Elliptical_Cylinder2_p', 'f', 0.3, 'firstFocusLength'],
    ['op_Elliptical_Cylinder2_q', 'f', 10.0, 'focalLength'],
    ['op_Elliptical_Cylinder2_ang', 'f', 0.00698131700798, 'grazingAngle'],
    ['op_Elliptical_Cylinder2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_Elliptical_Cylinder2_size_tang', 'f', 0.1, 'tangentialSize'],
    ['op_Elliptical_Cylinder2_size_sag', 'f', 0.001, 'sagittalSize'],
    ['op_Elliptical_Cylinder2_nvx', 'f', 0.0, 'normalVectorX'],
    ['op_Elliptical_Cylinder2_nvy', 'f', -0.999975630705, 'normalVectorY'],
    ['op_Elliptical_Cylinder2_nvz', 'f', -0.00698126029796, 'normalVectorZ'],
    ['op_Elliptical_Cylinder2_tvx', 'f', 0.0, 'tangentialVectorX'],
    ['op_Elliptical_Cylinder2_tvy', 'f', 0.00698126029796, 'tangentialVectorY'],
    ['op_Elliptical_Cylinder2_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Elliptical_Cylinder2_y', 'f', 0.0, 'verticalOffset'],

    # Elliptical_Cylinder2_CRL: drift
    ['op_Elliptical_Cylinder2_CRL_L', 'f', 0.01217310525, 'length'],

    # CRL_s6_on_3rd_crystal: drift
    ['op_CRL_s6_on_3rd_crystal_L', 'f', 0.00782689475, 'length'],

    # s6_on_3rd_crystal_Crystal3: drift
    ['op_s6_on_3rd_crystal_Crystal3_L', 'f', 1.00044417195e-11, 'length'],

    # Crystal3: crystal
    ['op_Crystal3_hfn', 's', '', 'heightProfileFile'],
    ['op_Crystal3_dim', 's', 'x', 'orientation'],
    ['op_Crystal3_d_sp', 'f', 3.13557135638, 'dSpacing'],
    ['op_Crystal3_psi0r', 'f', -5.11321939903e-05, 'psi0r'],
    ['op_Crystal3_psi0i', 'f', 3.58887560736e-06, 'psi0i'],
    ['op_Crystal3_psiHr', 'f', -2.71341902391e-05, 'psiHr'],
    ['op_Crystal3_psiHi', 'f', 2.50565837652e-06, 'psiHi'],
    ['op_Crystal3_psiHBr', 'f', -2.71341902391e-05, 'psiHBr'],
    ['op_Crystal3_psiHBi', 'f', 2.50565837652e-06, 'psiHBi'],
    ['op_Crystal3_tc', 'f', 0.01, 'crystalThickness'],
    ['op_Crystal3_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_Crystal3_nvx', 'f', 0.0, 'nvx'],
    ['op_Crystal3_nvy', 'f', -0.893336362055, 'nvy'],
    ['op_Crystal3_nvz', 'f', -0.449388633846, 'nvz'],
    ['op_Crystal3_tvx', 'f', -1.61321218547e-09, 'tvx'],
    ['op_Crystal3_tvy', 'f', -0.449388633846, 'tvy'],
    ['op_Crystal3_ang', 'f', 0.466080858374, 'grazingAngle'],
    ['op_Crystal3_amp_coef', 'f', 1.0, 'heightAmplification'],

    # Crystal3_s7_on_4th_crystal: drift
    ['op_Crystal3_s7_on_4th_crystal_L', 'f', 0.01788538036, 'length'],

    # Crystal4: crystal
    ['op_Crystal4_hfn', 's', '', 'heightProfileFile'],
    ['op_Crystal4_dim', 's', 'x', 'orientation'],
    ['op_Crystal4_d_sp', 'f', 3.13557135638, 'dSpacing'],
    ['op_Crystal4_psi0r', 'f', -5.11321939903e-05, 'psi0r'],
    ['op_Crystal4_psi0i', 'f', 3.58887560736e-06, 'psi0i'],
    ['op_Crystal4_psiHr', 'f', -2.71341902391e-05, 'psiHr'],
    ['op_Crystal4_psiHi', 'f', 2.50565837652e-06, 'psiHi'],
    ['op_Crystal4_psiHBr', 'f', -2.71341902391e-05, 'psiHBr'],
    ['op_Crystal4_psiHBi', 'f', 2.50565837652e-06, 'psiHBi'],
    ['op_Crystal4_tc', 'f', 0.01, 'crystalThickness'],
    ['op_Crystal4_ang_as', 'f', 0.366519142919, 'asymmetryAngle'],
    ['op_Crystal4_nvx', 'f', 0.0, 'nvx'],
    ['op_Crystal4_nvy', 'f', 0.995027350408, 'nvy'],
    ['op_Crystal4_nvz', 'f', -0.0996020679548, 'nvz'],
    ['op_Crystal4_tvx', 'f', 0.0, 'tvx'],
    ['op_Crystal4_tvy', 'f', 0.0996020679548, 'tvy'],
    ['op_Crystal4_ang', 'f', 0.099767492435, 'grazingAngle'],
    ['op_Crystal4_amp_coef', 'f', 1.0, 'heightAmplification'],

    # Crystal4_on_detector: drift
    ['op_Crystal4_on_detector_L', 'f', 9.22428772488, 'length'],

#---Propagation parameters
    ['op_Crystal_pp', 'f',                              [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal'],
    ['op_Crystal_s2_on_2nd_crystal_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal_s2_on_2nd_crystal'],
    ['op_s2_on_2nd_crystal_Crystal2_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 's2_on_2nd_crystal_Crystal2'],
    ['op_Crystal2_pp', 'f',                             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal2'],
    ['op_Crystal2_s3_on_1st_mirror_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal2_s3_on_1st_mirror'],
    ['op_s3_on_1st_mirror_Elliptical_Cylinder_pp', 'f', [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 's3_on_1st_mirror_Elliptical_Cylinder'],
    ['op_Elliptical_Cylinder_pp', 'f',                  [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Elliptical_Cylinder'],
    ['op_CRL_s4_on_slit_pp', 'f',                       [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL_s4_on_slit'],
    ['op_s4_on_slit_Aperture_pp', 'f',                  [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 's4_on_slit_Aperture'],
    ['op_Aperture_pp', 'f',                             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture'],
    ['op_Aperture_s5_on_2nd_mirror_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture_s5_on_2nd_mirror'],
    ['op_Elliptical_Cylinder2_pp', 'f',                 [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Elliptical_Cylinder2'],
    ['op_Elliptical_Cylinder2_CRL_pp', 'f',             [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Elliptical_Cylinder2_CRL'],
    ['op_CRL_s6_on_3rd_crystal_pp', 'f',                [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL_s6_on_3rd_crystal'],
    ['op_s6_on_3rd_crystal_Crystal3_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 's6_on_3rd_crystal_Crystal3'],
    ['op_Crystal3_pp', 'f',                             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal3'],
    ['op_Crystal3_s7_on_4th_crystal_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal3_s7_on_4th_crystal'],
    ['op_Crystal4_pp', 'f',                             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal4'],
    ['op_Crystal4_on_detector_pp', 'f',                 [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Crystal4_on_detector'],
    ['op_fin_pp', 'f',                                  [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],

    #[ 0]: Auto-Resize (1) or not (0) Before propagation
    #[ 1]: Auto-Resize (1) or not (0) After propagation
    #[ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
    #[ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
    #[ 6]: Horizontal Resolution modification factor at Resizing
    #[ 7]: Vertical Range modification factor at Resizing
    #[ 8]: Vertical Resolution modification factor at Resizing
    #[ 9]: Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
    #[12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
    #[13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
    #[14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
    #[15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
    #[16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate
])


def main():
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    op = set_optics(v)
    v.si = True
    v.si_pl = 'xy'
    v.ws = True
    v.ws_pl = 'xy'
    mag = None
    if v.rs_type == 'm':
        mag = srwlib.SRWLMagFldC()
        mag.arXc.append(0)
        mag.arYc.append(0)
        mag.arMagFld.append(srwlib.SRWLMagFldM(v.mp_field, v.mp_order, v.mp_distribution, v.mp_len))
        mag.arZc.append(v.mp_zc)
    srwl_bl.SRWLBeamline(_name=v.name, _mag_approx=mag).calc_all(v, op)

main()
