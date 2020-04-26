import numpy as np

# OE thickness functions
def Perfect_lens(x,y,n,f,xoff=0,yoff=0):
    delta_d = f - np.sqrt(f**2-np.square(x-xoff)-np.square(y-yoff))
    z = -delta_d/(n-1)
    return z

def Perfect_lens_array(x,y,n,f,nlx,nly):
    nx,ny = x.shape
    z = np.zeros((nx,ny))
    npts = int(nx/nlx)
    dl = (x.max()-x.min())/(nx-1) * npts    # length of each micro lens
    if npts != nx/nlx:
        print('need integer multiples of small lens dimension')
    else:
        for il in range(nlx):
            xi = il * npts
            xf = (il+1)*npts
            for jl in range(nly):
                yi = jl * npts
                yf = (jl+1)*npts
                xl = x[xi:xf,yi:yf]
                yl = y[xi:xf,yi:yf]
                xlcent = xl.mean()
                ylcent = yl.mean()
                z[xi:xf,yi:yf] = Perfect_lens(xl,yl,n,f,xlcent,ylcent)
    return z

def Pinhole(x,y,r_hole,xoff=0.,yoff=0.):
    z = np.ones(x.shape)*1e30
    window = (np.square(x-xoff)+np.square(y-yoff))<=r_hole**2
    z[window]=0
    return z

def Pinhole_array(x,y,r_hole,nhx,nhy):
    nx,ny = x.shape
    z = np.ones((nx,ny))*1e30
    npts = int(nx/nhx)
    dl = (x.max()-x.min())/(nx-1) * npts
    if npts != nx/nhx:
        print('need integer multiples of small lens dimension')
    else:
        for il in range(nhx):
            xi = il * npts
            xf = (il+1)*npts
            for jl in range(nhy):
                yi = jl * npts
                yf = (jl+1)*npts
                xl = x[xi:xf,yi:yf]
                yl = y[xi:xf,yi:yf]
                xlcent = xl.mean()
                ylcent = yl.mean()
                z[xi:xf,yi:yf] = Pinhole(xl,yl,r_hole,xlcent,ylcent)
    return z

# function to calculate optical path differnece and amplitude transmission
def Calc_OPD_and_AmpTr(srwTr, thicknessProfData, n, d_abs):
    N = len(thicknessProfData[0])
    auxMesh = srwTr.mesh
    nx = auxMesh.nx
    ny = auxMesh.ny
    if thicknessProfData.shape==(nx,ny):
        # amplitude transmission
        Tr = np.exp(-thicknessProfData/d_abs)
        # OPD
        OPD = thicknessProfData*(n-1)
        if thicknessProfData.max() >= 1e30:
            # pinhole
            Tr[thicknessProfData>=1e30] = 0.
            OPD = np.zeros(OPD.shape)
        #for iy in range(ny):
        #    for ix in range(nx):
        for iy in (np.arange(int(ny*0.4))+int(ny*0.3)):
            for ix in (np.arange(int(nx*0.4))+int(nx*0.3)):
                ofst = 2*ix + 2*nx*iy
                srwTr.arTr[ofst]=Tr[iy,ix]
                srwTr.arTr[ofst+1]=OPD[iy,ix]
    else:
        print('OE shape not matched')

