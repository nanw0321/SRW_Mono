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
        print('need integer multiples of small pinhole dimension')
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
'''
def Prism(x,y,slope,intersection,xoff=0.,yoff=0.):
    # input x, y: the x&y axis for the a single prism (meshgrid).
    # input xoff, yoff: x,y coordinate of the prism center (default to 0,0)
    nx,ny = x.shape
    window_size = (np.max(x) - np.min(x))*0.3/2
    # z1 is prism1 
    z1 = np.ones((nx,ny)) * 1e30
    left1 = np.where( (x[0,:]-xoff) <= - window_size)
    window1 = np.where(( (x[0,:]-xoff) < window_size) & (-window_size < (x[0,:]-xoff)))
    right1 = np.where(window_size <= (x[0,:]-xoff) )

    for i in range(ny):
        z1[left1,i] = - slope * (x[0,left1]-xoff) + intersection
        z1[right1,i] = slope * (x[0,right1]-xoff) + intersection
        z1[window1,i] = 0 

    #z2 is prism2
    z2 = np.ones((nx,ny)) * 1e30
    left2 = np.where((y[:,0]-yoff) <= - window_size)
    window2 = np.where(((y[:,0]-yoff) < window_size) & (-window_size < (y[:,0]-yoff)))
    right2 = np.where(window_size <= (y[:,0]-yoff))

    for i in range(nx):
        z2[i,left2] = - slope * (y[left2,0]-yoff) + intersection
        z2[i,right2] = slope * (y[right2,0]-yoff) + intersection
        z2[i,window2] = 0 
    z = z1 + z2
    return z
'''
def Prism(x,y,slope,d_hole,xoff=0.,yoff=0.):
    nx,ny = x.shape
    # For prism in x:
    z1 = np.zeros((nx,ny))          # initialize prism with infinite thickness
    index1 = (x[0]-xoff)>=d_hole/2. # one edge
    nn = index1.sum()
    z1[:,-nn:] = slope * (x[:,-nn:]-xoff)
    z1 = z1 + np.fliplr(z1)
    
    # For prism in y:
    z2 = z1.T
    z = z1+z2
    return z

def Prism_array(x,y,slope,d_hole,npx,npy):
    # input npx, npy: # of individual prisms in each dimension.
    # output z: height profile of the entire optical element.
    nx,ny = x.shape                         # dimension of the input 
    z = np.ones((nx,ny))*1e30               # initialize height profile
    npts = int(nx/npx)                      # dimension of a single prism
    dl = (x.max()-x.min())/(nx-1) * npts    # width of a single prism
    if npts != nx/npx:
        print('need integer multiples of small prism dimension')
    else:
        for il in range(npx):
            xi = il * npts              # indices for a single prism
            xf = (il+1)*npts
            for jl in range(npy):
                yi = jl * npts
                yf = (jl+1)*npts
                xl = x[xi:xf,yi:yf]     # x&y axis for a single prism
                yl = y[xi:xf,yi:yf]
                xlcent = xl.mean()      # single prism center position
                ylcent = yl.mean()
                z[xi:xf,yi:yf] = Prism(xl,yl,slope,d_hole,xlcent,ylcent)
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
            OPD[thicknessProfData>=1e30] = 0.
        #for iy in range(ny):
        #    for ix in range(nx):
        for iy in (np.arange(int(ny*0.4))+int(ny*0.3)):
            for ix in (np.arange(int(nx*0.4))+int(nx*0.3)):
                ofst = 2*ix + 2*nx*iy
                srwTr.arTr[ofst]=Tr[iy,ix]
                srwTr.arTr[ofst+1]=OPD[iy,ix]
    else:
        print('OE shape not matched')

