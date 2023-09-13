import numpy as np
import sys


def get_M(Omegas,disc_len):
    """
        calculates the onsite stiffness matrix (m=1)
    """
    shape = np.shape(Omegas)
    N = shape[0]*shape[1]
    Cov = np.zeros([3,3])
    for i in range(len(Omegas)):
        for j in range(len(Omegas[0])):
            Cov+=np.outer(Omegas[i,j],Omegas[i,j])
    Cov = Cov/N
    Stiff = np.linalg.inv(Cov)/disc_len
    return Stiff
    

def get_mStep(Omegas,disc_len,m_max):
    """
        calculates the mstep stiffness matrix considering all possible m-step sums
    """
    shape = np.shape(Omegas)
    n_frames=shape[0]
    n_bps   =shape[1]
    if (m_max>n_bps):
        m_max=n_bps
        
    mCov    = np.zeros([m_max,3,3])
    count   = np.zeros(m_max)
    for i in range(n_frames):
        if (i%100==0):
            print("frame %d"%i)
        for j in range(10 , n_bps-10):
       #for j in range(n_bps):
            Gamma = np.zeros([3])
            for m in range(np.min([n_bps-j,m_max])):
                Gamma       += Omegas[i,j+m]
                mCov[m]     += np.outer(Gamma,Gamma)
                count[m]    += 1
    
    mStep = np.zeros([m_max,3,3])
    for m in range(m_max):
        mCov[m] /= count[m]
        mStep[m] = (m+1)*np.linalg.inv(mCov[m])/disc_len
        print("m=%d"%(m+1))
        print(mStep[m])
        filename = 'mstep'+ sys.argv[1]
        with open (filename , "a+") as f:
        	print >> f , '{:^8}{:^8}{:^8}{:^8}{:^8}{:^8}{:^8}'.format((m+1) , round(mStep[m][0,0],2) , round(mStep[m][1,1],2), round(mStep[m][2,2],2), round(mStep[m][2,1],2), round(mStep[m][1,0],2), round(mStep[m][2,0],2))
    return mStep[m]
 

#calculating m-step from correlations eq-3.31 enrico thesis

def get_M_correlations(Omegas,disc_len):
    """
        calculates the onsite stiffness matrix (m=1)
    """
    shape = np.shape(Omegas)
    N = shape[0]*shape[1]
    Cov = np.zeros([3,3])
    for i in range(len(Omegas)):
        for j in range(len(Omegas[0])):
            Cov+=np.outer(Omegas[i,j],Omegas[i,j])
    Cov = Cov/N
    Stiff = (Cov)
    return Stiff

    
def get_mStep_correlations(Omegas,disc_len,m_max,M_corr):
    """
        calculates the mstep stiffness matrix considering all possible m-step sums
    """
    shape = np.shape(Omegas)
    n_frames = shape[0]
    n_bps   = shape[1]
    if (m_max>n_bps):
        m_max=n_bps
        
    mCov    = np.zeros([m_max,3,3])
    count   = np.zeros(m_max)
    for i in range(n_frames):
        if (i%100==0):
            print("frame %d"%i)
        for m in range (m_max):
            for j in range (10, n_bps-m-1-10):
                mCov[m]        += np.outer(Omegas[i,j], Omegas[i,j+ m +1])
                count[m]       += 1
    
    mStep   = np.zeros([m_max,3,3])
    mCov1   = np.zeros([m_max,3,3])
    for m in range(m_max):
        mCov[m] /= count[m]
        for i in range (m+1):
		mCov1[m]  += (2+m-i-1)*(mCov[i] + np.transpose(mCov[i]))
        mStep[m]  = (m+2)*(np.linalg.inv(mCov1[m] + (m+2)*M_corr))/disc_len
        print("m=%d"%(m+2))
        print(mStep[m])
        filename = 'Correlations'
        with open (filename , "a+") as f:
        	print >> f , '{:^8}{:^8}{:^8}{:^8}{:^8}{:^8}{:^8}{:^8}{:^8}'.format((m+2) , round(mStep[m][0,0],2) , round(mStep[m][1,1],2), round(mStep[m][2,2],2), round(mStep[m][2,1],2), round(mStep[m][1,0],2), round(mStep[m][0,1],2) ,round(mStep[m][2,0],2),round(mStep[m][0,2],2))
    return mStep[m]
    


if __name__ == "__main__":
    """
    The dimensions of the Omegas array should be organized in the following order: snapshot, bp_id, dimension
    """
    if len(sys.argv) < 4:
        print("Wrong input!")
        print("usage: %s binary_filename discretization_length m_max n_snapshots (optional)"%sys.argv[0])
        sys.exit()
    OmegaBinary = sys.argv[1]
    disc_len    = float((sys.argv[2]))
    m_max       = int((sys.argv[3]))
    Omegas = np.load(OmegaBinary)
    shape= np.shape(Omegas)
    if len(sys.argv) > 4:
        n_snapshots = int((sys.argv[4]))
        if (n_snapshots<shape[0]):
            Omegas = Omegas[:n_snapshots]

    #~ M = get_M(Omegas,disc_len)
    mStep = get_mStep(Omegas,disc_len,m_max)
    #mStep = get_mStep_simple(Omegas,disc_len,m_max)
    #M_corr = get_M_correlations(Omegas,disc_len)
   # mStep_corr = get_mStep_correlations(Omegas,disc_len,m_max,M_corr)
    
