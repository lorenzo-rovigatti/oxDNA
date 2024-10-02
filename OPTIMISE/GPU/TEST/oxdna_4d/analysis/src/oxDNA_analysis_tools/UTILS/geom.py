import numpy as np

#I bet there's a builtin for this...
def fit_plane(points):
    #fits plane through points, return normal to the plane:
    rc= np.array(np.zeros(3))
    A = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
    for point in points:
        rc += point
    rc /= len(points)

    for point in points:
        point  = point - rc
        for i in range(3):
            for j in range(3):
                A[i][j] += point[i]*point[j]

    #this could be made faster by using a symmetric eigensolver
    vals, vecs = np.linalg.eigh(A) 		
    #print vals, vecs		
    return vecs[:,0]

def get_RNA_axis(particles, d):
    """
        Returns the axis of a RNA duplex
        
        Parameters:
            particles (oxpy.config_info): The positions/orientations of the particles
            start1 (int): Particle ID of the first particle in the first strand
            end1 (int): Particle ID of the last particle in the first strand
            end2 (int): Particle ID of the particle complimentary to start1
            start2 (int): Particle ID of the particle complimentary to end1    
    """

    nucA = particles[d.start1]
    nucB = particles[d.end2]
    midA = (nucA.backbone_site() + nucB.backbone_site()) / 2
    nucAc = particles[d.end1]
    nucBc = particles[d.start2]
    midAc = (nucAc.backbone_site() + nucBc.backbone_site()) / 2
    guess = (midAc - midA)
    guess /= np.linalg.norm(guess)

    vectorsA = []
    posAs = []
    posBs = []
    back_poses = []	

    for (i, j) in zip(range(d.start1, d.end1), range(d.end2, d.start2, -1)):
        nucA = particles[i]
        nucB = particles[j]
        nucAc = particles[i+1]
        nucBc = particles[j-1]
        midpointA = (nucA.backbone_site() + nucB.backbone_site()) / 2
        midpointB = (nucBc.backbone_site() + nucAc.backbone_site()) / 2
        vectorsA.append(midpointB - midpointA)
        posAs.append(nucA.backbone_site())
        posBs.append(nucB.backbone_site())
        if (i == d.end1 - 1):
            posAs.append(nucAc.backbone_site())
            posBs.append(nucBc.backbone_site())

        back_poses.append(nucAc.backbone_site() - nucA.backbone_site())
        back_poses.append(nucBc.backbone_site() - nucB.backbone_site())

    plane_vector = fit_plane(back_poses)
    
    if(np.dot(guess,plane_vector) < 0):
        plane_vector = -1 * plane_vector

    #if(np.rad2deg(np.arccos(np.dot(plane_vector,guess))) > 20):
    #    print('Warning, guess vector and plane vector have angles:', np.rad2deg(np.arccos(np.dot(guess,plane_vector))))

    #now we find the point where the helical vecotr originates
    hel_pos = []
    for i in range(len(posAs)-1):
        #project to the plane
        apos = posAs[i]
        bpos = posBs[i]
        apos = apos - np.dot(apos,plane_vector) * plane_vector
        bpos = bpos - np.dot(bpos,plane_vector) * plane_vector
        bp_vec = -apos + bpos
        if np.linalg.norm(bp_vec) == 0:
            continue
        midpointA = 0.5 * (apos + bpos)
        perpendicular_vecA = np.cross(bp_vec,plane_vector)
        perpendicular_vecA /= np.linalg.norm(perpendicular_vecA)
        apos = posAs[i+1]
        bpos = posBs[i+1]
        apos = apos - np.dot(apos,plane_vector) * plane_vector
        bpos = bpos - np.dot(bpos,plane_vector) * plane_vector
        bp_vec = -apos + bpos
        if np.linalg.norm(bp_vec) == 0:
            continue
        midpointB = 0.5 * (apos + bpos)
        perpendicular_vecB = np.cross(bp_vec,plane_vector)
        perpendicular_vecB /= np.linalg.norm(perpendicular_vecB)
        mat = np.array([perpendicular_vecA, -perpendicular_vecB ]).transpose()
        y = midpointB - midpointA
        t,c = np.linalg.lstsq(mat,y)[0]
        if np.linalg.norm(midpointA + t*perpendicular_vecA  -  (midpointB + c *perpendicular_vecB)) > 1.e-6 :
            print ('Error in finding common intersection point',midpointA + t*perpendicular_vecA , midpointB + c *perpendicular_vecB)
        hel_position = midpointA + t*perpendicular_vecA	
        #print ('Hel position from BP',i,'is ',hel_position)
        hel_pos.append(hel_position)	
                
    final_hel_pos = np.zeros(3)
    for pos in hel_pos:
        final_hel_pos += pos
    final_hel_pos /= len(hel_pos)	

    return plane_vector, final_hel_pos

def get_DNA_axis (particles, d):
    """
        Returns the axis of a DNA duplex
        
        Parameters:
            particles (oxpy.config_info): The positions/orientations of the particles
            start1 (int): Particle ID of the first particle in the first strand
            end1 (int): Particle ID of the last particle in the first strand
            end2 (int): Particle ID of the particle complimentary to start1
            start2 (int): Particle ID of the particle complimentary to end1    
    """
    vec = np.empty((d.end1-d.start1, 3))
    for i, j in zip(range(d.start1, d.end1), range(d.end2, d.start2, -1)):
        nucA = particles[i]
        nucAc= particles[j]
        midA = (nucA.base_site() + nucAc.base_site() ) / 2.
        nucB = particles[i+1]
        nucBc = particles[j-1] 
        midB = (nucB.base_site() + nucBc.base_site()) / 2.
        vec[i-d.start1] = (midB - midA)/(np.linalg.norm(midB - midA))

    #print(s._nucleotides[cfirst_base].index, s._nucleotides[cfirst_base].get_pos_base())

    pos = (particles[d.start1].base_site() + particles[d.end2].base_site() ) / 2.
    vector = np.mean(vec, axis=0)
    return(vector, pos)