import base
import numpy as np
import sys
import subprocess
import pickle
import os
import tempfile

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

def get_pos_midpoint(r1, r2, box):
    """
    return the midpoint of two vectors r1 and r2
    
    use the minimum image: i.e. make sure we get a sensible answer if the                                                                     
    positions r1 and r2 correspond to nucleotides which were put at opposite                                                                  
    ends of the box due to pbc's                                                                                                            
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    return r1 - min_distance(r1, r2, box)/2

def array_index(mylist, myarray):
    """
    Return the index in a list mylist of a numpy array myarray
    """
    return map(lambda x: (myarray == x).all(), mylist).index(True)

def min_distance (r1, r2, box):
    """
    return the minimum image distance in going from r1 to r2, in a box of size box

    stolen from base.py Nucleotide.distance()
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    dr = r2 - r1
    dr -= box * np.rint (dr / box)
    return dr

def vecs2spline(vecs, per):
    import scipy.interpolate
    # interpolate vecs by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in vecs]
    yy = [vec[1] for vec in vecs]
    zz = [vec[2] for vec in vecs]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = per)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = per)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = per)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def get_base_spline(strand, reverse = False):
    """
    return a cartesian spline that represents a fit through the bases for the strand 'strand'

    args:
    strand: base.Strand object
    """
    import scipy

    base_pos = []
    for nuc in strand._nucleotides:
        base_pos.append(nuc.get_pos_base())

    if reverse:
        base_pos.reverse()

    if strand._circular:
        if reverse:
            base_pos.append(strand._nucleotides[-1].get_pos_base())
        else:
            base_pos.append(strand._nucleotides[0].get_pos_base())

    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in base_pos]
    yy = [vec[1] for vec in base_pos]
    zz = [vec[2] for vec in base_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = strand._circular)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = strand._circular)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = strand._circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def get_sayar_twist(s1, s2, smin, smax, npoints = 1000, circular = False, integral_type = "simple"):
    """
    return the twist for a given pair of spline fits, one through the bases of each strand

    from Sayar et al. 2010 Phys. Rev. E

    Just need to integrate along the contour parameter s that is common to both splines. We need the normalised tangent vector to the spline formed by the midpoint of the two input splines t(s), the normalised normal vector formed by the vectors between the splines u(s), and the derivative of the normalised normal vector between the splines d/ds (u(s)). NB, the normal u(s) vector should be orthogonal to the tangent vector t(s); we ensure this by using only the component orthogonal to t(s).

    Using integral_type = 'simple' and npoints = 200, it will give a correct twist, or at least one that gives a conserved linking number when combined with get_sayar_writhe

    args:
    s1: list of 3 splines corresponding to 3-D spline through strand 1's bases (e.g. use get_base_spline())
    s2: list of 3 splines corresponding to 3-D spline through strand 2's bases -- NB the splines should run in the same direction, i.e. one must reverse one of the splines if they come from get_base_spline (e.g. use get_base_spline(reverse = True))
    smin: minimum value for s, which parameterises the splines
    smax: maximum value for s, which parameterises the splines
    npoints: number of points for the discrete integration
    """

    import scipy.interpolate
    import scipy.integrate
    
    s1xx, s1yy, s1zz = s1
    s2xx, s2yy, s2zz = s2

    # bpi is the base pair index parameter that common to both splines
    bpi = np.linspace(smin, smax, npoints)

    # find the midpoint between the input splines, as a function of base pair index
    mxx = (scipy.interpolate.splev(bpi, s1xx) + scipy.interpolate.splev(bpi, s2xx)) / 2
    myy = (scipy.interpolate.splev(bpi, s1yy) + scipy.interpolate.splev(bpi, s2yy)) / 2
    mzz = (scipy.interpolate.splev(bpi, s1zz) + scipy.interpolate.splev(bpi, s2zz)) / 2

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((mxx[ii+1]-mxx[ii])**2+(myy[ii+1]-myy[ii])**2+(mzz[ii+1]-mzz[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    # get the midpoint spline as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, mxx, k = 3, s = 0, per = circular)
    msyy = scipy.interpolate.splrep(contour_len, myy, k = 3, s = 0, per = circular)
    mszz = scipy.interpolate.splrep(contour_len, mzz, k = 3, s = 0, per = circular)

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of normalised tangent vectors; __call__(xxx, 1) returns the first derivative
    # the tangent vector is a unit vector
    dmxx = scipy.interpolate.splev(ss, msxx, 1)
    dmyy = scipy.interpolate.splev(ss, msyy, 1)
    dmzz = scipy.interpolate.splev(ss, mszz, 1)
    tt = range(len(ss))
    for ii in range(len(ss)):
        tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])

    # we also need the 'normal' vector u(s) which points between the base pairs. (or between the spline fits through the bases in this case)
    # n.b. these uxx, uyy, uzz are not normalised
    uxx_bpi = scipy.interpolate.splev(bpi, s2xx) - scipy.interpolate.splev(bpi, s1xx)
    uyy_bpi = scipy.interpolate.splev(bpi, s2yy) - scipy.interpolate.splev(bpi, s1yy)
    uzz_bpi = scipy.interpolate.splev(bpi, s2zz) - scipy.interpolate.splev(bpi, s1zz)

    # get the normal vector spline as a function of contour length
    suxx = scipy.interpolate.splrep(contour_len, uxx_bpi, k = 3, s = 0, per = circular)
    suyy = scipy.interpolate.splrep(contour_len, uyy_bpi, k = 3, s = 0, per = circular)
    suzz = scipy.interpolate.splrep(contour_len, uzz_bpi, k = 3, s = 0, per = circular)

    # evaluate the normal vector spline as a function of contour length
    uxx = scipy.interpolate.splev(ss, suxx)
    uyy = scipy.interpolate.splev(ss, suyy)
    uzz = scipy.interpolate.splev(ss, suzz)

    uu = range(len(ss))
    for ii in range(len(ss)):
        uu[ii] = np.array([uxx[ii], uyy[ii], uzz[ii]])
        uu[ii] = uu[ii] - np.dot(tt[ii], uu[ii]) * tt[ii]
        # the normal vector should be normalised
        uu[ii] = norm(uu[ii])

    # and finally we need the derivatives of that vector u(s). It takes a bit of work to get a spline of the normalised version of u from the unnormalised one
    nuxx = [vec[0] for vec in uu]
    nuyy = [vec[1] for vec in uu]
    nuzz = [vec[2] for vec in uu]
    nusxx = scipy.interpolate.splrep(ss, nuxx, k = 3, s = 0, per = circular)
    nusyy = scipy.interpolate.splrep(ss, nuyy, k = 3, s = 0, per = circular)
    nuszz = scipy.interpolate.splrep(ss, nuzz, k = 3, s = 0, per = circular)
    duxx = scipy.interpolate.splev(ss, nusxx, 1)
    duyy = scipy.interpolate.splev(ss, nusyy, 1)
    duzz = scipy.interpolate.splev(ss, nuszz, 1)
    duu = range(len(ss))
    for ii in range(len(ss)):
        duu[ii] = np.array([duxx[ii], duyy[ii], duzz[ii]])

    ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
    # do the integration w.r.t. s
    if circular:
        srange = range(len(ss)-1)
    else:
        srange = range(len(ss))
    if integral_type == "simple":
        integral = 0
        for ii in srange:
            #print np.dot(uu[ii], tt[ii])
            triple_scalar_product = np.dot(tt[ii], np.cross(uu[ii], duu[ii]))
            integral += triple_scalar_product * ds
    elif integral_type == "quad":
        assert False, "not currently supported; shouldn't be difficult to implement if wanted"
        integral, err = scipy.integrate.quad(twist_integrand, ss[0], ss[-1], args = (msxx, msyy, mszz, nusxx, nusyy, nuszz), limit = 500)
        print >> sys.stderr, "error estimate:", err
        
    twist = integral/(2 * np.pi)

    return twist

def get_sayar_writhe(splines1, smin, smax, splines2 = False, npoints = 1000, debug = False, circular = False, integral_type = "simple"):
    """
    return the writhe for a 3D spline fit through a set of duplex midpoints

    from Sayar et al. 2010 Phys. Rev. E

    Using integral_type = 'simple' and npoints = 200, it will give a correct writhe, or at least one that gives a conserved linking number when combined with get_sayar_twist

    args:
    splines1: list of 3 splines corresponding to either (if not splines2) a 3D spline through the duplex or (if splines2) strand 1's bases
    smin: minimum value for s, which parameterises the splines
    smax: maximum value for s, which parameterises the splines
    splines2: optionally, (see splines1) list of 3 splines corresponding to a 3D spline through strand2's bases
    npoints: number of points for the discrete integration
    debug: print a load of debugging information
    """
    import scipy.integrate

    # bpi is the base pair index parameter that common to both strands' splines
    bpi = np.linspace(smin, smax, npoints)

    ## get midpoint splines sxx, syy, szz
    if not splines2:
        # splines1 is the midpoint 3D spline as a function of base pair index
        sxx_bpi, syy_bpi, szz_bpi = splines1
        xx_bpi = scipy.interpolate.splev(bpi, sxx_bpi)
        yy_bpi = scipy.interpolate.splev(bpi, syy_bpi)
        zz_bpi = scipy.interpolate.splev(bpi, szz_bpi)
    else:
        # take splines1 and splines2 to be the splines through the bases of each strand; in that case we need to find the midpoint here first
        s1xx_bpi, s1yy_bpi, s1zz_bpi = splines1
        s2xx_bpi, s2yy_bpi, s2zz_bpi = splines2

        # find the midpoint as a function of base pair index between the input splines 
        xx_bpi = (scipy.interpolate.splev(bpi, s1xx_bpi) + scipy.interpolate.splev(bpi, s2xx_bpi)) / 2
        yy_bpi = (scipy.interpolate.splev(bpi, s1yy_bpi) + scipy.interpolate.splev(bpi, s2yy_bpi)) / 2
        zz_bpi = (scipy.interpolate.splev(bpi, s1zz_bpi) + scipy.interpolate.splev(bpi, s2zz_bpi)) / 2

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((xx_bpi[ii+1]-xx_bpi[ii])**2+(yy_bpi[ii+1]-yy_bpi[ii])**2+(zz_bpi[ii+1]-zz_bpi[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    sxx = scipy.interpolate.splrep(contour_len, xx_bpi, k = 3, s = 0, per = circular)
    syy = scipy.interpolate.splrep(contour_len, yy_bpi, k = 3, s = 0, per = circular)
    szz = scipy.interpolate.splrep(contour_len, zz_bpi, k = 3, s = 0, per = circular)
    xx = scipy.interpolate.splev(ss, sxx)
    yy = scipy.interpolate.splev(ss, syy)
    zz = scipy.interpolate.splev(ss, szz)

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of tangent vectors; __call__(xxx, 1) returns the first derivative
    dxx = scipy.interpolate.splev(ss, sxx, 1)
    dyy = scipy.interpolate.splev(ss, syy, 1)
    dzz = scipy.interpolate.splev(ss, szz, 1)
    tt = range(len(ss))
    for ii in range(len(ss)):
        tt[ii] = np.array([dxx[ii], dyy[ii], dzz[ii]])
    
    # do the double integration w.r.t. s and s'
    if integral_type == "simple":
        integral = 0
        if circular:
            srange = range(len(ss)-1)
            ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
        else:
            srange = range(len(ss))
            ds = float(contour_len[-1] - contour_len[0]) / npoints
        for ii in srange:
            for jj in srange:
                # skip ii=jj and use symmetry in {ii, jj}
                if ii > jj:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = np.dot(np.cross(tt[ii], tt[jj]), diff_frac)
                    integral += triple_scalar_product * ds * ds
        # multiply by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
        integral *= 2
    elif  integral_type == "dblquad":
        # contour_len[0] to contour[-1] SHOULD be a complete integral on the closed curve; i.e. sxx(contour_len[0]) = sxx(contour_len[-1]) etc.
        val, err = scipy.integrate.dblquad(writhe_integrand, ss[0], ss[-1], lambda x: ss[0], lambda x: ss[-1], args = (sxx, syy, szz, ss[-1]))
        print >> sys.stderr, err
        integral = val
    elif integral_type == "chopped dblquad":
        integral = 0
        for ss_coarse in np.linspace(ss[0], ss[-1], 10):
            for ss_coarse_prime in np.linspace(ss[0], ss[-1], 10):
                val, err = scipy.integrate.dblquad(writhe_integrand, ss_coarse, ss_coarse + float(ss[-1]-ss[0])/9, lambda x: ss_coarse_prime, lambda x: ss_coarse_prime + float(ss[-1]-ss[0])/9, args = (sxx, syy, szz, contour_len[-1]))
                print err
                integral += val
    elif integral_type == "quad":
        integral, err = scipy.integrate.quad(writhe_integrand2, ss[0], ss[-1], args = (sxx, syy, szz, ss[0], ss[-1]), limit = 100, epsabs = 1e-5, epsrel = 0)
    elif integral_type == "simps":
        srange = range(len(ss))
        integrand = [[] for ii in srange]
        for ii in srange:
            for jj in srange:
                # skip ii=jj
                if ii == jj:
                    triple_scalar_product = 0
                else:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = np.dot(np.cross(tt[ii], tt[jj]), diff_frac)
                integrand[ii].append(triple_scalar_product)
        integral = scipy.integrate.simps(scipy.integrate.simps(integrand, ss), ss)
    else:
        assert False
        
    writhe = float(integral) / (4*np.pi)

    return writhe
    
def get_vhelix_vis(fname):
    """
    read the 'virtual helix visibility' from a text file and ignore any virtual helices set to invisible
    the syntax is the same as for the base.py visibility
    modified from the base.py strand visibility parsing
    """
    actions = {'vis' : True, 'inv' : False}
    visibility_list = []
    path = fname

    try:
        inp = open (path, 'r')
    except:
        base.Logger.log ("Origami visibility file `" + path + "' not found. Assuming default visibility", base.Logger.INFO)
        return False

    base.Logger.log("Using origami visibility file " +path, base.Logger.INFO)
    # read the visibility from a file; if we got here the file was opened
    lines = []
    for line in inp.readlines():
        line = line.strip().lower()
        # remove everything that comes after '#'
        line = line.split('#')[0]
        if len(line) > 0: lines.append(line)

    for line in lines:
        if '=' in line:
            sp = line.split("=")
            one, two, three = [p.strip() for p in sp[0], '=', sp[1]]
        else:
            one, two, three = line, "", ""

        if two != '=':
            base.Logger.log ("Lines in visibility must begin with one of inv=, vis= and default=. Skipping this line: --" + line + "--", base.Logger.WARNING)
            continue

        if one == 'default':
            if three not in ['inv', 'vis']:
                base.Logger.log ("Wrong default in visibility file. Assuming visible as default", base.Logger.WARNING)
                three = 'vis'
            if three == 'inv':
                mode = 'inv'
            else:
                mode = 'vis'
        else:
            # filter removes all the empty strings
            arr = [a.strip() for a in filter(None, three.split(','))]
            for a in arr:
                try:
                    ind = int(a)
                except:
                    base.Logger.log ("Could not cast '%s' to int. Assuming 0" % a, base.Logger.WARNING)
                    ind = 0
                try:
                    visibility_list.append(ind)
                except:
                    base.Logger.log ("vhelix %i does not exist in system, cannot assign visibility. Ignoring" % ind, base.Logger.WARNING)

    return mode, visibility_list

def norm(vec):
    return vec / np.sqrt(np.dot(vec,vec))

def parse_scanfile(infile):
    try:
        f = open(infile, "r")
    except IOError:
        base.Logger.log("could not open file %s" % infile, base.Logger.CRITICAL)
        sys.exit()
    for line in f.readlines():
        if line.startswith('begin_vb'):
            begin_vb = int(line.split()[2])
        elif line.startswith('end_vb'):
            end_vb = int(line.split()[2])
        elif line.startswith('vh_list'):
            vh_list = [int(x) for x in (line.split()[2]).split(",")]
        elif line.startswith('trim'):
            trim = int(line.split()[2])
    try:
        return begin_vb, end_vb, vh_list, trim
    except NameError:
        base.Logger.log("error while reading scan file %s, dying now" % infile, base.Logger.CRITICAL)
        sys.exit()

def get_scaffold_index(system):
    # find the index of the scaffold strand in the system
    strand_lengths = [strand.get_length() for strand in system._strands]
    return strand_lengths.index(max(strand_lengths))

def get_bb_midpoint(system, strand, n_index, interaction_list):
    base.Logger.log("origami_utils.get_bb_midpoint: deprecated function, use origami_utils.Origami.get_bb_midpoint", base.Logger.WARNING)
    # get midpoint vector between 2 hybridised bases
    r1 = strand._nucleotides[n_index].get_pos_base()
    r2 = system._nucleotides[interaction_list[strand._nucleotides[n_index].index]].get_pos_base()
    vec = (r1+r2)/2
    return vec

def get_nucleotide(vhelix, vbase, vhelix_indices):
    # find the system nucleotide index of a nucleotide given a position on the origami
    if vhelix % 2 == 0:
        dir = 1
    else:
        dir = -1
    return vhelix_indices[vhelix] + vbase * dir

def parse_vh_data(filename, origami):
    """
    get vhelices data from file - format is either <auto,[number of vhelices]> or, if a region of an origami is to be analysed, <[region width],[list of starting nucleotide index for the region for each vhelix]>
    """
    scaf_index = get_scaffold_index(origami._sys)
    vhelix_def_file = open(filename, "r")
    data = [x for x in vhelix_def_file.readline().replace("\n","").replace(" ","").split(",")]
    if data[0] == "auto":
        origami.num_vh = int(data[1])
        origami.width = origami._sys._strands[scaf_index].get_length() / origami.num_vh
        origami.vhelix_indices = []
        start_nuc_ind = -1
        for i in range(origami.num_vh):
            if i % 2 == 0:
                start_nuc_ind += 1
            else:
                start_nuc_ind += origami.width*2 - 1
            origami.vhelix_indices.append(start_nuc_ind)
    else:
        origami.width = int(data[0])
        origami.vhelix_indices = [int(x) for x in data[1:]]
        origami.num_vh = len(origami.vhelix_indices)
    base.Logger.log("using file data.vhd, %d virtual helices found" % origami.num_vh, base.Logger.INFO)

def print_arrow_debug_line(begin, end, file):
    if file:
        file.write("draw arrow {%f %f %f} {%f %f %f}\n" % (begin[0], begin[1], begin[2], end[0], end[1], end[2]))
    return 0

def open_arrow_debug_file(filename, type="w"):
    f_arrow = open(filename, type)
    f_arrow.write(
        "color Display Background white\n" +
        "set mName [mol new]\n" +
        "proc vmd_draw_arrow {mol start end} {\n" +
        "# an arrow is made of a cylinder and a cone\n"
        "set middle [vecadd $start [vecscale 0.65 [vecsub $end $start]]]\n" +
        "graphics $mol cylinder $start $middle radius 0.05\n" +
        "graphics $mol cone $middle $end radius 0.15\n" +
        "}\n")
    return f_arrow

def print_box_debug_line(vecs, file):
    if file:
        if len(vecs) == 4:
            file.write("draw box {%f %f %f} {%f %f %f} {%f %f %f} {%f %f %f}\n" % (vecs[0][0], vecs[0][1], vecs[0][2], vecs[1][0], vecs[1][1], vecs[1][2], vecs[2][0], vecs[2][1], vecs[2][2], vecs[3][0], vecs[3][1], vecs[3][2]))
        else:
            base.Logger.log("drawing boxes only works for 4 vertices at the moment", base.Logger.WARNING)

def open_box_debug_file(filename):
    f_boxes = open(filename, "w")
    f_boxes.write(
        "color Display Background white\n" +
        "set mName [mol new]\n" +
        "proc vmd_draw_box {mol vert1 vert2 vert3 vert4} {\n" +
        "# a 'box' is a plane made of 2 triangles here\n" +
        "graphics $mol triangle $vert1 $vert2 $vert3\n" +
        "graphics $mol triangle $vert1 $vert4 $vert3\n" +
        "}\n")
    return f_boxes

def angle_sense(v1, v2, axis):
    # return the angle between two vectors, using a 3rd vector to determine the sense of the angle
    v1 = norm(v1)
    v2 = norm(v2)
    axis = norm(axis)
    angle = np.arccos(np.dot(v1,v2))
    if np.dot(norm(np.cross(v1,v2)),axis) < 0:
        angle *= -1 # attempt to check the 'sense' of the angle w.r.t. an axis
    return angle

def dihedral_angle_sense(v1, v2, axis, sanity_check = False):
    # return the dihedral angle between two vectors, using a 3rd vector to determine the sense of the angle and to define the axis normal to the plane we want the vectors in
    v1n = norm(v1)
    v2n = norm(v2)
    axis = norm(axis)
    if sanity_check and (abs(np.dot(v1n,axis)) > 0.6 or abs(np.dot(v2n,axis)) > 0.6):
        return False
    # in plane
    v1p = v1n - np.dot(v1n,axis)*axis
    v2p = v2n - np.dot(v2n,axis)*axis
    v1p = norm(v1p)
    v2p = norm(v2p)
    angle = np.arccos(np.dot(v1p,v2p))
    if np.dot(norm(np.cross(v1p,v2p)),axis) < 0:
        angle *= -1 # attempt to check the 'sense' of the angle w.r.t. an axis
    return angle

class vhelix_vbase_to_nucleotide(object):
    # at the moment squares with skips in have entries in the dicts but with the nucleotide list empty (rather than having no entry) - I'm not sure whether or not this is desirable. It's probably ok
    def __init__(self):
        self._scaf = {}
        self._stap = {}
        self.nuc_count = 0 # record the nucleotide count, updated only after a whole strand is added
        self.strand_count = 0

    def add_scaf(self, vh, vb, strand, nuc):
        self._scaf[(vh, vb)] = (strand, nuc)

    def add_stap(self, vh, vb, strand, nuc):
        self._stap[(vh, vb)] = (strand, nuc)

    # these methods use a reference vhvb2n object to make the final vhvb2n object
    def add_scaf_strand(self, add_strand, reference, continue_join = False):
        count = 0
        size = len(self._scaf)
        for (vh, vb), [strand_ind, nuc] in reference._scaf.iteritems():
            if strand_ind == add_strand:
                self.add_scaf(vh, vb, self.strand_count, [x + self.nuc_count for x in nuc])
                count += len(nuc)
        self.nuc_count += count
        if len(self._scaf) == size:
            return 1
        else:
            if continue_join == False:
                self.strand_count += 1
            return 0

    def add_stap_strand(self, add_strand, reference, continue_join = False):
        count = 0
        size = len(self._stap)
        for (vh, vb), [strand_ind, nuc] in reference._stap.iteritems():
            if strand_ind == add_strand:
                self.add_stap(vh, vb, self.strand_count, [x + self.nuc_count for x in nuc])
                count += len(nuc)
        self.nuc_count += count
        if len(self._stap) == size:
            return 1
        else:
            if continue_join == False:
                self.strand_count += 1
            return 0

    def add_strand(self, add_strand, reference, continue_join = False):
        if self.add_scaf_strand(add_strand, reference, continue_join) and self.add_stap_strand(add_strand, reference, continue_join):
            #base.Logger.log("while adding strand %s to vhelix_vbase_to_nucleotide object; either strand already present or strand not found in reference object" % add_strand, base.Logger.WARNING)
            return 1
        else:
            return 0

class Origami(object):

    def __init__(self, system=False, cad2cuda_file = False, visibility = False):
        self.width = 0
        self.num_vh = 0
        if not system:
            base.Logger.log("origami_utils: Origami.__init__: system is False, this is no longer supported", base.Logger.CRITICAL)
            sys.exit(1)
        self._sys = system
        self.interaction_list = [-1 for x in range(self._sys._N)]
        self.vhelix_indices = []
        self.vbase_indices = []
        self.vec_long = np.array([0., 0., 0.])
        self.vec_lat = np.array([0., 0., 0.])
        self._vhelix_pattern = False
        self.vh_midpoints = []
        if cad2cuda_file:
            self.get_cad2cudadna(cad2cuda_file, visibility = visibility)
            # build list of complementary nucleotides (according to cadnano scheme)
            self.complementary_list = ["na" for x in range(self._sys._N)]
            for (vhelix, vbase), (strand1, nucs1) in self._cad2cudadna._scaf.iteritems():
                try:
                    (strand2, nucs2) = self._cad2cudadna._stap[vhelix, vbase]
                    for i in range(len(nucs1)):
                        # I'm pretty sure these should always line up for any possible insertion/deletion scheme
                        self.complementary_list[nucs1[i]] = nucs2[i]
                        self.complementary_list[nucs2[i]] = nucs1[i]
                except KeyError:
                    pass

            # build a list of the vhelix indices contained in the cadnano design
            self.vhelix_indices = sorted(list(set([key[0] for key in self._cad2cudadna._scaf.keys()] + [key[0] for key in self._cad2cudadna._stap.keys()])))
            if len(self.vhelix_indices) == 0:
                print >> sys.stderr, "WARNING: vhelix_indices member variable list is empty, this probably means the virt2nuc file you are using is out of date"
            self.vbase_indices = sorted(list(set([key[1] for key in self._cad2cudadna._scaf.keys()] + [key[1] for key in self._cad2cudadna._stap.keys()])))

            # build a list of the occupied virtual bases in each virtual helix
            self.vh_vbase_indices = [[] for x in self.vhelix_indices]
            self.vh_vbase_indices_scaf = [[] for x in self.vhelix_indices]
            self.vh_vbase_indices_stap = [[] for x in self.vhelix_indices]
            self.vvib = [[] for x in self.vhelix_indices] # vhelix_vbase_indices_both
            # scaffold vbase occupation
            for (vh, vb) in iter(self._cad2cudadna._scaf):
                self.vh_vbase_indices_scaf[self.vhelix_indices.index(vh)].append(vb)
            for row in self.vh_vbase_indices_scaf:
                row.sort()
            # staples vbase occupation
            for (vh, vb) in iter(self._cad2cudadna._stap):
                self.vh_vbase_indices_stap[self.vhelix_indices.index(vh)].append(vb)
            for row in self.vh_vbase_indices_stap:
                row.sort()
            # occupation of both staple and scaffold strand
            for vhi in range(len(self.vh_vbase_indices_scaf)):
                for vb_scaf in self.vh_vbase_indices_scaf[vhi]:
                    if vb_scaf in self.vh_vbase_indices_stap[vhi]:
                        self.vvib[vhi].append(vb_scaf)
            for row in self.vvib:
                row.sort()
            # occupation of either staple or scaffold strand
            for vhi in range(len(self.vh_vbase_indices_scaf)):
                for vb_scaf in self.vh_vbase_indices_scaf[vhi]:
                    if vb_scaf not in self.vh_vbase_indices[vhi]:
                        self.vh_vbase_indices[vhi].append(vb_scaf)
                for vb_stap in self.vh_vbase_indices_stap[vhi]:
                    if vb_stap not in self.vh_vbase_indices[vhi]:
                        self.vh_vbase_indices[vhi].append(vb_stap)
            for row in self.vh_vbase_indices:
                row.sort()
            # nicer aliases
            self.vvi = self.vh_vbase_indices
            self.vvisc = self.vh_vbase_indices_scaf
            self.vvist = self.vh_vbase_indices_stap

            self.num_vh = len(self.vhelix_indices)
            self.scaf_index = get_scaffold_index(self._sys)
#            self.width = self._sys._strands[scaf_index].get_length() / self.num_vh
#            if self.width != len(self.vbase_indices):
#                pass #this warning got annoying base.Logger.log("not a rectangular origami!", base.Logger.WARNING)

        else:
            self._cad2cudadna = {}


#        print self.vvib[0]
#        for ii in self.vvib[0]:
#            print self._cad2cudadna._stap[0,ii]
#        exit(1)
        
    def update_system(self, system):
        self._sys = system

    def get_corners(self):
        if self._cad2cudadna == {}:
            base.Logger.log("get_corners: build cad2cudadna property first", base.Logger.CRITICAL)
            sys.exit()

        # make sure that neighbours in the list are neighbours in the origami
        a = self.get_nucleotides(self.vhelix_indices[0],self.vh_vbase_indices[0][0])[0]
        b = self.get_nucleotides(self.vhelix_indices[0],self.vh_vbase_indices[0][-1])[0]
        c = self.get_nucleotides(self.vhelix_indices[-1],self.vh_vbase_indices[-1][-1])[0]
        d = self.get_nucleotides(self.vhelix_indices[-1],self.vh_vbase_indices[-1][0])[0]
        return [a, b, c, d]

    def get_nucleotides(self, vhelix, vbase, type="default"):
        # tries to return scaffold strand nucleotide, failing that staple strand, failing that error
        if self._cad2cudadna:
            if type == "default" or type == "single":
                try:
                    strand, nucs = self._cad2cudadna._scaf[(vhelix, vbase)]
                except KeyError:
                    strand, nucs = self._cad2cudadna._stap[(vhelix, vbase)]
                return nucs
            elif type == "scaf":
                try:
                    strand, nucs = self._cad2cudadna._scaf[(vhelix, vbase)]
                except KeyError:
                    nucs = []
                return nucs
            elif type == "stap":
                try:
                    strand, nucs = self._cad2cudadna._stap[(vhelix, vbase)]
                except KeyError:
                    nucs = []
                return nucs
            elif type == "double":
                # return double strands
                strand, nucs1 = self._cad2cudadna._scaf[(vhelix, vbase)]
                strand, nucs2 = self._cad2cudadna._stap[(vhelix, vbase)]
                nucs = []
                nucs.extend(nucs1)
                nucs.extend(nucs2)
                return nucs
        else:
            base.Logger.log("no cadnano to cudadna file detected, using old and possibly wrong get_nucleotides function", base.Logger.WARNING)
            # find the system nucleotide index of a nucleotide given a position on the origami
            if vhelix % 2 == 0:
                dir = 1
            else:
                dir = -1
            return self.vhelix_indices[vhelix] + vbase * dir

    def get_vhelix_ds_length(self, vhi):
        # requires one continuous double strand - otherwise how is double strand length defined??
        nucleotide_count = 0
        for vb in self.vvib[vhi]:
            try:
                nucs = self.get_nucleotides(self.vhelix_indices[vhi], vb, type="double")
            except KeyError:
                continue
            for nuc in nucs:
                nucleotide_count += 1
        return nucleotide_count/2


    def get_flat_nucs(self, vh):
        vhi = self.vhelix_indices.index(vh)
        nucs_flat = []
        iterable = self.vh_vbase_indices[vhi]
        for vb in iterable:
            try:
                nucs = self.get_nucleotides(vh, vb)
            except:
                continue #pass
            nucs_flat.extend(nucs)

        return nucs_flat

    def vb2nuci(self, vhi, vb0):
        vh = self.vhelix_indices[vhi]
        nucs_count = []
        vbi = self.vvib[vhi].index(vb0)
        iterable = self.vvib[vhi][:vbi]
        for vb in iterable:
            try:
                nucs = self.get_nucleotides(vh,vb)
            except:
                continue
            nucs_count.extend(nucs)

        return len(nucs_count)

    def vb2vhelix_nucid(self, vh, vb, vb_nuci, mode='default'):
        """
        given a virtual helix, a virtual base, and a nucleotide id for that virtual base, return the index of that nucleotide in a scheme where the leftmost (lowest virtual base index) nucleotide has index 0 and all nucleotides are given a unique integer index

        vh: virtual helix
        vb: virtual base
        vb_nuci: index of particular nucleotide in the virtual base (if the virtual base does not have a loop there is always exactly one nucleotide, which has index 0)
        mode: gets given to get_nucleotides as type; 'default' will count a vb that has either a scaffold or staple strand, while 'double' will only count vbs that have both a scaffold and staple strand. 'default' is the default
        """
        vhi = self.vhelix_indices.index(vh)
        nucs_count = 0
        this_vb = self.vvib[vhi][0]
        # loop over all virtual bases, adding up the number of nucleotides for each virtual base, until we get to the virtual base of interest
        while this_vb < vb:
            nucs_count += len(self.get_nucleotides(vh, this_vb, type=mode))
            this_vb += 1

        # this line may look pointless, but it allows correct behaviour even if vb_nuci is -1
        last_vb_nuc_index = self.get_nucleotides(vh, this_vb, type=mode).index(self.get_nucleotides(vh, this_vb, type=mode)[vb_nuci])
        # add the virtual-base-nucleotide-index of the last virtual base to our running total to get the final result
        vhelix_nucid = nucs_count + last_vb_nuc_index

        return vhelix_nucid
                           
    def prepare_principal_axes_calc():
        # get the central nucleotides that will (may) be used as reference points for the principal axes calculations
        com = np.array([0.,0.,0.])
        for nuc in self._sys._nucleotides:
            com += nuc.cm_pos
        com /= self._sys._N

        displ = range(self._sys._N)
        for ii in range(self._sys._N):
            rr = self._sys._nucleotides[ii].cm_pos - com
            displ[ii] = np.dot(rr,rr)

        minnuc = displ.index(min(displ))

        self._sys.map_nucleotides_to_strands()
        if minnuc+1 < self._sys._N and self._sys._nucleotide_to_strand[minnuc] == self._sys._nucleotide_to_strand[minnuc+1]:
            minnucs = (minnuc, minnuc+1)
        elif minnuc > 0 and self._sys._nucleotide_to_strand[minnuc] == self._sys._nucleotide_to_strand[minnuc-1]:
            minnucs = (minnuc, minnuc-1)
        else:
            base.Logger.die("something went wrong when trying to find the reference nucleotides for the principal axes calculation")
        return minnucs

    def get_principal_axes(self, approxaxes, minnucs):
        print "origami_utils.py: Origami.get_principal_axes: unsupported function, dying"
        sys.exit()
        # find moment of inertia tensor I, then eigenvectors are the principal axes
        # first get centre of mass
        com = np.array([0.,0.,0.])
        for nuc in self._sys._nucleotides:
            com += nuc.cm_pos
        com /= self._sys._N

        # get global rotation, the rotation to ensure that a particular vector always lies on [1,0,0] for every configuration
        vecref = norm(self._sys._nucleotides[minnucs[0]].cm_pos - self._sys._nucleotides[minnucs[1]].cm_pos)
        globalrotcg = 0#mat3().fromToRotation([vecref[0], vecref[1], vecref[2]], [1.,0.,0.])
        globalrot = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
        # convert to numpy array
        for ii in range(3):
            for jj in range(3):
                globalrot[ii][jj] = globalrotcg[ii][jj]
        # find I wrt centre of mass
        I = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        for nuc in self._sys._nucleotides:
            # rotate system so that vecref along x axis
            rk = np.dot(globalrot,(nuc.cm_pos - com))
            I[0,0] += rk[1]*rk[1] + rk[2]*rk[2]
            I[1,1] += rk[0]*rk[0] + rk[2]*rk[2]
            I[2,2] += rk[0]*rk[0] + rk[1]*rk[1]
            I[0,1] -= rk[0]*rk[1]
            I[0,2] -= rk[0]*rk[2]
            I[1,2] -= rk[1]*rk[2]
            I[1,0] -= rk[0]*rk[1]
            I[2,0] -= rk[0]*rk[2]
            I[2,1] -= rk[1]*rk[2]

        eigvals, eigvecs = np.linalg.eig(I)

        # order eigenvectors by size of eigenvalue
        i1 = np.where(eigvals == max(eigvals))[0][0]
        i3 = np.where(eigvals == min(eigvals))[0][0]
        i2 = 3 - i1 - i3
        v1 = eigvecs[:,i1]
        v2 = eigvecs[:,i2]
        v3 = eigvecs[:,i3]

        # order eigenvectors by how much they overlap with the approximate axes
        eigvecsl = [eigvecs[:,i] for i in range(len(eigvecs))] # make a list of 1d arrays from a 2d array
        eigvecs_ordered = []
        for ii in range(len(approxaxes)):
            res = range(len(eigvecsl))
            for jj in range(len(eigvecsl)):
                res[jj] = abs(np.dot(approxaxes[ii], norm(eigvecsl[jj]))) # to allow antiparallel vectors to count as aligned - not sure whether this is the best method

            eigvecs_ordered.append(eigvecsl.pop(res.index(max(res))))

        return np.array([eigvecs_ordered[0], eigvecs_ordered[1], eigvecs_ordered[2]])

    def get_1d_vecs(self, discard_unbonded = True):
        # returns orthonormal vectors based on the orientation of the strand
        # "invariant" co-ordinates x', y', z'. x': along double helix; yy: average backbone-base vector for the 0th strand; z': x' cross yy; y': z' cross x'

        trim = 3 # ignore ends of the double strand; for very short strands this may not be desirable

        vh = self.vhelix_indices[0]
        nucs_flat = self.get_flat_nucs(vh)

        # x'
        if len(nucs_flat) > trim * 2:
            n1 = nucs_flat[trim]
            n2 = nucs_flat[-trim]
        else:
            n1 = nucs_flat[0]
            n2 = nucs_flat[0]

        # check they exist # surely we already know they do??? 13/06/12
        bbm1 = self.get_bb_midpoint(n1)
        bbm2 = self.get_bb_midpoint(n2)
        if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
            xprime = bbm2 - bbm1

        # yy
        yy = np.zeros(3)
        if len(nucs_flat) > trim * 2:
            nuciter = nucs_flat[trim:-trim]
        else:
            nuciter = nucs_flat
        for nuci in nuciter:
            nuc = self._sys._nucleotides[nuci]
            yy += nuc._a1

        zprime = np.cross(xprime, yy)
        yprime = np.cross(zprime, xprime)

        xprime = norm(xprime)
        yprime = norm(yprime)
        zprime = norm(zprime)
        m = np.array( [xprime,yprime,zprime] )

        return m

    def get_3d_vecs(self, vhelix_extent = False, discard_unbonded = True):
        # returns orthonormal vectors based on the orientation of the origami.
        # z: average vector along helices; y: average vec between vertical helices; x: y cross z
        # vhelix_extent currently unused.... uses extent of double helices (assuming there is one duplex per virtual helix)
        # assumes the skips are lined up across vhelices - in particular it will presumably break if there is a vb with a skip adjacent (in vhelix dimension) to a vb that represents the end of a double helix

        if not self._vhelix_pattern:
            base.Logger.log("origami_utils: Origami.get_3d_vecs(): no virtual helix pattern found; either missing virt2nuc file in this directory, or virt2nuc file is out of date; Aborting now", base.Logger.CRITICAL)
            sys.exit(1)

        if vhelix_extent:
            vh_begin = vhelix_extent[0]
            vh_end = vhelix_extent[1]

        y = np.zeros(3, dtype = "float64")
        z = np.zeros(3, dtype = "float64")
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            # z: vec along helices
            if vhelix_extent:
                #print vhi, len(vhelix_extent)
                vb1 = vh_begin[vhi]
                vb2 = vh_end[vhi]
            else:
                vb1 = self.vvib[vhi][0]
                vb2 = self.vvib[vhi][-1]
            n1 = self.get_nucleotides(vh, vb1)[0]
            n2 = self.get_nucleotides(vh, vb2)[-1]
            if (self.interaction_list[n1] == -1 or self.interaction_list[n2] == -1) and discard_unbonded:
                base.Logger.log("unbonded base when calculating interaction between %d,%d and %d,%d; interaction will not be used" % (vh, 0, vh, self.width - 1), base.Logger.WARNING)
            else:
                bbm1 = self.get_bb_midpoint(n1)
                bbm2 = self.get_bb_midpoint(n2)
                if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                    z += bbm2 - bbm1

            # y: vec across virtual helices that are arranged vertically on the honeycomb lattice
            row = self._vhelix_pattern[vh][0]
            col = self._vhelix_pattern[vh][1]
            # check for a virtual helix above the current one
            result = [vh1 for vh1, (row1,col1) in self._vhelix_pattern.iteritems() if row1 == row - 1 and col1 == col]
            if len(result) > 0:
                # a virtual helix above the current one was found
                vhin = result[0]
                #print row, col, vhi, vhin
                vb1 = max(self.vvib[vhi][0], self.vvib[vhin][0])
                vb2 = min(self.vvib[vhi][-1], self.vvib[vhin][-1])
                for vb in (vb1, vb2):
                    lnucs = len(self.get_nucleotides(vh, vb))
                    if lnucs > 0:
                        if vb in self.vvib[vhin]:
                            n1 = self.get_nucleotides(self.vhelix_indices[vhi],vb)[0]
                            n2 = self.get_nucleotides(self.vhelix_indices[vhin],vb)[0]
                            bbm1 = self.get_bb_midpoint(n1)
                            bbm2 = self.get_bb_midpoint(n2)
                            if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                                y += bbm2 - bbm1

        x = np.cross(y,z)
        # redefine y to make sure it's perpendicular to z
        y = np.cross(z,x)

        x /= np.sqrt(np.dot(x,x))
        y /= np.sqrt(np.dot(y,y))
        z /= np.sqrt(np.dot(z,z))
        m = np.array([x,y,z])

        return m

    def get_plane_vecs(self, discard_unbonded = True):
        # first find average bbm-bbm vector along helices
        av_long = np.zeros(3, dtype = "float64")
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            vb1 = self.vvib[vhi][0]
            vb2 = self.vvib[vhi][-1]
            n1 = self.get_nucleotides(vh, vb1)[0]
            n2 = self.get_nucleotides(vh, vb2)[-1]
            '''
            for vbi in range(len(self.vvib[vhi][:-1])):
                vb = self.vvib[vhi][vbi]
                nucs = self.get_nucleotides(vh, vb)
                for nuci in range(len(nucs)):
                    # find the next nucleotide along
                    n1 = nucs[nuci]
                    try:
                        n2 = nucs[nuci+1]
                    except IndexError:
                        ii = vbi
                        while ii < len(self.vvib[vhi]):
                            try:
                                n2 = self.get_nucleotides(vh, vb+1)[0]
                                break
                            except:
                                ii += 1

            '''
            if (self.interaction_list[n1] == -1 or self.interaction_list[n2] == -1) and discard_unbonded:
                base.Logger.log("unbonded base when calculating interaction between %d,%d and %d,%d; interaction will not be used" % (vh, vb1, vh, vb2), base.Logger.WARNING)
            else:
                bbm1 = self.get_bb_midpoint(n1)
                bbm2 = self.get_bb_midpoint(n2)
                if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                    av_long += bbm2 - bbm1

        # next find average bbm-bbm vector across helices
        av_lat = np.zeros(3, dtype = "float64")
        for vb in self.vbase_indices:
            # get a list of all vhelices with that vbase occupied
            vhs = []
            for x in iter(self._cad2cudadna._scaf):
                if x[1] == vb:
                    vhs.append(x[0])
            vhs.sort()
            if len(vhs) < 2:
                continue
            # if there is a pattern of skips and loops that is not the same across vhelices, this code tries to cope with it but it's not the most rigorous way...
            nucs = self.get_nucleotides(vhs[0], vb)
            for nuci in range(len(nucs)):
                n1 = nucs[nuci]
                skip = False
                for jj in range(len(vhs))[::-1]:
                    try:
                        n2 = self.get_nucleotides(vhs[jj], vb)[nuci]
                        break
                    except:
                        if jj in [0,1]:
                            skip = True
                            break
                        else:
                            pass
                if skip:
                    continue
                if (self.interaction_list[n1] == -1 or self.interaction_list[n2] == -1) and discard_unbonded:
                    base.Logger.log("unbonded base when calculating interactions between %d,%d and %d,%d; interaction will not be used" % (0, vb, self.num_vh - 1, vb), base.Logger.WARNING)
                else:
                    bbm1 = self.get_bb_midpoint(n1)
                    bbm2 = self.get_bb_midpoint(n2)
                    if isinstance(bbm1, np.ndarray) and isinstance(bbm2, np.ndarray):
                        av_lat += bbm2 - bbm1

        av_long /= np.sqrt(np.dot(av_long, av_long))
        av_lat /= np.sqrt(np.dot(av_lat, av_lat))

        self.vec_long = av_long
        self.vec_lat = av_lat

    def get_bb_midpoint(self, n_index, pbc = True):
        # get midpoint 2 base vectors that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_base()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_base()
            if pbc:
                # make sure we get a sensible answer if the nucleotides got put at opposite ends of the box due to pbc's
                vec = r1 + min_distance(r1, r2, self._sys._box)/2
            else:
                vec = (r1+r2)/2
            return vec

    def get_backback_midpoint(self, n_index):
        # get midpoint of 2 backbone vectors corresponding to nucleotides that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_back()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_back()
            vec = (r1+r2)/2
            return vec

    def get_bb_vec(self, n_index):
        # get vector between 2 bases that are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_base()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_base()
            vec = r1-r2
            return vec

    def get_backback_vec(self, n_index):
        # get vector between 2 backbones whose nucleotides are hybridised according to cadnano scheme
        if self.complementary_list[n_index] == "na":
            return False
        else:
            r1 = self._sys._nucleotides[n_index].get_pos_back()
            r2 = self._sys._nucleotides[self.complementary_list[n_index]].get_pos_back()
            vec = r1-r2
            return vec

    def get_flat_units(self, discard_unbonded):
        # find width, height units of the origami if it had been laid out flat
        flat_width = 0
        flat_height = 0
        for vh in range(self.num_vh):
            for vb in range(self.width - 1):
                n1 = self.get_nucleotide(vh, vb)
                n2 = self.get_nucleotide(vh, vb + 1)
                if (self.interaction_list[n1] != -1 and self.interaction_list[n2] != -1) or not discard_unbonded:
                    dist = self.get_bb_midpoint(n1) - self.get_bb_midpoint(n2)
                dist = np.sqrt(np.dot(dist, dist))
                flat_width += dist
        flat_width /= self.num_vh * (self.width - 1)

        for vh in range(self.num_vh - 1):
            for vb in range(self.width):
                n1 = self.get_nucleotide(vh, vb)
                n2 = self.get_nucleotide(vh + 1, vb)
                if (self.interaction_list[n1] != -1 and self.interaction_list[n2] != -1) or not discard_unbonded:
                    dist = self.get_bb_midpoint(n1) - self.get_bb_midpoint(n2)
                dist = np.sqrt(np.dot(dist, dist))
                flat_height += dist
        flat_height /= (self.num_vh - 1) * self.width
        return flat_width, flat_height

    def get_com(self):
        # unsupported/may no longer work - e.g. take care in using self.width and get_bb_midpoint
        com = np.zeros(3)
        for vh in range(self.num_vh):
            for vb in range(self.width):
                com += self.get_bb_midpoint(self.get_nucleotide(vh, vb))
        com /= self.num_vh * self.width
        return com

    def get_h_bond_list_output_bonds(self, infile, conffile, conf_num):
        """
        the old (slow) version using output bonds
        """
        system = self._sys

        system.map_nucleotides_to_strands()
        try:
            open(infile)
        except:
            base.Logger.log("unable to find file %s, exit" % infile, base.Logger.CRITICAL)
            sys.exit()
        launchargs = [PROCESSDIR + 'output_bonds',infile,conffile,str(conf_num)]
        myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        system.read_H_bonds_output_bonds(myinput.stdout.readlines())

        self.interaction_list = [[] for x in range(self._sys._N)]
        for nucleotide in system._nucleotides:
            for i in nucleotide.interactions:
                self.interaction_list[nucleotide.index].append(i)

    def get_h_bond_list(self, infile):
        """
        Fill in self.interaction_list for the current system self._sys; self.interaction_list[nucid1] is a list of nucleotides that nucid1 has hydrogen bonds with
        """

        # DNAnalysis command
        command_for_data =  'analysis_data_output_1 = { \n name = stdout \n print_every = 1 \n col_1 = { \n type=pair_energy \n} \n}'
        PROCESSPROGRAM = os.path.join(os.path.dirname(__file__), "../build/bin/DNAnalysis")
        tempfile_obj = tempfile.NamedTemporaryFile()
        launchargs = [PROCESSPROGRAM,infile ,'trajectory_file='+tempfile_obj.name,command_for_data]

        self._sys.map_nucleotides_to_strands()
        
        try:
            open(infile)
        except:
            base.Logger.log("unable to find file %s, exit" % infile, base.Logger.CRITICAL)
            sys.exit()
            
        # print system to temporary file so we don't have to scan through a huge trajectory file...
	self._sys.print_lorenzo_output(tempfile_obj.name,'/dev/null')
	tempfile_obj.flush()

        # run DNAnalysis
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout,stderr = myinput.communicate()
	linewise = stdout.split('\n')
	self._sys.read_H_bonds(linewise[:-2])
        # check for errors in DNAnalysis
	for line in stderr.split('\n'):
      	  if "CRITICAL" in line:
              	  print line

        # fill in self.interaction_list
        self.interaction_list = [[] for x in range(self._sys._N)]
        for nucleotide in self._sys._nucleotides:
            for i in nucleotide.interactions:
                self.interaction_list[nucleotide.index].append(i)

    def get_cad2cudadna(self, infile, visibility = False):
        f = open(infile, "r")
        data = pickle.load(f)
        if type(data) is tuple:
            self._cad2cudadna = data[0]
            self._vhelix_pattern = data[1]
            # set the visibility if necessary
            if visibility:
                vis_res = get_vhelix_vis(visibility)
                if vis_res:
                    vis_mode, vhelix_vis = vis_res
                else:
                    vhelix_vis = False
                if vhelix_vis:
                    vis_c2cdna = vhelix_vbase_to_nucleotide()
                    # all invisible except those listed; if they are in the list, make them visible
                    if vis_mode == 'inv':
                        for vhid in self._vhelix_pattern.keys():
                            if vhid not in vhelix_vis:
                                del self._vhelix_pattern[vhid]
                        for vis_id in vhelix_vis:
                            for vh, vb in self._cad2cudadna._scaf.keys():
                                if vh == vis_id:
                                    vis_c2cdna._scaf[(vh, vb)] = self._cad2cudadna._scaf[(vh, vb)]
                            for vh, vb in self._cad2cudadna._stap.keys():
                                if vh == vis_id:
                                    vis_c2cdna._stap[(vh, vb)] = self._cad2cudadna._stap[(vh, vb)]
                    # all visible except those listed; if they are not in the list, make them visible
                    else:
                        assert vis_mode == 'vis'
                        for vhid in self._vhelix_pattern.keys():
                            if vhid in vhelix_vis:
                                del self._vhelix_pattern[vhid]
                        for vh, vb in self._cad2cudadna._scaf.keys():
                            if vh not in vhelix_vis:
                                vis_c2cdna._scaf[(vh, vb)] = self._cad2cudadna._scaf[(vh, vb)]
                        for vh, vb in self._cad2cudadna._stap.keys():
                            if vh not in vhelix_vis:
                                vis_c2cdna._stap[(vh, vb)] = self._cad2cudadna._stap[(vh, vb)]
                    self._cad2cudadna = vis_c2cdna
        else:
            self._cad2cudadna = data
        f.close()

    def update_system(self, system):
        self._sys = system

    def map_base_to_vbase(self, vhelices_len, begin_bases, end_bases):
        # not used right now
        # returns virtual base AND the index of the nucleotide at that virtual base (ie to deal with the possiblity of a loop)
        self._base_to_vbase = [[] for x in range(vhelices_len)]
        for vhelix in range(vhelices_len):
            for mybase in range(begin_bases[vhelix], end_bases[vhelix]+1):
                current_base = -1
                current_vbase = -1
                while current_base < mybase:
                    current_vbase += 1
                    try:
                        strand, nucs = self._cad2cudadna._stap[(vhelix, current_vbase)]
                    except KeyError:
                        nucs = []
                    for i in range(len(nucs)):
                        current_base += 1
                        if current_base == mybase:
                            self._base_to_vbase[vhelix].append((current_vbase, i))
                            break
                    if current_vbase > (end_bases[vhelix]+1) * 2:
                        base.Logger.die("Unable to find virtual base in map_base_to_vbase, dying now")

    def base_to_vbase(self, vhelix, base):
        # not used right now
        return self._base_to_vbase[(vhelix, base)]

    def vbase_to_base(self, vhelix, vbase):
        # not used right now
        # return base corresponding to first base at a particular virtual base i.e. taking into account skips and loops (counting from the left of a cadnano diagram)
        base = 0
        for i in range(vbase - 1):
            try:
                strand, nucs = self._cad2cudadna._stap[(vhelix, i)]
            except KeyError:
                nucs = []
            base += len(nucs)
        base += 1
        return base

    def square_contains_exactly_one_nucleotide(self, vh, vb, strand_type):
        """
        Does the square vh, vb contain exactly one nucleotide?
        """
        if strand_type == "staple":
            try:
                strand, nucs = self._cad2cudadna._stap[(vh,vb)]
            except KeyError:
                return False
        else:
            assert strand_type == "scaffold"
            try:
                strand, nucs = self._cad2cudadna._scaf[(vh,vb)]
            except KeyError:
                return False
        return len(nucs) == 1

    def is_crossover(self, vh, vhn, vb, strand_type):
        """
        Does a crossover of type strand_type span the squares (vh, vb) and (vhn, vb)?
        """
        # for a crossover we have no insertions/deletions
        if (self.square_contains_exactly_one_nucleotide(vh, vb, strand_type) and self.square_contains_exactly_one_nucleotide(vhn, vb, strand_type)):
            if strand_type == "staple":
                cstrand, [cnuc] = self._cad2cudadna._stap[(vh,vb)]
                nstrand, [nnuc] = self._cad2cudadna._stap[(vhn,vb)]
            else:
                assert strand_type == "scaffold"
                cstrand, [cnuc] = self._cad2cudadna._scaf[(vh,vb)]
                nstrand, [nnuc] = self._cad2cudadna._scaf[(vhn,vb)]
            # if there is a crossover between the squares (vh,vb) and (vhn,vb), the strand will be the same and the nucleotides will be consecutive
            adjacent_nucs = abs(cnuc - nnuc) == 1 or abs(cnuc - nnuc) == self._sys._strands[get_scaffold_index(self._sys)].get_length() - 1
            return (cstrand == nstrand and adjacent_nucs)
        else:
            return False

    def get_holliday_junctions(self, single=False, single_only=False):
        """
        returns list of hjs in format (vh, vb, vh_neighbour, vb_next)

        set single=True to get a list of all single crossovers rather than holliday junctions ie pairs of crossovers
        set single_only = True to get a list of just single crossovers (double crossovers excluded). If single_only is true, the state of single should not matter

        if somehow there were two holliday junctions between the same pairs of vbases (one for scaffold and one for staple strand), the 2nd would not be detected here
        assumes no insertions/deletions when using cad2cudadna
        """
        vh_pattern = self._vhelix_pattern
        vh_neighbours_below = [[] for x in range(max(vh_pattern.keys())+1)]
        for vh, (row,col) in vh_pattern.iteritems():
            neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and (col1 - col) == 1) or ((row1 - row) == 1 and col1 == col)]
            for vh1 in neighbours:
                vh_neighbours_below[vh].append(vh1)

        hjs = []
        skip_scaffold = False
        skip_staple = False
        lone_crossover = 0
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            if len(vh_neighbours_below[vh]) > 0:
                vhn = vh_neighbours_below[vh][0]
                for ii in range(len(self.vvib[vhi])):
                    vb = self.vvib[vhi][ii]
                    vbn = vb + 1
                    # check for staple single/double crossover
                    if self.is_crossover(vh, vhn, vb, "staple") and not skip_staple:
                        if self.is_crossover(vh, vhn, vbn, "staple"):
                            if not single_only:
                                hjs.append([vh, vb, vhn, vbn])
                            skip_staple = True # skip the next vbase since it is part of a double crossover
                        else:
                            lone_crossover += 1
                            if single or single_only:
                                hjs.append([vh, vhn, vb])
                    else:
                        skip_staple = False

                    # check for scaffold single/double crossover
                    if self.is_crossover(vh, vhn, vb, "scaffold") and not skip_scaffold:
                        if self.is_crossover(vh, vhn, vbn, "scaffold"):
                            if not single_only:
                                hjs.append([vh, vb, vhn, vbn])
                            skip_scaffold = True # skip the next vbase since it is part of a double crossover
                        else:
                            lone_crossover += 1
                            if single or single_only:
                                hjs.append([vh, vhn, vb])
                    else:
                        skip_scaffold = False

        if not single and not single_only and lone_crossover > 0:
            base.Logger.log("%d lone crossovers found, they will not be used in analysis" % lone_crossover, base.Logger.INFO)
        return hjs

    def get_vhelix_neighbours(self, type):
        # checked for square lattice, not yet double checked for honeycomb lattice
        if type not in ("he", "sq"):
            base.Logger.log("origami_utils: get_vhelix_neighbours(): error while building neighbour list; unknown lattice type %s, dying now" % type, base.Logger.CRITICAL)
            sys.exit()
        vh_pattern = self._vhelix_pattern
        vh_neighbour_list = [[] for x in range(max(vh_pattern.keys())+1)]
        for vh, (row,col) in vh_pattern.iteritems():
            if type == "sq":
                neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and abs(col1 - col) == 1) or (abs(row1 - row) == 1 and col1 == col)]
            elif type == "he":
                if vh % 2 == 0:
                    neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and abs(col1 - col) == 1) or (row1 == row - 1 and col1 == col)]
                if vh % 2 == 1:
                    neighbours = [vh1 for vh1, (row1,col1) in vh_pattern.iteritems() if (row1 == row and abs(col1 - col) == 1) or (row1 == row + 1 and col1 == col)]
            for vh1 in neighbours:
                vh_neighbour_list[vh].append(vh1)

        return vh_neighbour_list

    def get_local_twist(self, vh_id, vvib_id, type, conf=False):
        # find local twist between given vhelix, vbase and the next one along
        # assumes no skip/loop / insertion/deletion
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id+1])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_local_twist: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]
        # get the 3 normalised vectors we need: the base-base (or back-back) vectors and the vector from one base-base (or back-back) midpoint to the other
        if type == "base":
            r1 = ss._nucleotides[n1].get_pos_base()
            r2 = ss._nucleotides[n2].get_pos_base()
            r3 = ss._nucleotides[n3].get_pos_base()
            r4 = ss._nucleotides[n4].get_pos_base()
        elif type == "back":
            r1 = ss._nucleotides[n1].get_pos_back()
            r2 = ss._nucleotides[n2].get_pos_back()
            r3 = ss._nucleotides[n3].get_pos_back()
            r4 = ss._nucleotides[n4].get_pos_back()
        else:
            base.Logger.die("origami_utils.Origami.get_local_twist: unknown type %s; use either base or back" % type)

        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        n = mid34 - mid12
        v12 = r2 - r1
        v34 = r4 - r3
        n = norm(n)
        v12 = norm(v12)
        v34 = norm(v34)
        # find the component of the base-base (or back-back) vectors that is in the plane normal to n
        v12prime = v12 - np.dot(v12, n) * n
        v34prime = v34 - np.dot(v34, n) * n
        v12prime = norm(v12prime)
        v34prime = norm(v34prime)


        bp_twist = np.arccos(np.dot(v12prime, v34prime))
        bp_twist *= 180./np.pi
        xp = np.cross(v12prime, v34prime)
        if np.dot(xp,n) < 0:
            # this is interesting to note: happens fairly often in edge bp (of course we should be discarding those ones anyway)
            pass #print "not aligned", bp_twist, vh_id, vvib_id, conf
        return bp_twist

    def get_radius(self, vh_id, vvib_id, new_method=False):
        # find radius of a base pair (i.e. to find radius of helix)
        # assumes no skip/loop
        if new_method:
            # gives the same result for an ideal config but gives greater radii for an equilibrated config
            ss = self._sys
            n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
            n2 = self.complementary_list[n1]
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id+1])[0]
            n4 = self.complementary_list[n3]
            # use backbone centres
            back1 = ss._nucleotides[n1].get_pos_back()
            back2 = ss._nucleotides[n2].get_pos_back()
            r1 = ss._nucleotides[n1].get_pos_base()
            r2 = ss._nucleotides[n2].get_pos_base()
            r3 = ss._nucleotides[n3].get_pos_base()
            r4 = ss._nucleotides[n4].get_pos_base()
            # rr is the vector from one base-base midpoint to the next one
            bbm1 = (r1+r2)/2
            bbm2 = (r3+r4)/2
            rr = bbm2 - bbm1
            rr = rr/np.sqrt(np.dot(rr,rr))
            # vectors d1, d2 from bbmidpoint to backbone sites
            d1 = back1 - bbm1
            d2 = back2 - bbm1
            d1_in_plane = d1 - np.dot(d1,rr)
            d2_in_plane = d2 - np.dot(d2,rr)
            # average over these two values of centre-to-backbone distance
            d1_in_plane = np.sqrt(np.dot(d1_in_plane,d1_in_plane))
            d2_in_plane = np.sqrt(np.dot(d2_in_plane,d2_in_plane))
            d = (d1_in_plane + d2_in_plane)/2
            # include the excluded volume of the backbone site
            d += base.EXCL_S1/2
            
        else:
            ss = self._sys
            n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
            n2 = self.complementary_list[n1]
            # use backbone centres
            r1 = ss._nucleotides[n1].get_pos_back()
            r2 = ss._nucleotides[n2].get_pos_back()
            d = r1-r2
            d = np.sqrt(np.dot(d,d))/2 # radius
            d += base.EXCL_S1/2

        return d

    def get_rise(self, vh_id, vvib_id):
        # find rise for a base pair
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id+1])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_rise: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34
        d = np.sqrt(np.dot(d,d))

        return d

    def get_bpwise_weave(self, vh_id, vvib_id):
        # find weave for a pair of base pairs located at (vh_id, vvib_id) and (vh_id + 1, vvib_id)
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(self.vhelix_indices[vh_id], self.vvib[vh_id][vvib_id])[0]
        n2 = self.complementary_list[n1]
        try:
            n3 = self.get_nucleotides(self.vhelix_indices[vh_id+1], self.vvib[vh_id+1][vvib_id])[0]
        except KeyError:
            base.Logger.die("origami_utils.Origami.get_bpwise_weave: KeyError - probably the virtual base index argument was too high")
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34
        d = np.sqrt(np.dot(d,d))

        return d

    def get_weave(self, verbose=False):
        """
        Return the weave pattern as a list of lists. Each pair of virtual 
        helices has a list associated with it. Each of those lists contains the
        interhelical distance as a function of base-pair index along the virtual
        helix.
        
        This method copes with insertions/deletions by creating a list of base-
        base midpoints for each virtual helix, and assuming that each pair of 
        virtual helices has the same number of base-base midpoints. This should
        mean that the base pairs are fairly well lined up and results in a 
        reasonable definition of the weave pattern.
        """
        out = [[] for xx in range(len(self.vhelix_indices)-1)]
        for vhi in range(len(self.vhelix_indices)-1):
            # make a list of base-base midpoints for each helix
            bbms1 = self.get_vh_midpoints(self.vhelix_indices[vhi])
            bbms2 = self.get_vh_midpoints(self.vhelix_indices[vhi + 1])

            if len(bbms1) != len(bbms2):
                if verbose:
                    print >> sys.stderr, "INFO: skipping virtual helix pair %d, %d: number of base pairs in virtual helix number %d (%d) different to number of base pairs in virtual helix number %d (%d)" % (self.vhelix_indices[vhi], self.vhelix_indices[vhi+1], self.vhelix_indices[vhi], len(bbms1), self.vhelix_indices[vhi+1], len(bbms2))
            else:
                # compute the weave as the inter-base-base-midpoint distance
                out[vhi] = [-1 for xx in bbms1]
                for ii in range(len(bbms1)):
                    displ = min_distance(bbms2[ii], bbms1[ii], self._sys._box)
                    out[vhi][ii] = np.sqrt(np.dot(displ, displ))

        return out

    def compute_vh_midpoints(self):
        """
        compute the base pair midpoints for each virtual helix just once to 
        save time
        """
        if self._sys:
            self.vh_midpoints = []
            for vhi in range(len(self.vhelix_indices)):
                self.vh_midpoints.append(self.get_vh_midpoints(vhi))

    def get_vh_midpoints(self, vhi):
        """
        Return a list of base pair midpoints for a virtual helix. Effectively 
        'flattens' a virtual helix, dealing with any skips or loops.
        """
        bbms = []
        for vb in self.vvib[vhi]:
            bbms.extend([self.get_bb_midpoint(nuc, pbc = True) for nuc in self.get_nucleotides(self.vhelix_indices[vhi], vb, "default")])
        return bbms

    def get_alignment(self, trim, period):
        # find alignment of a pair of nucleotides, ie to see if ones that are supposed to line up do.
        # assumes no skip/loop / insertion/deletion
        pair_count = len(self.vvib[0]) - trim * 2 - period
        twists = range(pair_count)
        for ii in range(pair_count):
            ss = self._sys
            n1 = self.get_nucleotides(self.vhelix_indices[0], self.vvib[0][trim+ii])[0]
            n2 = self.complementary_list[n1]
            n3 = self.get_nucleotides(self.vhelix_indices[0], self.vvib[0][trim+period+ii])[0]
            n4 = self.complementary_list[n3]
            r1 = ss._nucleotides[n1].get_pos_base()
            r2 = ss._nucleotides[n2].get_pos_base()
            r3 = ss._nucleotides[n3].get_pos_base()
            r4 = ss._nucleotides[n4].get_pos_base()

            ##
            mid12 = (r1+r2)/2
            mid34 = (r3+r4)/2
            n = mid34 - mid12
            a11 = ss._nucleotides[n1]._a1
            a13 = ss._nucleotides[n3]._a1
            n = norm(n)
            a11 = norm(a11)
            a13 = norm(a13)
            # find the component of the a1 vectors that is in the plane with normal n
            a11prime = a11 - np.dot(a11, n) * n
            a13prime = a13 - np.dot(a13, n) * n
            a11prime = norm(a11prime)
            a13prime = norm(a13prime)

            twist = np.arccos(np.dot(a11prime, a13prime))
            twist *= 180./np.pi
            xp = np.cross(a11prime, a13prime)
            if np.dot(xp,n) < 0:
                twist *= -1

            twists[ii] = twist

        return twists

    def get_hj_alignment(self, start, period):
        # find alignment of a pair of nucleotides, ie to see if ones that are supposed to line up do.
        # assumes exactly 2 virtual helices
        # assumes no skip/loop / insertion/deletion
        # start should be virtual base index of the first crossover (reading left to right)
        ss = self._sys
        twists = range(2)
        for ii in (0,1):
            # when ii = 0
            #                period
            #               <------->
            # vhelix 1 ====2===  ===4=====
            #              ||       ||
            # vhelix 0 ====1========3=====
            #
            # or when ii = 1
            # vhelix 1 =====2===  ===4====
            #              ||       ||
            # vhelix 0 =====1========3====
            n1 = self.get_nucleotides(self.vhelix_indices[0], start+ii)[0]
            n2 = self.get_nucleotides(self.vhelix_indices[1], start+ii)[0]
            n3 = self.get_nucleotides(self.vhelix_indices[0], start+period+ii)[0]
            n4 = self.get_nucleotides(self.vhelix_indices[1], start+period+ii)[0]

            bbm1 = self.get_bb_midpoint(n1)
            bbm2 = self.get_bb_midpoint(n2)
            bbm3 = self.get_bb_midpoint(n3)
            bbm4 = self.get_bb_midpoint(n4)

            ##
            n = bbm3 - bbm1
            align12 = bbm2 - bbm1
            align34 = bbm4 - bbm3
            n = norm(n)
            align12 = norm(align12)
            align34 = norm(align34)
            # find the component of the align vectors that is in the plane with normal n
            align12prime = align12 - np.dot(align12, n) * n
            align34prime = align34 - np.dot(align34, n) * n
            align12prime = norm(align12prime)
            align34prime = norm(align34prime)


            twist = np.arccos(np.dot(align12prime, align34prime))
            twist *= 180./np.pi
            xp = np.cross(align12prime, align34prime)
            if np.dot(xp,n) < 0:
                twist *= -1
            twists[ii] = twist

        return twists

    def get_arm_vec(self, vhc,vbc,vbn):
        system = self._sys
        nuc0scid = self.get_nucleotides(vhc,vbc)[0]
        nuc0stid = self.complementary_list[nuc0scid]
        nuc1scid = self.get_nucleotides(vhc,vbn)[0]
        nuc1stid = self.complementary_list[nuc1scid]

        nuc0sc = system._nucleotides[nuc0scid]
        nuc0st = system._nucleotides[nuc0stid]
        nuc1sc = system._nucleotides[nuc1scid]
        nuc1st = system._nucleotides[nuc1stid]

        v0 = nuc0sc.get_pos_base() + nuc0st.get_pos_base()
        v1 = nuc1sc.get_pos_base() + nuc1st.get_pos_base()

        v = norm(v1-v0)
        return v

    def get_junction_normal(self, hj, armlen=3):
        # hj should come from hjs[ii] where hjs=origami.get_holliday_junctions()
        # assume double crossover, no skip/loop (not sure if matters)
        vh1 = hj[0]
        vb1 = hj[1]
        vh2 = hj[2]
        vb2 = hj[3]
        # arm vectors - in cadnano representation:
        # A====>C
        #    ||
        # D<====B

        vA = self.get_arm_vec(vh1,vb1,vb1-armlen)
        vB = self.get_arm_vec(vh2,vb2,vb2+armlen)
        vC = self.get_arm_vec(vh1,vb2,vb2+armlen)
        vD = self.get_arm_vec(vh2,vb1,vb1-armlen)

        n = np.cross(vA,vB) + np.cross(vB,vC) + np.cross(vC,vD) + np.cross(vD,vA)
        n = norm(n)

        return n

    def get_junction_normal2(self, hj):
        # hj should come from hjs[ii] where hjs=origami.get_holliday_junctions()
        # assume double crossover
        vh1 = hj[0]
        vb1 = hj[1]
        vh2 = hj[2]
        vb2 = hj[3]
        # arm vectors - in cadnano representation:
        # A====>C
        #    ||
        # D<====B

        ss = self._sys
        n1_id = self.get_nucleotides(vh2, vb1)[0]
        n2_id = self.get_nucleotides(vh1, vb1)[0]
        n3_id = self.get_nucleotides(vh1, vb2)[0]

        r1 = ss._nucleotides[n1_id].get_pos_base()
        r2 = ss._nucleotides[n2_id].get_pos_base()
        r3 = ss._nucleotides[n3_id].get_pos_base()

        v1 = r2 - r1
        v2 = r2 - r3
        n = np.cross(v1,v2)
        n = norm(n)

        return n

    def intra_helix_autocorr(self):
        """
        Get the autocorrelation (i.e. the dot product) between the intra-helix 
        vectors along each pair of virtual helices in the origami.

        This method does assume that the base pairs on adjacent virtual helices
        'line up,' in the sense that the list index of the closest base pair 
        on one helix is the same as the index of the closest base pair on the
        neighbouring helix. A weak check for this is the condition that the
        lengths of the bbm lists (and so the number of base pairs in each virtual
        helix) are the same.

        Note that the list this method returns is indexed by the actual virtual
        helix number given by cadnano, not the 'virtual helix index' which
        indexes the self.vhelix_indices list.

        run compute_vh_midpoints() first
        """
        corrs = [[] for xx in range(max(self.vhelix_indices))]
        for vhi in range(len(self.vhelix_indices) - 1):
            try:
                bbms1 = self.vh_midpoints[vhi]
            except IndexError:
                print >> sys.stderr, "make sure Origami.compute_vh_midpoints() has run"
                raise
            bbms2 = self.vh_midpoints[vhi + 1]
            # weak check for base-pair aligment across virtual helices
            if len(bbms1) != len(bbms2):
                continue
            nbp = len(bbms1)
            corr = [0 for xx in range(nbp)]
            count = [0 for xx in range(nbp)]
            for bpidi in range(nbp):
                intra_i = norm(self.min_distance(bbms2[bpidi], bbms1[bpidi]))
                for bpidj in range(bpidi, nbp):
                    intra_j = norm(self.min_distance(bbms2[bpidj], bbms1[bpidj]))
                    corr[bpidj - bpidi] += np.dot(intra_i, intra_j)
                    count[bpidj - bpidi] += 1

            corrs[self.vhelix_indices[vhi]] = [float(corr[ii])/count[ii] for ii in range(nbp)]

        return corrs

    def corrugation(self, hj, span=16):
        """
        Find the corrugation a distance span around the junction hj. The 
        corrugation pattern is quantified as the angle from the dot product between the average
        of the two intra-helix vectors at the junction and each of the other
        intra-helix vectors.

        run compute_vh_midpoints() first
        """
        vh1 = hj[0]
        vb1 = hj[1]
        vh2 = hj[2]
        vb2 = hj[3]

        # label the base pairs like this
        # ===A-B==>
        #    | | 
        # <==C-D===
        bbmA = self.get_bb_midpoint(self.get_nucleotides(vh1, vb1)[0], pbc=True)
        bbmB = self.get_bb_midpoint(self.get_nucleotides(vh1, vb2)[0], pbc=True)
        bbmC = self.get_bb_midpoint(self.get_nucleotides(vh2, vb1)[0], pbc=True)
        bbmD = self.get_bb_midpoint(self.get_nucleotides(vh2, vb2)[0], pbc=True)

        # average of the helix axes at the junction
        av_axis = norm(norm(self.min_distance(bbmA, bbmB)) + norm(self.min_distance(bbmC, bbmD)))

        # find component of reference intrahelix vector that is in the plane 
        # normal to av_axis
        ref_intra = self.min_distance(bbmC, bbmA) + self.min_distance(bbmD, bbmB)
        ref_intra_in_plane = norm(ref_intra - np.dot(ref_intra, av_axis) * av_axis)

        # bbmsX are lists of base-pair midpoints for each virtual helix
        bbms1 = self.vh_midpoints[vh1]
        bbms2 = self.vh_midpoints[vh2]
        # find the location of the junction within the base-pair midpoint lists
        idA = array_index(bbms1, bbmA)
        idB = array_index(bbms1, bbmB)
        idC = array_index(bbms2, bbmC)
        idD = array_index(bbms2, bbmD)

        # now loop through the revelant intra-strand pairs and compute the dot
        # products of each intra-strand vector with the reference vector
        ret = []
        for xpos in range(-span, 0):
            ret.append(self.corrug_angle(xpos, bbms1, bbms2, idA, idC, av_axis, ref_intra_in_plane))
        for xpos in range(1, span+1):
            ret.append(self.corrug_angle(xpos, bbms1, bbms2, idB, idD, av_axis, ref_intra_in_plane))

        return ret

    def corrug_angle(self, xpos, bbms1, bbms2, id1, id2, av_axis, ref_intra_in_plane):
        """
        compute the corrugation angle as the dot product of the intra-strand
        vector at xpos with the reference vector, after both have been projected
        into the plane perpendicular to the average helix axis at the junction.
        """
        # the if condition is here to ensure we don't have a negative 
        # index (which with python logic gets treated as the penultimate
        # list element) and that we don't run off the end of the list
        if id1 + xpos > 0 and id1 + xpos < len(bbms1) and id2 + xpos > 0 and id2 + xpos < len(bbms2):
            this_intra = self.min_distance(bbms2[id2 + xpos], bbms1[id1 + xpos])
            this_intra_in_plane = norm(this_intra - np.dot(this_intra, av_axis) * av_axis)
            # give the dot product a sign corresponding to the sense of the
            # angle between the vectors, given by the triple product between
            # the two intra-helix vectors and the average helix axis.
            fac = 1
            if np.dot(np.cross(this_intra_in_plane, ref_intra_in_plane), av_axis) < 0:
                fac = -1
            return fac * np.arccos(np.dot(this_intra_in_plane, ref_intra_in_plane))
        else:
            return "NA"

    def old_get_corrug_weave(self, vh1, vh2, vb, n, norm2, corrug_signed=False):
        # find distance between adjacent helices and decompose into weave and corrugation
        # for a pair of base pairs located at (vh1, vb) and (vh2, vb)
        # assumes no skip/loop
        ss = self._sys
        n1 = self.get_nucleotides(vh1, vb)[0]
        n2 = self.complementary_list[n1]
        n3 = self.get_nucleotides(vh2, vb)[0]
        n4 = self.complementary_list[n3]

        r1 = ss._nucleotides[n1].get_pos_base()
        r2 = ss._nucleotides[n2].get_pos_base()
        r3 = ss._nucleotides[n3].get_pos_base()
        r4 = ss._nucleotides[n4].get_pos_base()
        mid12 = (r1+r2)/2
        mid34 = (r3+r4)/2
        d = mid12 - mid34

        # caution - wd is scalar, d and cd are vectors here
        wd = np.dot(n, d) # weave_distance, corrugation_distance
        cd = d - wd*n

        # corrugation defined negative if opposite dir to plane normal
        # ^ (code deleted) not meaningful with corrugation defined as component of distance between dhds normal to local plane, rather than distance from local plane, as it currently is

        if corrug_signed:
            if np.dot(cd,norm2) > 0:
                fac = 1
            else:
                fac = -1
        else:
            fac = 1

        cd = np.sqrt(np.dot(cd,cd))

        return wd, cd*fac, d

    def get_nearest_native_bonded(self, vhi, vb, increment, bound):
        """
        Return nearest natively bonded base-pair.

        args:
        vhi: vhelix_indices index
        vb: virtual base number (i.e. the actual number corresponding to that virtual base in cadnano)
        increment: must be +1 or -1, the amount to increment the virtual base by when searching for a virtual base with a native bond
        bound: if we get to this virtual base and we still haven't found anything, we give up and throw and error

        example usage:
        get_nearest_native_bonded(vhi, vh_begin[vhi], 1, vh_end[vhi])
        """
        assert increment == 1 or increment == -1
        n1, n2 = self.get_nucleotides(self.vhelix_indices[vhi], vb, "double")[:2]
        #print n1, n2, self.interaction_list[n1], self.interaction_list[n2]
        while n2 not in self.interaction_list[n1]:
            try:
                # we can't have the vector starting at the same place it ends or after its end
                if vb == bound:
                    raise IndexError
                n1, n2 = self.get_nucleotides(self.vhelix_indices[vhi], vb, "double")[:2]
                #print n1, n2, self.interaction_list[n1], self.interaction_list[n2]
                vb += increment
            except IndexError:
                print >> sys.stderr, "This error is likely to be caused by not running get_h_bond_list first, or (much less likely) because there are no native bonds for this virtual helix."
                print vhi, vb
                raise
        #print "nnb", vhi, vb, n1, n2
        return n1, n2
    
    def get_av_helix_axis(self, vhelix_extent = False, discard_unbonded = True):
        """
        find the average helix axis. Optionally discard unbonded (but designed) 'base-pairs', and ignore the regions outside the 'vhelix_extent'
        """
        av_helix_axis = np.zeros(3)

        # set the extent to be considered of each virtual helix
        if vhelix_extent:
            vh_begin = vhelix_extent[0]
            vh_end = vhelix_extent[1]
        else:
            vh_begin = [self.vvib[vhi][0] for vhi in range(len(self.vhelix_indices))]
            vh_end = [self.vvib[vhi][-1] for vhi in range(len(self.vhelix_indices))]

        for vhi in range(len(self.vhelix_indices)):
            n1, n2 = self.get_nearest_native_bonded(vhi, vh_begin[vhi], +1, vh_end[vhi])
            n3, n4 = self.get_nearest_native_bonded(vhi, vh_end[vhi], -1, vh_begin[vhi])

            bbm1 = self.get_bb_midpoint(n1, pbc = True)
            bbm2 = self.get_bb_midpoint(n3, pbc = True)

            # find the vector between the two base-pairs
            av_helix_axis += min_distance(bbm1, bbm2, self._sys._box)
            
        av_helix_axis /= len(self.vhelix_indices)
        return av_helix_axis
    
    def get_3d_twist(self, f_edge_vh, vhelix_extent = False, discard_unbonded = True):
        """
        return the 3d global twist (top-bottom and left-right)

        RUN get_h_bond_list FIRST!!
        """

        ## Bring first strand to centre...!
        # save a copy of the system
        sys_copy = base.System(self._sys._box)
        sys_copy = sys_copy.join(self._sys)

        to_translate = -self._sys._strands[0].get_cm_pos()
        for strandid, strand in enumerate(self._sys._strands):
            strand.translate(to_translate)

        av_helix_axis = self.get_av_helix_axis(vhelix_extent, discard_unbonded)
        
        vhelices1, vhelices2 = self.parse_vh_edge(f_edge_vh)

        if vhelix_extent:
            vh_begin = vhelix_extent[0]
            vh_end = vhelix_extent[1]
        else:
            vh_begin = [self.vvib[vhi][0] for vhi in range(len(self.vhelix_indices))]
            vh_end = [self.vvib[vhi][-1] for vhi in range(len(self.vhelix_indices))]

        top = np.zeros(3)
        bottom = np.zeros(3)
        top_left = [np.zeros(3) for xx in range(len(vhelices1))]
        top_right = [np.zeros(3) for xx in range(len(vhelices1))]
        bottom_left = [np.zeros(3) for xx in range(len(vhelices2))]
        bottom_right = [np.zeros(3) for xx in range(len(vhelices2))]
        # top
        for id, vh in enumerate(vhelices1):
            vhi = self.vhelix_indices.index(vh)
            # should work with loops
            n1, n2 = self.get_nearest_native_bonded(vhi, vh_begin[vhi], 1, vh_end[vhi])
            n3, n4 = self.get_nearest_native_bonded(vhi, vh_end[vhi], -1, vh_begin[vhi])

            vec1 = self.get_bb_midpoint(n1, pbc = True)
            vec2 = self.get_bb_midpoint(n3, pbc = True)
            #print self._sys._nucleotides[n1].get_pos_base(), self._sys._nucleotides[self.complementary_list[n1]].get_pos_base(), vec1


            # left and right correspond to left (low virtual base index) and right (high virtual base index) in cadnano
            top_left[id] = vec1
            #print "setting top_right[", id, "] to", vec2, "and top_left to", vec1
            top_right[id] = vec2
            top += min_distance(vec1, vec2, self._sys._box)
            #top += vec2-vec1
        # bottom
        for id, vh in enumerate(vhelices2):
            vhi = self.vhelix_indices.index(vh)
            # should work with loops
            n1, n2 = self.get_nearest_native_bonded(vhi, vh_begin[vhi], 1, vh_end[vhi])
            n3, n4 = self.get_nearest_native_bonded(vhi, vh_end[vhi], -1, vh_begin[vhi])

            vec1 = self.get_bb_midpoint(n1, pbc = True)
            vec2 = self.get_bb_midpoint(n3, pbc = True)

            bottom_left[id] = vec1
            #print "setting bottom_right[", id, "] to", vec2, "and bottom_left to", vec1
            bottom_right[id] = vec2
            bottom += min_distance(vec1, vec2, self._sys._box)
            #bottom += vec2-vec1

        top /= len(vhelices1)
        bottom /= len(vhelices2)

        left = np.zeros(3)
        for ii in range(len(bottom_left)):
            left += min_distance(bottom_left[ii], top_left[ii], self._sys._box)
            #left += top_left[ii] - bottom_left[ii]
        left /= len(bottom_left)

        right = np.zeros(3)
        for ii in range(len(bottom_right)):
            right += min_distance(bottom_right[ii], top_right[ii], self._sys._box)
            #right += top_right[ii] - bottom_right[ii]
            #print "TL,BL", top_left[ii], bottom_left[ii]
            #print "TLBL", top_left[ii] - bottom_left[ii]
            #print "TR,BR", top_right[ii], bottom_right[ii]
            #print "TRBR", top_right[ii] - bottom_right[ii]
        right /= len(bottom_right)

        #print left, right

        # Schematic of vectors in cadnano representation. Each ===> arrow is a double helix in the cadnano representation
        #
        #       top   ------------>
        #
        #  ^          ============>   ^
        #  | left     <============   | right
        #  |          ============>   |
        #  |          <============   |
        #
        #      bottom ------------>

        twist1 = dihedral_angle_sense(left,right,av_helix_axis,sanity_check=True)
        if twist1 == False:
            print >> sys.stderr, "\n\nsanity check for left/right twist failed"
            print >> sys.stderr, "left is %f %f %f" % (left[0], left[1], left[2])
            print >> sys.stderr, "right is %f %f %f" % (right[0], right[1], right[2])
            print >> sys.stderr, "axis is %f %f %f" % (av_helix_axis[0], av_helix_axis[1], av_helix_axis[2])
            print >> sys.stderr, "projections onto axis are %f %f" % (np.dot(norm(left),norm(av_helix_axis)), np.dot(norm(right),norm(av_helix_axis)))
            print >> sys.stderr, "dihedral angle is %f\n\n" % (dihedral_angle_sense(left,right,av_helix_axis) * 180./np.pi)
            raise RuntimeError
        twist2 = dihedral_angle_sense(top,bottom,left+right,sanity_check=True)
        if twist2 == False:
            print >> sys.stderr, "sanity check for top/bottom twist failed"
            print >> sys.stderr, "left is %f %f %f" % (left[0], left[1], left[2])
            print >> sys.stderr, "right is %f %f %f" % (right[0], right[1], right[2])
            print >> sys.stderr, "axis is %f %f %f" % (left[0]+right[0], left[1]+right[1], left[2]+right[2])
            print >> sys.stderr, "dihedral angle is %f" % abs(dihedral_angle_sense(top,bottom,left+right))
            raise RuntimeError
        twist1 *= 180./np.pi
        twist2 *= 180./np.pi

        # restore the copy of the system
        self._sys = sys_copy

        return twist1, twist2

    def is_strand_end(self, vh, vb, strand_type):
        """
        return string first or last for first or last nucleotide of the strand, otherwise false
        this assumes we don't have any 1-nucleotide-long strands....
        """

        if strand_type == "scaffold":
            strandid, nucids = self._cad2cudadna._scaf[(vh, vb)]
        else:
            assert strand_type == "staple"
            strandid, nucids = self._cad2cudadna._stap[(vh, vb)]

        strand = self._sys._strands[strandid]
        for nucid in nucids:
            # I hope that a single virtual base never contains both the beginning and end of a strand...
            if nucid == strand._first:
                return "first"
            elif nucid == strand._last:
                return "last"
            else:
                return False

    def get_n_nicks(self):
        """
        return the number of nicks (staple or scaffold) in an origami
        Here a nick must be between adjacent strands (e.g. if a nucleotide is missing that isn't counted as a nick at all)
        """
        nicks = self.get_nicks()
        return len(nicks[0]) + len(nicks[1])
    
    def get_nicks(self):
        """
        return a list [scaffold_nicks, staple_nicks], each of which is a list of nicked neighbours in the origami
        Here a nick must be between adjacent strands (e.g. if a nucleotide is missing that isn't counted as a nick at all)
        An element of each nick list has virtual helix number (NOT array index), virtual base 1, virtual base 2
        """
        scaf_nicks = []
        stap_nicks = []
        for vhi in range(len(self.vhelix_indices)):
            vh = self.vhelix_indices[vhi]
            for vb in self.vvib[vhi]:
                # check we have another virtual base after this one
                if vb+1 in self.vvib[vhi]:
                    # check for scaffold nick
                    if (self.is_strand_end(vh, vb, "scaffold") == "first" and self.is_strand_end(vh, vb + 1, "scaffold") == "last") or (self.is_strand_end(vh, vb, "scaffold") == "last" and self.is_strand_end(vh, vb + 1, "scaffold") == "first"):
                        scaf_nicks.append([vh, vb, vb + 1])
                    # check for staple nick
                    if (self.is_strand_end(vh, vb, "staple") == "first" and self.is_strand_end(vh, vb + 1, "staple") == "last") or (self.is_strand_end(vh, vb, "staple") == "last" and self.is_strand_end(vh, vb + 1, "staple") == "first"):
                        stap_nicks.append([vh, vb, vb + 1])
                        
        return scaf_nicks, stap_nicks

    def is_nicked(self, vh, vb1, vb2):
        """
        return true if there is a nick across vb1 and vb2, false otherwise
        """
        # check for scaffold nick
        if (self.is_strand_end(vh, vb1, "scaffold") == "first" and self.is_strand_end(vh, vb2, "scaffold") == "last") or (self.is_strand_end(vh, vb1, "scaffold") == "last" and self.is_strand_end(vh, vb2, "scaffold") == "first"):
            return True
        
        # check for staple nick
        if (self.is_strand_end(vh, vb1, "staple") == "first" and self.is_strand_end(vh, vb2, "staple") == "last") or (self.is_strand_end(vh, vb1, "staple") == "last" and self.is_strand_end(vh, vb2, "staple") == "first"):
            return True
            
        return False

    def check_branch_migration_strict(self, vh1, vb1, vh2, vb2):
        """
        check for branch migration using strict definition; return true if there is branch migration, false otherwise

        strict definition: at least one of the 4 junction base pairs is missing
        self.interaction_list must be filled, using some version of get_h_bond_list
        the two virtual helices and two virtual bases that define the junction's position are given as the four arguments (the actual numbers should be given, not the indices in their respective self.xxx arrays)
        """
        assert len(self.interaction_list) > 0
        for square in ([vh1, vb1], [vh1, vb2], [vh2, vb1], [vh2, vb2]):
            (nuc1, nuc2) = self.get_nucleotides(square[0], square[1], type="double")
            if nuc2 not in self.interaction_list[nuc1]:
                return True
        return False

    def has_full_hbond(self, vh, vb):
        """
        return true if all cadnano base-pairs at this virtual base have their correct base-pairs in the configuration

        assumes self.interaction_list is already filled
        vh: virtual helix number that appears in cadnano
        vb: virtual base number that appears in cadnano
        """
        assert len(self.interaction_list) > 0
        nucs = self.get_nucleotides(vh, vb, "single")
        for nuc in nucs:
            nuc2 = self.complementary_list[nuc]
            if nuc2 not in self.interaction_list[nuc]:
                # return false if any base-pairs are missing
                return False

        # return true if we got to the end
        return True

    def is_junction(self, vh, vb1, vb2, hjs_double, hjs_single):
        """
        return 'double' if double crossover, 'single' if single crossover, False otherwise

        vh: virtual helix number that appears in cadnano
        vb1: virtual base number that appears in cadnano (left square)
        vb2: virtual base number that appears in cadnano (right square)
        hjs_double: array of holliday junctions (double crossovers)
        hjs_single: array of single crossovers
        """
        # test for double crossover: 2 possible pairs of cadnano squares
        vh1s = [[xx[0], xx[1], xx[3]] for xx in hjs_double]
        vh2s = [[xx[1], xx[1], xx[3]] for xx in hjs_double]
        if [vh, vb1, vb2] in vh1s or [vh, vb1, vb2] in vh2s:
            return "double"

        # test for single crossover: 2 possible cadnano squares (in this case vb2 isn't used)
        vhvb1 = [[xx[0], xx[1]] for xx in hjs_single]
        vhvb2 = [[xx[2], xx[1]] for xx in hjs_single]
        if [vh, vb1] in vhvb1 or [vh, vb1] in vhvb2:
            return "single"

        return False

    def parse_vhelix_extent_file(self, vhelix_extentf):
        """
        parse the vhelix_extent file and return a list [vh_begin, vh_end], where vh_begin is a list of length
        <number of virtual helices> and each element being the left virtual base, and analogous for vh_end.

        ARGS
        vhelix_extentf: path to vhelix_extent file

        vhelix_extent file example:
        10, 30

        This function could be extended to parse a file with different vh_begin and vh_end for each virtual helix;
        the machinery that uses this function should be able to deal with that without modification.
        """
        vhelix_extent = []
        try:
            f = open(vhelix_extentf, "rb")
        except IOError:
            base.Logger.die("could not find file %s, dying" % vhelix_extentf)

        for str in f.readline().replace(" ", "").split(","):
            i = int(str)
            vhelix_extent.append(i)
            
        assert len(vhelix_extent) == 2

        vh_begin_raw = vhelix_extent[0]
        vh_end_raw = vhelix_extent[1]
        vh_begin = [vh_begin_raw for x in self.vvib]
        vh_end = [vh_end_raw for x in self.vvib]
        return (vh_begin, vh_end)

    def parse_vh_edge(self, f_edge_vh):
        """
        return a tuple of lists, each element corresponds to a vhelix number that is either in the top (vhelices1) or bottom (vhelices2) of the origami
        """
        try:
            f = open(f_edge_vh, "r")
        except IOError:
            base.Logger.die("could not find file %s, dying" % f_edge_vh)

        str1, str2 = f.readline().replace(" ","").split(";")
        vhelices1 = [int(x) for x in str1.split(",")]
        vhelices2 = [int(x) for x in str2.split(",")]

        return vhelices1, vhelices2

    def get_vh_spline(self, vh, mode = "midpoint", vhelix_extent = False, discard_unbonded = True, force_circular = False):
        """
        return a cartesian spline for the base-base midpoints of a virtual helix

        args:
        force_circular: force a closed curve; the ends of the midpoint curve (which is open) will be joined together by a straight line
        vh: virtual helix number (NOT the corresponding index from the vhelix_indices array!)
        """
        import scipy.interpolate

        vhi = self.vhelix_indices.index(vh)
        if vhelix_extent:
            vh_begin = vhelix_extent[0]
            vh_end = vhelix_extent[1]
        else:
            vh_begin = [self.vvib[vhi][0] for vhi in range(len(self.vhelix_indices))]
            vh_end = [self.vvib[vhi][-1] for vhi in range(len(self.vhelix_indices))]

        if discard_unbonded:
            assert self.interaction_list[0] != -1, "If using discard_unbonded, make sure that you run Origami.get_h_bond_list() or similar first to fill in the interaction list"
            
        if mode == "midpoint":
            bbms = []
            for vb in self.vvib[vhi]:
                nucs = self.get_nucleotides(vh, vb, "default")
                for nuc in nucs:
                    # check for correct base pairing if necessary
                    if (discard_unbonded and self.complementary_list[nuc] in self.interaction_list[nuc]) or not discard_unbonded:
                        bbms.append(self.get_bb_midpoint(nuc))
            if force_circular:
                bbms.append(bbms[0])
            ret = vecs2spline(bbms, force_circular)
        else:
            assert mode == "both strands"
            bbms1 = []
            bbms2 = []
            for vb in self.vvib[vhi]:
                nucs1 = self.get_nucleotides(vh, vb, "scaf")
                nucs2 = self.get_nucleotides(vh, vb, "stap")
                for nuc in nucs1:
                    # check for correct base pairing if necessary
                    if (discard_unbonded and self.complementary_list[nuc] in self.interaction_list[nuc]) or not discard_unbonded:
                        bbms1.append(self._sys._nucleotides[nuc].get_pos_base())
                for nuc in nucs2:
                    # check for correct base pairing if necessary
                    if (discard_unbonded and self.complementary_list[nuc] in self.interaction_list[nuc]) or not discard_unbonded:
                        bbms2.append(self._sys._nucleotides[nuc].get_pos_base())
            if force_circular:
                bbms1.append(bbms1[0])
                bbms2.append(bbms2[0])
            ret = [vecs2spline(bbms1, force_circular), vecs2spline(bbms2, force_circular)]

        return ret

    def get_n_bp(self, vh):
        """
        Return the number of base pairs in the virtual helix with number 
        'vh'
        """
        nucs = []
        for vb in self.vvib[self.vhelix_indices.index(vh)]:
            nucs.extend(self.get_nucleotides(vh, vb, "default"))
        return len(nucs)

    def min_distance (self, r1, r2):
        """
        return the minimum image distance in going from r1 to r2, in this system's box
        """
        return min_distance(r1, r2, self._sys._box)
