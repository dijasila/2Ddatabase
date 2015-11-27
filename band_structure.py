from math import sqrt, pi
import numpy as np
from ase.units import Hartree, m, kJ, _hbar, _me, Bohr
from gpaw import GPAW
from gpaw.kpt_descriptor import to1bz
from fractions import gcd

def get_ibz_kpt_path(calc, wpts_xc):
    kd = calc.wfs.kd
    op_scc = kd.symmetry.op_scc
    cell_cv = calc.atoms.get_cell()
    rcell_cv = 2. * pi * calc.atoms.get_reciprocal_cell()
    N_c = calc.wfs.kd.N_c
    ibzk_kc = calc.get_ibz_k_points()

    x_x = []
    k_xc = []
    k_x = []
    x = 0.
    X = []
    for nwpt in range(1, len(wpts_xc)):
        X.append(x)
        to_c = wpts_xc[nwpt]
        from_c = wpts_xc[nwpt - 1]
        vec_c = to_c - from_c
        Nv_c = (vec_c * N_c).astype(int)
        Nv = abs(gcd(gcd(Nv_c[0], Nv_c[1]), Nv_c[2]))
        dv_c = vec_c / Nv
        dv_v = np.dot(dv_c, rcell_cv)
        dx = np.linalg.norm(dv_v)
        if nwpt == len(wpts_xc) - 1:
            X.append(x + Nv * dx)
            Nv += 1
        for n in range(Nv):
            k_c = from_c + n * dv_c
            bzk_c = to1bz(np.array([k_c]), cell_cv)[0]
            ibzkpt, _, _, _ = kd.find_ibzkpt(op_scc, ibzk_kc, bzk_c)
            x_x.append(x)
            k_xc.append(k_c)
            k_x.append(ibzkpt)
            x += dx
    return np.array(x_x), np.array(X), np.array(k_xc), np.array(k_x)

def get_kpt_path(calc, wpts_xc):
    """
    Gets the bz k-points included in the calculator along a path between a given
    set of reciprocal space waypoints (typically high symmetry points). Useful for
    extracting a band  structure from a calculation with only Monkhorst-Pack
    k-point sampling.

    Keyword arguments:
    calc: The PAW calculator object.
    wpts_n: A list of indices of the bz k-points to use as waypoints.
    """

    rcell_cv = 2 * np.pi * calc.atoms.get_reciprocal_cell()
    wpts_xv = np.dot(wpts_xc, rcell_cv)
    bz_kpts = calc.get_bz_k_points()
    bz_rkpts = np.dot(bz_kpts, rcell_cv)

    pts = []
    x = 0.

    path = []

    x = 0
    # Loop through each of the subpaths
    for nwpt in range(1, len(wpts_xc)):
        to_c = wpts_xc[nwpt]
        frm_c = wpts_xc[nwpt-1]
        to_v = np.dot(to_c, rcell_cv)
        frm_v = np.dot(frm_c, rcell_cv)
        vec_v = to_v - frm_v
        subpath = []
        for n, k_v in enumerate(bz_rkpts):
            kvec_v = k_v - frm_v
            knorm = np.linalg.norm(kvec_v)
            if knorm < 1e-9:
                subpath.append((n, x + knorm))
            else:
                cos = np.dot(vec_v, kvec_v) / np.linalg.norm(vec_v) / \
                  np.sqrt(np.sum(kvec_v**2))
                # Check if the point is on the line through to and frm
                if (abs(cos - 1) < 1e-9 and
                    np.linalg.norm(kvec_v) < np.linalg.norm(vec_v) + 1e-9):
                    subpath.append((n, x + knorm))

        # Add next waypoint to path
        subpath.sort(key=lambda elem: elem[1])
        if len(path) > 0 and subpath[0][0] == path[-1][0]:
            path += subpath[1:]
        else:
            path += subpath
        x += np.linalg.norm(vec_v)

    return zip(*path)

def get_ibz_values(calc, kpts_kc, values_k):
    ibzk_kc = calc.get_ibz_k_points()
    ibzvalues_k = [None] * len(ibzk_kc)
    kibzk_k = []

    cell_cv = calc.wfs.gd.cell_cv
    kd = calc.wfs.kd
    op_scc = kd.symmetry.op_scc
    for k, k_c in enumerate(kpts_kc):
        bz1k_c = to1bz(np.array([k_c]), cell_cv)[0]
        ibzkpt, _, _, _ = kd.find_ibzkpt(op_scc, ibzk_kc, bz1k_c)
        if not ibzkpt in kibzk_k:
            kibzk_k.append(ibzkpt)
            ibzvalues_k[ibzkpt] = values_k[k]
    
    return ibzvalues_k

def get_homo_lumo_states(calc, epsilon=1e-2):
    """
    Finds the homo and lumo states of a calculation and return their energies
    as well as their spin, k-point and band indices.

    calc: Calculator object or string representing the path to a .gpw file.
    epsilon: Occupation threshold. Minimum occupation to consider as occupied.
    """
    if isinstance(calc, (str, unicode)):
        calc = GPAW(calc, txt=None)

    homo = (float('-inf'), 0, 0, 0)
    lumo = (float('inf'), 0, 0, 0)
    for kpt in calc.wfs.kpt_u:
        occ_n = np.amax(np.where(kpt.f_n/kpt.weight > epsilon)[0])
        unocc_n = np.amin(np.where(kpt.f_n/kpt.weight < epsilon)[0])
        eps_occ = kpt.eps_n[occ_n] * Hartree
        eps_unocc = kpt.eps_n[unocc_n] * Hartree
        if eps_occ > homo[0]:
            homo = (eps_occ, kpt.s, kpt.k, occ_n)
        if eps_unocc < lumo[0]:
            lumo = (eps_unocc, kpt.s, kpt.k, unocc_n)

    return homo, lumo

def get_band_gap(calc, epsilon=1e-2):
    """
    Calculates the band gap. If the system is a metal 0 is returned.

    calc: Calculator object or string representing the path to a .gpw file.
    epsilon: Occupation threshold. Minimum occupation to consider as occupied.
    """
    if isinstance(calc, (str, unicode)):
        calc = GPAW(calc, txt=None)

    homo, lumo = get_homo_lumo_states(calc, epsilon)
    val_s = homo[1]
    val_n = homo[3]
    cond_s = lumo[1]
    cond_n = lumo[3]

    metal = False
    occ_val = np.array([kpt.f_n[val_n]/kpt.weight for kpt in calc.wfs.kpt_u
                        if kpt.s == val_s])
    occ_cond = np.array([kpt.f_n[cond_n]/kpt.weight for kpt in calc.wfs.kpt_u
                         if kpt.s == cond_s])
    #if np.amin(occ_val) < epsilon and np.amax(occ_cond) > epsilon:
    if np.amax(occ_cond) > epsilon:
        return 0
    else:
        return lumo[0] - homo[0]

def get_direct_band_gap(calc, epsilon=1e-2):
    if isinstance(calc, (str, unicode)):
        calc = GPAW(calc, txt=None)

    gap = get_band_gap(calc, epsilon)
    if gap == 0:
        return 0

    nkpts = len(calc.get_ibz_k_points())
    nspins = calc.get_number_of_spins()
    
    homo, lumo = get_homo_lumo_states(calc)
    svbm = homo[1]
    nvbm = homo[3]
    scbm = lumo[1]
    ncbm = lumo[3]
    
    e_skn = np.array([[calc.get_eigenvalues(kpt=k, spin=svbm)
                      for k in range(nkpts)]
                      for s in range(nspins)])
    dir_gap = np.amin(e_skn[scbm, :, ncbm] - e_skn[svbm, :, ncbm])
    return dir_gap

def get_effective_mass(calc, kpt, spin, band, direction, npts=1):
    """
    Calculates the effective mass at a given k-point for the spin and band
    specified.

    calc: string or PAW object
        Calculator object or string representing the path to a .gpw file.
    kpt: int
        k-point index.
    spin: int
        Spin index.
    band: int
        Band index.
    direction: ndarray
        Direction given in reciprocal lattice vectors.
    """
    if isinstance(calc, basestring):
        calc = GPAW(calc, txt=None)

    dir_c = direction

    cell_cv = calc.atoms.get_cell()
    rcell = 2 * np.pi * np.linalg.inv(cell_cv).transpose()
    kd = calc.wfs.kd
    symm = kd.symmetry

    kpt_kc = kd.ibzk_kc[kpt]
    dir_v = np.dot(dir_c, rcell)

    kdst_k = np.zeros(2*npts + 1)
    eps_k = np.zeros(2*npts + 1)
    for k, d in enumerate(range(-npts, npts+1)):
        bzk_c = kpt_kc + d*dir_c

        # Find the corresponding k-point in 1st BZ
        bz1k_c = to1bz(np.array([bzk_c]), cell_cv)[0]
        # Use symmetry to find the related k-point in the IBZ
        ibzkpt, iop, timerev, diff_c = kd.find_ibzkpt(symm.op_scc,
                                                      kd.ibzk_kc,
                                                      bz1k_c)

        # Save rec. space distance and energy
        kdst_k[k] = d * np.sqrt(np.sum(dir_v**2))
        eps_k[k] = calc.get_eigenvalues(kpt=ibzkpt, spin=spin)[band]

    # Fit a parabola to the data points
    p = np.polyfit(kdst_k * Bohr, eps_k / Hartree, 2)

    # Return the effective mass
    return 0.5 / p[0]
        
