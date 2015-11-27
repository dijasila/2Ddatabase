import numpy as np

def get_kpt_path(calc, wpts_n):
    """Gets the bz k-points included in the calculator along a path between a given
    set of reciprocal space waypoints (typically high symmetry points). Useful for
    extracting a band  structure from a calculation with only Monkhorst-Pack
    k-point sampling.

    Keyword arguments:
    calc -- The PAW calculator object
    wpts_n -- A list of indices of the bz k-points to use as waypoints
    """

    bz_kpts = calc.get_bz_k_points()
    bz_rkpts = np.dot(bz_kpts, calc.atoms.get_reciprocal_cell())

    pts = []
    x = 0.

    def create_path_point(n, x):
        """Creates a k-point data object that holds the index, k-point and path
        distance"""
        bz_n = n
        ibz_n = calc.wfs.kd.bz2ibz_k[bz_n]
        kpt = bz_kpts[bz_n]
        rkpt = bz_rkpts[bz_n]
        return {'bz_n': bz_n,
                'ibz_n': ibz_n,
                'kpt': kpt,
                'rkpt': rkpt,
                'x': x}

    # Add the first waypoint to the path list
    pts.append(create_path_point(wpts_n[0], x))
    wpts_x = [x]

    # Loop through each of the subpaths
    for nwpt in range(1, len(wpts_n)):
        to_n = wpts_n[nwpt]
        frm_n = wpts_n[nwpt-1]
        to = bz_rkpts[to_n]
        frm = bz_rkpts[frm_n]
        vec = to - frm
        for n, k in enumerate(bz_rkpts):
            kvec = k - frm
            # Check if the point is on the line through to and frm
            if (np.linalg.norm(np.cross(k - frm, to - frm)) < 1e-13
                and not n in wpts_n):
                kdst = np.dot(vec, kvec) / np.linalg.norm(vec)**2
                # Check if the point is between to and frm
                if kdst > 0 and kdst < 1:
                    pts.append(create_path_point(n, x + np.linalg.norm(kvec)))

        # Add next waypoint to path
        x += np.linalg.norm(vec)
        pts.append(create_path_point(to_n, x))
        wpts_x.append(x)

    # Sort point according to their distance along the path
    pts.sort(key=lambda pnt: pnt['x'])
    return pts, wpts_x
