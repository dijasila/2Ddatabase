import numpy as np
import os
import pickle
from math import pi
from optparse import OptionParser
from ase.parallel import paropen
from ase.dft.kpoints import monkhorst_pack
from ase.units import Hartree, Bohr
from gpaw import GPAW
from gpaw.kpt_descriptor import KPointDescriptor, to1bz
from gpaw.mpi import serial_comm
from gpaw.response.qp import GWQuasiParticleCalculator
from gpaw.response.selfenergy import QuadQPointIntegration
from band_structure import get_homo_lumo_states
from qpt_grids import two_density_grid_2d, voronoi_cell_weights_2d

parser = OptionParser()
parser.add_option('-c', '--ecut', type='float', dest='ecut', default=50.)
parser.add_option('-x', '--noxp', action='store_true',
                  dest='noxp', default=False)
parser.add_option('-g', '--cgrid', type='int', dest='cgrid', default=1)
parser.add_option('-r', '--rad', type='float', dest='rad', default=0.0)
parser.add_option('-n', '--nblocks', type='int', dest='nblocks', default=1)
parser.add_option('-i', '--iso', action='store_true', dest='iso', default=False)
parser.add_option('-s', '--suffix', dest='suffix', default='')

opts, args = parser.parse_args()

name = args[0]
gsfile = args[1]

if not opts.noxp:
    # Extrapolate
    ecut = opts.ecut
else:
    # No extrapolation
    ecut = [opts.ecut]

calc = GPAW(gsfile, communicator=serial_comm, txt=None)
atoms = calc.get_atoms()
xc = calc.get_xc_functional()
ibzk_kc = calc.get_ibz_k_points()

filename = '%s_%s_g0w0' % (name, xc)
filename += opts.suffix

homo, lumo = get_homo_lumo_states(calc)

nvb = homo[3]
ncb = lumo[3]

na = max(0, nvb - 4)
nb = ncb + 5
bandrange = (na, nb)


kd = calc.wfs.kd

offset_c = 0.5 * ((kd.N_c + 1) % 2) / kd.N_c
bzq_qc = monkhorst_pack(kd.N_c) + offset_c
qd = KPointDescriptor(bzq_qc)
qd.set_symmetry(calc.atoms, kd.symmetry)
cell_cv = calc.wfs.gd.cell_cv
rcell_cv = 2 * pi * np.linalg.inv(cell_cv).T

bzq_qv = np.dot(bzq_qc, rcell_cv)

#frad = opts.rad
b1 = np.linalg.norm(rcell_cv[0])
b2 = np.linalg.norm(rcell_cv[1])
frad = opts.rad * min(b1, b2)
cgrid = opts.cgrid

grid_qpts = two_density_grid_2d(qd, rcell_cv, cgrid, frad)

gridq_qc = bzq_qc[grid_qpts]
gridq_qv = bzq_qv[grid_qpts]

Nqpts = np.prod(qd.N_c)
fweight = 1. / np.prod(qd.N_c[0:2])
cweight = 1. / np.prod(1.0 * qd.N_c[0:2] / cgrid)
weights_q = voronoi_cell_weights_2d(gridq_qv, rcell_cv, cweight)

anisotropic = not opts.iso
qptint = QuadQPointIntegration(qd, cell_cv, gridq_qc, weight_q=weights_q,
                               anisotropic=anisotropic)

qpcalc = GWQuasiParticleCalculator(calc=calc,
                                   filename=filename,
                                   txt=filename + '.txt',
                                   temp=True,
                                   savechi0=False,
                                   savepair=False,
                                   bandrange=bandrange,
                                   ecut=ecut,
                                   domega0=0.1,
                                   omega2=15.,
                                   truncation='2D',
                                   qptint=qptint,
                                   nblocks=opts.nblocks)

qpcalc.load_iteration()
qp_skn = qpcalc.get_qp_energies()
qpcalc.save_restart_file()
data = {#'grid_qpts': grid_qpts,
        'bandrange': bandrange,
        'ibzk_kc': ibzk_kc,
        'rcell_cv': rcell_cv * Bohr,
        'eps_skn': qpcalc.eps_skn * Hartree,
        'qp_skn': qp_skn,
        'vxc_skn': qpcalc.vxc_skn * Hartree,
        'exx_skn': qpcalc.exx_skn * Hartree,
        'sigma_skn': qpcalc.sigma_skn * Hartree,
        'sigerr_skn': qpcalc.selfenergy.sigerr_skn * Hartree,
        'sigma_iskn': qpcalc.selfenergy.sigma_iskn * Hartree,
        'Z_skn': qpcalc.Z_skn}
pickle.dump(data, paropen(filename + '_result.pckl', 'wb'))
