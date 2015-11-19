import os
from math import pi, sqrt, acos
import numpy as np
import ase.db
import ase.io
import subprocess
from gpaw import GPAW
from materials import honeycombs, hexenes, hexanes, icsds, specials

#nkx = 36
#nky = 36
xc = 'PBE'
ecut = 600.
vacuum = 15
non_symmorphic = False

#a0 = 2.5
a0 = 3.2
N0relax = 18
N0gs = 24
l0 = 1. / a0

queue = 'medium'
nodes = 2
ppn = 4
cpu = 'opteron4'

root = '/home/niflheim2/fras/2D_database/'
script = root + 'gs.py'

#names = honeycombs + hexenes + hexanes + icsds
names = honeycombs + hexenes + hexanes
#rmlist = ['SrIF', 'Sb2Ge2Te5', 'FePSe3', 'FeTe', 'PFeLi', 'FeBr3', 'P2CuSe6Bi',
#          'CrSiTe3', 'Cu2S', 'P2AgSe6Bi', 'FeS', 'KC6FeO3N3', 'Ag2ReCl6',
#          'LiFeAs', 'Ni2Te2Sb', 'ZnIn2S4', 'FeSe', 'Bi14Te13S8', 'ScP2AgSe6']

#for rmname in rmlist:
#    if rmname in names:
#        names.remove(rmname)
#names = ['Graphane', 'Silicene', 'Silicane', 'Germanene', 'Germanane',
#         'Tinene', 'Tinane', 'SiC', 'GeC', 'SnGe', 'SnGe', 'SnSi', 'SnC',
#         'BN', 'AlN', 'GaN', 'InN', 'InP', 'InAs', 'InSb', 'GaAs', 'BP',
#         'BAs', 'GaP', 'AlSb', 'BSb']

db = ase.db.connect(root + '2ddb.db')
#for row in db.select():
#    if not row.name in names:
#        names.append(row.name)
nkmin = 6

names = ['Phosphorene']
for name in names[:]:
    rows = list(db.select(name=name))
    if len(rows) > 0:
        row = rows[0]
    else:
        continue
    
    atoms = row.toatoms()
    cell = atoms.get_cell()
    icell = np.linalg.inv(cell).T
    b1 = icell[0]
    b2 = icell[1]
    l1 = np.linalg.norm(b1)
    l2 = np.linalg.norm(b2)
    print('l1=%s, l2=%s' % (l1, l2))
    ang = 180 / pi * acos(np.dot(b1, b2) / (l1 * l2))
    nkx_relax = 18
    nky_relax = 18
    if (abs(ang - 120) < 0.5 or abs(ang - 60) < 0.5) and abs(l1 - l2) < tol:
        print('Hexagonal')
        # Hexagonal unit cell
        nkx_gs = max(int(round(N0gs * l1 / l0 / 6) * 6), nkmin)
        nky_gs = max(int(round(N0gs * l2 / l0 / 6) * 6), nkmin)
    else:
        nkx_gs = max(int(round(N0gs * l1 / l0 / 2) * 2), nkmin)
        nky_gs = max(int(round(N0gs * l2 / l0 / 2) * 2), nkmin)
    print('%s: %dx%d -> %dx%d' % (name, nkx_relax, nky_relax, nkx_gs, nky_gs))
    #print(name)
    nkx_relax = 18
    nky_relax = 18
    #nkx_gs = nkmin
    #nky_gs = nkmin
    
    gsdir = root + 'data/%s/%dx%dx1/' % (name, nkx_gs, nky_gs)
    gsgpw = gsdir + '%s_PBE_gs.gpw' % name
    if os.path.isfile(gsgpw):
        print('has gs')
        #continue
    
    relaxdir = root + 'data/%s/%dx%dx1/' % (name, nkx_relax, nky_relax)
    relaxtraj = relaxdir + '%s_PBE_relax.traj' % name
    relaxgpw = relaxdir + '%s_PBE_relax.gpw' % name
    
    spinpol = True

    atoms = None
    atomsfile = None
    try:
        calc = GPAW(relaxgpw, txt=None)
        atoms = calc.get_atoms()
        if calc.get_number_of_spins() == 1:
            spinpol = False
        atomsfile = relaxgpw
    except Exception as e:
        #print(e)
        pass

    try:
        atoms = ase.io.read(relaxtraj)
        atomsfile = relaxtraj
    except Exception as e:
        #print(e)
        pass

    """
    if atoms is None or atomsfile is None:
        relaxdir = root + 'data/%s/18x18x1/' % name
        relaxtraj = relaxdir + '%s_PBE_relax.traj' % name
        relaxgpw = relaxdir + '%s_PBE_relax.gpw' % name
        try:
            calc = GPAW(relaxgpw, txt=None)
            atoms = calc.get_atoms()
            #atoms = ase.io.read(relaxgpw)
            atomsfile = relaxgpw
        except Exception as e:
            #print(e)
            pass

        try:
            atoms = ase.io.read(relaxtraj)
            atomsfile = relaxtraj
        except Exception as e:
            #print(e)
            pass
    """

    if atoms is None or atomsfile is None:
        print('%s has no atomsfile' % name)
        continue

    if spinpol:
        maxmagmom = max(np.abs(atoms.get_magnetic_moments()))
        print('%s, max magmom=%s' % (name, maxmagmom))
        if max(np.abs(atoms.get_magnetic_moments())) < 0.01:
            spinpol = False

    if not os.path.isdir(gsdir):
        os.makedirs(gsdir)

    jobname = '%s_%s_gs' % (name, xc)
    job = script
    job += ' -k %d,%d' % (nkx_gs, nky_gs)
    job += ' -x %s' % xc
    job += ' -c %s' % ecut
    #job += ' -v %s' % vacuum
    if spinpol:
        job += ' -s'
    if non_symmorphic:
        job += ' -n'

    cmd = 'gpaw-qsub0 -N %s -m ae' % jobname
    cmd += ' -q %s' % queue
    cmd += ' -l nodes=%d:ppn=%d:%s' % (nodes, ppn, cpu)
    cmd += ' %s %s %s' % (job, name, atomsfile)
    
    os.chdir(gsdir)
    sp = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
    sp.communicate()
