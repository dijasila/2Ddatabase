from __future__ import print_function
import os
from math import pi, sqrt, acos
import numpy as np
import ase.db
import ase.io
from ase.dft.bandgap import get_band_gap
import subprocess
from gpaw import GPAW, restart
#from materials import honeycombs, hexenes, hexanes, icsds, specials

xc = 'PBE'
ecut = 600
vacuum = 15
non_symmorphic = False

a0 = 3.2
N0relax = 18
N0gs = 24
l0 = 1. / a0
tol = 0.01

queue = 'verylong'
nodes = 2
ppn = 4
cpu = 'opteron4'

root = '/home/niflheim2/kiran/2ddatabase/'
relaxroot = '/home/niflheim2/mohpa/2D_Halides/MX2/Monolayer/'

script = root + 'gs.py'
names = list(np.loadtxt(root + 'structures.txt', dtype=str))
magnetic = list(np.loadtxt(root + 'magneticlist.txt', dtype=str))
metallic = list(np.loadtxt(root + 'metalliclist.txt', dtype=str))
done = ['ZnBr2']

names = [x for x in names if x not in magnetic and x not in metallic and x not in done]


#db = ase.db.connect(root + '2ddb.db')

#for row in db.select():
#    if not row.name in names:
#        names.append(row.name)
nkmin = 6
fdmetallic = open('metalliclist.txt', 'a')
fdmagnetic=open('magneticlist.txt', 'a')

for n, name in enumerate(magnetic):
    #rows = list(db.select(name=name))
    #if len(rows) > 0:
    #    row = rows[0]
    #else:
    #    continue
    
    #atoms = row.toatoms() 
    print(name + ':' )
    spinpol=True
    relaxdir = relaxroot + name
    relaxtraj = relaxdir + '/%s.traj' % name
    relaxgpw = relaxdir + '/%s.gpw' % name

    gsdir = root + 'data/magnetic/%s/' % (name)
    gsgpw = gsdir + 'PBE_gs.gpw'
    gstraj = gsdir + 'PBE_gs.traj'
    gsgpwspinpol = gsdir + 'PBE_gs_spinpol.gpw'
    gstrajspinpol = gsdir + 'PBE_gs_spinpol.traj'
    
    if os.path.isfile(gsgpw):
        print('    has gs')
        calc =  GPAW(relaxgpw, xc='PBE', txt=None)
        atoms = ase.io.read(gstraj)
        maxmagmom = max(np.abs(atoms.get_magnetic_moments()))
        print('    max magmom = %1.2f' % (maxmagmom))
        if maxmagmom < 0.01: # check if magnetic
            print('    %s is not magnetic!' % name)
            os.rename(gsdir, root + 'data/' + name)

        # Check band gap:
        bandgap = get_band_gap(calc, output=None)
        print('    bandgap = %1.2f' % bandgap[0])
        if bandgap[0] == 0: 
            print('    %s is metallic!' % name)
            print(name, file=fdmetallic)
            os.rename(gsdir, root + 'data/metallic/' + name)

        """
        os.rename(gsgpw, gsgpwspinpol)
        os.rename(gstraj, gstrajspinpol)
        if os.path.isfile(gsdir + 'PBE_gs0.gpw'):
            os.remove(gsdir + 'PBE_gs0.gpw')
            os.remove(gsdir + 'PBE_gs0.traj')
        """
        continue
    
    if os.path.isfile(relaxtraj):
        atoms = ase.io.read(relaxtraj)
    else:
        print('1st step not done!')
        continue
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
    nkx_relax = 18
    nky_relax = 18

    atomsfile=relaxgpw

    #if spinpol:
    #    maxmagmom = max(np.abs(atoms.get_magnetic_moments()))
    #    print('%s, max magmom=%s' % (name, maxmagmom))
    #    if max(np.abs(atoms.get_magnetic_moments())) < 0.01:
    #        spinpol = False

    if not os.path.isdir(gsdir):
        os.makedirs(gsdir)

    jobname = '%s_%s_gs' % (name, xc)
    job = script
    job += ' -k %d,%d' % (nkx_gs, nky_gs)
    job += ' -x %s' % xc
    job += ' -c %s' % ecut
    job += ' -v %s' % vacuum

    if spinpol:
        job += ' -s'
    if non_symmorphic:
        job += ' -n'

    cmd = 'gpaw-qsub -N %s -m ae' % jobname
    cmd += ' -q %s' % queue
    cmd += ' -l nodes=%d:ppn=%d:%s' % (nodes, ppn, cpu)
    cmd += ' %s %s %s' % (job, name, atomsfile)

    os.chdir(gsdir)
    sp = subprocess.call(cmd, shell=True)


#    sp = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
#    sp.communicate()
