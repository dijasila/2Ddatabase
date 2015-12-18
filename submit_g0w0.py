import os
import subprocess
from math import pi
import numpy as np
import ase.db
from ase.units import Bohr
#from materials import *

rrad = 0
cgrid = 1
ecut = 100.
extrapolate = True#False

#isotropic = True
#isotropic = False
#suffix = '_iso'
suffix = None

queue = 'verylong'
nodes = 8
ppn = 8
cpu = 'xeon8'
nblocks = 8

#gs_nks = [42, 58]

#db = ase.db.connect('2ddb.db')

root = '/home/niflheim2/kiran/2ddatabase/'
script = root + 'g0w0.py'

#names = np.loadtxt(root + 'structures.txt', dtype=str)
todo =  list(np.loadtxt(root + 'doneGsFull.txt', dtype=str))


#for row in db.select(has_gs=True):
for name in todo[6:]:
    #name = row.name
    #if name not in todo:
    #    continue

    #gs_nks = [int(x) for x in row.gs_nks.split('x')]
    gsdir = root + 'data/%s/' % (name)
    gsfile0 = gsdir + 'PBE_gs.gpw'
    gsfile = gsdir + 'PBE_gs_full.gpw'
    
    if not os.path.isfile(gsfile):
        print('%s does not have full gs' % name)
        continue
    
    gwdir = gsdir + 'gw_ecut%03.f/' % ecut
    if not os.path.isdir(gwdir):
        os.makedirs(gwdir)

    gwfile = gwdir + 'PBE_g0w0%s.pckl' % (suffix)
    if os.path.isfile(gwfile):
        print('%s has G0W0' % name)
        continue

    #atoms = row.toatoms()
    #rcell_cv = 2 * pi * atoms.get_reciprocal_cell() * Bohr
    #b1 = np.linalg.norm(rcell_cv[0])
    #b2 = np.linalg.norm(rcell_cv[1])
    #frad = rrad * min(b1, b2)

    job = script
    job += ' -c %s' % ecut
    #job += ' -g %d' % cgrid
    #job += ' -r %s' % rrad
    job += ' -n %d' % nblocks
    if not extrapolate:
        job += ' -x'
    #if isotropic:
    #    job += ' -i'
    if suffix:
        job += ' -s %s' % suffix
    job += ' %s %s' % (name, gsfile)

    jobname = '%s_g0w0_ecut%03.f' % (name, ecut)

    cmd = 'gpaw-qsub'
    cmd += ' -N %s' % jobname
    cmd += ' -q %s' % queue
    cmd += ' -l nodes=%d:ppn=%d:%s' % (nodes, ppn, cpu)
    cmd += ' -m ae'
    cmd += ' ' + job

    os.chdir(gwdir)
    sp = subprocess.call(cmd, shell=True)

    #sp = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
    #sp.communicate()
