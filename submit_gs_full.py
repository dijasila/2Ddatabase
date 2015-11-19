import os
import ase.db
import ase.io
import subprocess
from materials import honeycombs, hexenes, hexanes, specials, icsds

nkx0 = 16
nky0 = 24
nkx = 42
nky = 58
xc = 'PBE'
nbecut = 150
non_symmorphic = False
scalapack = True

queue = 'verylong'
nodes = 4
ppn = 8
cpu = 'xeon8'

root = '/home/niflheim2/fras/2D_database/'
script = root + 'gs_full.py'

#todo = honeycombs + hexenes + hexanes

todo = ['Phosphorene']

"""
names = ['Graphane', 'Silicene', 'Silicane', 'Germanene', 'Germanane',
         'Tinene', 'Tinane', 'SiC', 'GeC', 'SnGe', 'SnGe', 'SnSi', 'SnC',
         'BN', 'AlN', 'GaN', 'InN', 'InP', 'InAs', 'InSb', 'GaAs', 'BP',
         'BAs', 'GaP', 'AlSb', 'BSb']
"""
#db = ase.db.connect(root + '2ddb.db')
#for row in db.select(relaxed=True):
for name in todo:
    #print('%s, nks=%s, has_gs=%s' % (row.name, row.relax_nks, row.has_gs))
    #name = row.name
    #if name not in todo:
    #    continue
    #nks = [int(x) for x in row.gs_nks.split('x')]
    #nkx = nks[0]
    #nky = nks[1]
    #nkx = 18
    #nky = 18
    
    mydir0 = root + 'data/%s/%dx%dx1/' % (name, nkx0, nky0)
    gsfile = mydir0 + '%s_%s_gs.gpw' % (name, xc)
    if not os.path.isfile(gsfile):
        print('Does not have gs .gpw file for %s' % name)
        continue

    mydir = root + 'data/%s/%dx%dx1/' % (name, nkx, nky)
    if not os.path.isdir(mydir):
        os.makedirs(mydir)

    jobname = '%s_%s_gs_full' % (name, xc)
    job = script
    if nbecut:
        job += ' -b %s' % nbecut
    if nkx != nkx0 or nky != nky0:
        job += ' -k %d,%d' % (nkx, nky)
    if non_symmorphic:
        job += ' -s'
    if scalapack:
        job += ' -p'

    cmd = 'gpaw-qsub0 -N %s -m ae' % jobname
    cmd += ' -q %s' % queue
    cmd += ' -l nodes=%d:ppn=%d:%s' % (nodes, ppn, cpu)
    cmd += ' %s %s' % (job, gsfile)
    
    os.chdir(mydir)
    
    sp = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
    sp.communicate()
