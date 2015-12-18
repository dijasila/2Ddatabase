from __future__ import print_function
import os
import ase.db
import ase.io
import subprocess
import numpy as np
#from materials import honeycombs, hexenes, hexanes, specials, icsds

nkx0 = 18
nky0 = 18
nkx = 18
nky = 18

xc = 'PBE'
nbecut = 150
non_symmorphic = False
scalapack = True

queue = 'long'
nodes = 1
ppn = 16
cpu = 'xeon16' #'xeon8'

root = '/home/niflheim2/kiran/2ddatabase/'

script = root + 'gs_full.py'
names = np.loadtxt(root + 'structures.txt', dtype=str)
magnetic = list(np.loadtxt(root + 'magneticlist.txt', dtype=str))
metallic = list(np.loadtxt(root + 'metalliclist.txt', dtype=str))
done = list(np.loadtxt(root + 'doneGsFull.txt', dtype=str))
names = [x for x in names if x not in magnetic and x not in metallic
         and x not in done]

#db = ase.db.connect(root + '2ddb.db')
#for row in db.select(relaxed=True):

for name in names[:]:
    #print('%s, nks=%s, has_gs=%s' % (row.name, row.relax_nks, row.has_gs))
    #name = row.name
    #if name not in todo:
    #    continue
    #nks = [int(x) for x in row.gs_nks.split('x')]
    #nkx = nks[0]
    #nky = nks[1]
    #nkx = 18
    #nky = 18
    print(name)
    mydir = root + 'data/%s/' % (name)
    gsfile = mydir + '%s_gs.gpw' % (xc)
    gsfullfile = mydir + '%s_gs_full.gpw' % (xc)

    if not os.path.isfile(gsfile):
        print('Does not have gs .gpw file for %s' % name)
        continue
    
    if os.path.isfile(gsfullfile):
        print('full gs .gpw file already done!')
        continue

    #mydir = mydir0 #root + 'data/%s/' % (name)
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

    cmd = 'gpaw-qsub -N %s -m ae' % jobname
    cmd += ' -q %s' % queue
    cmd += ' -l nodes=%d:ppn=%d:%s' % (nodes, ppn, cpu)
    cmd += ' %s %s' % (job, gsfile)
    
    os.chdir(mydir)
    sp = subprocess.call(cmd, shell=True)
#    sp = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
#    sp.communicate()
