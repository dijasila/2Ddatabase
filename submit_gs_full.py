import os
import ase.db
import ase.io
import subprocess
from materials import honeycombs, hexenes, hexanes, specials, icsds

nkx0 = 18
nky0 = 18
nkx = 18
nky = 18

xc = 'PBE'
nbecut = 150
non_symmorphic = False
scalapack = True

queue = 'verylong'
nodes = 4
ppn = 8
cpu = 'xeon8'

root = '/home/niflheim2/kiran/2ddatabase/'

script = root + 'gs_full.py'
names = np.loadtxt(root + 'structures.txt', dtype=str)

todo = names[:1] 

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
    
    mydir0 = root + '%s/' % (name)
    gsfile = mydir0 + '%s_gs.gpw' % (xc)
    if not os.path.isfile(gsfile):
        print('Does not have gs .gpw file for %s' % name)
        continue

    mydir = root + '%s/' % (name)
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
    
    sp = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
    sp.communicate()
