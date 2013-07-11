import readhalos.RSDataReaderv2 as RSDataReader
import readhalos.readsubf as readsubf
import readsnapshots.readsnapHDF5 as rsHD
import readsnapshots.readids as readids
import readsnapshots.readsnap as rs
import numpy as np
from matplotlib import *
import pylab as plt
import sys
from grifflib import *
from mayavi import mlab



# DEFINE PARAMS
boxwidth = 3. # width of zoomed region
nhalo = 3098
snapnum = 63
titlestr = "halo" + str(nhalo)

# SET PATHS
#halopath = '/n/scratch2/hernquist_lab/bgriffen/caterpillar/parent/RockstarData'
#filename = '/n/home01/bgriffen/data/caterpillar/halos/candidates.dat'

halopath = '/spacebase/data/AnnaGroup/caterpillar/parent/RockstarData'
filename = '/spacebase/data/bgriffen/code/caterpillar/candidates.dat'

# LOAD HALOS
halodata = RSDataReader.RSDataReader(halopath,snapnum,digits=2)
allhalos = halodata.get_hosts()

xpos = np.array(allhalos['posX'])
ypos = np.array(allhalos['posY'])
zpos = np.array(allhalos['posZ'])
rvir = np.array(allhalos['rvir'])

# LOAD CANDIDATES
listin = []
for line in open(filename,'r'):
     li=line.strip()
     if not li.startswith("#"):
         line = line.partition('#')[0]
         listin.append(np.array(line.split(' ')[0:6]))

listin = np.array(listin)
listin = listin.astype(np.float)

# SPECIFY CANDIDATE CHOICE
idcand = listin[:,0]
candindex = idcand[nhalo]
rvircand = allhalos.ix[candindex]['rvir']

print "Halo ID:",candindex
print "Halo Mass:",listin[nhalo,1]
print "Halo virial radius:",rvircand
print "Halo z = 0 position:",listin[nhalo,3:6]
xmid = listin[nhalo,3]
ymid = listin[nhalo,4]
zmid = listin[nhalo,5]
rvircand = rvircand/1000

xposclose = []
yposclose = []
zposclose = []
rvirclose = []


R = np.sqrt((xmid-xpos)**2 + (ymid-ypos)**2 + (zmid-zpos)**2)
for i in xrange(0,len(xpos)):
    if R[i] < boxwidth:
        xposclose.append(xpos[i])
        yposclose.append(ypos[i])
        zposclose.append(zpos[i])
        rvirclose.append(rvir[i])

xposclose = np.array(xposclose)
yposclose = np.array(yposclose)
zposclose = np.array(zposclose)
rvirclose = np.array(rvirclose)
boxwidth = boxwidth/2
extent = [xmid-boxwidth,xmid+boxwidth,ymid-boxwidth,ymid+boxwidth,zmid-boxwidth,zmid+boxwidth]
#mlab.points3d(xmid, ymid, zmid, rvircand)
mlab.points3d(xposclose, yposclose, zposclose, rvirclose/1000,colormap="copper")
#mlab.xlabel('x-pos')
#mlab.ylabel('y-pos')
#mlab.zlabel('z-pos')
#mlab.colorbar(orientation='vertical',title='rvir')
mlab.outline(extent = extent)
mlab.axes(extent = extent) 
mlab.orientation_axes(xlabel='x-pos',ylabel='x-pos',zlabel='z-pos')
#mlab.orientation_axes()
mlab.show()

#plt.show()
