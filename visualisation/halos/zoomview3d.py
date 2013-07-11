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
haloid = 80609
snapnum = 63
titlestr = "halo" + str(nhalo)

ictype = 'box'
res = 11
nvir = 9
pad = 7
hubble = 0.6711
basedir = '/spacebase/data/bgriffen/data/caterpillar/halos/' + titlestr + '/' + ictype + '/l' + str(res) + '/p' + str(pad) + '/nvir' + str(nvir)

tmppath = basedir + '/' + ictype + '/l' + str(res) + '/p' + str(pad) + '/nvir' + str(nvir) + '/outputs'

gadget = basedir + '/outputs'
halopath = basedir + 'RockstarData'
filename = '/spacebase/data/bgriffen/code/caterpillar/candidates.dat'

# LOAD HALOS
s = readsubf.subfind_catalog(tmppath, snapnum)
id = readsubf.subf_ids(tmppath, snapnum, 0, 0, read_all=1)
rvir = s.sub_halfmassrad
xpos = s.sub_pos[:,0]
ypos = s.sub_pos[:,1]
zpos = s.sub_pos[:,2]
m200 = s.sub_mass*10**10/hubble

mask = np.argsort(m200)
m200big = m200[mask[-topNhalos:]]

for inx in xrange(0,len(m200big)):
                if m200big[inx] > 1.15E12 and m200big[inx] < 1.25E12:
                    if xposbig[inx] > xcen-deltahalo and xposbig[inx] < xcen+deltahalo:
                        if yposbig[inx] > ycen-deltahalo and yposbig[inx] < ycen+deltahalo:
                            if zposbig[inx] > zcen-deltahalo and zposbig[inx] < zcen+deltahalo:
                                xmid = xposbig[inx]
                                ymid = yposbig[inx]
                                zmid = zposbig[inx]
                                rvircand = rvirbig[inx]
                                #mvirsubhost = m200big[inx]

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

plt.show()
