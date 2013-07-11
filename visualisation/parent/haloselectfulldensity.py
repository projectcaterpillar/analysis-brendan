import brendanlib.gadgetlibomp as glib
from brendanlib.grifflib import makecolormap,drawcircle,getcandidatelist
import readsnapshots.readidsHDF5
import readsnapshots.readsnapHDF5 as rs
import numpy as np
import pylab as plt
import readhalos.RSDataReaderv2 as RSDataReader
import sys
import matplotlib 

dim = 1000
Noutputs = 64
delta = 3.0
boxsize = 100.0

datadir = '/spacebase/data/AnnaGroup/caterpillar/parent/512Parent/outputs/'
halopath = '/spacebase/data/AnnaGroup/caterpillar/parent/RockstarData' # Spacebase
#candfile = '/n/home01/bgriffen/data/caterpillar/halos/candidates.dat' #odyssey

candfile = '/spacebase/data/bgriffen/projects/caterpillar/analysis/selection/candidates.dat' #spacebase
candlist = getcandidatelist(candfile)
candid = candlist[:,0]

halodata = RSDataReader.RSDataReader(halopath,63,digits=2)
allhalos = halodata.get_hosts()
candhalos = allhalos.ix[candid]
xhalos = np.array(candhalos['posX'])
yhalos = np.array(candhalos['posY'])

sampleid = np.array([190897,208737,140666,28221,28188,147273,78411,131988,19910,147419])
subsamphalos = allhalos.ix[sampleid]
xsamp = np.array(subsamphalos['posX'])
ysamp = np.array(subsamphalos['posY'])

fig = plt.figure(figsize=(15,6))
fig.subplots_adjust(hspace=0.15,wspace=0.08)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.plot(xhalos,yhalos,'bo',markeredgewidth=0.0,markersize=1)
ax1.plot(xsamp,ysamp,'r*',markeredgewidth=0.0,markersize=12)
ax1.set_xlim([0,100])
ax1.set_ylim([0,100])
ax1.set_xlabel('x-pos [Mpc/h]',fontsize=14)
ax1.set_ylabel('y-pos [Mpc/h]',fontsize=14)

interplist = ['bilinear','nearest','bicubic','spline16','hanning','hamming','hermite','kaiser','quadric','catrom','gaussian','bessel','mitchell','sinc','lanczos']

for i in range(63,Noutputs):
    snapshot = '%.3i' % i
    
    basedir = 'snapdir_' + str(snapshot)
    filename =  '/snap_' + str(snapshot)

    totname = datadir + basedir + filename
    print totname
    print "Reading from directory...", totname
    print "#####################################"
    print "            LOADING HEADER           "
    print "#####################################"
    header = rs.snapshot_header(totname)
    print "  Nparttot:",'{0:.2e}'.format(float(sum(header.nall))) 
    print "     Mpart:",'{0:.2e}'.format(sum(header.massarr)*10**10)
    print " Exp. Fact:",header.time
    print "   Boxsize:",header.boxsize
    print "  Redshift:",'{:.2f}'.format(float(header.redshift))
    figname = './images/dim' + str(dim) + '_snap' + str(snapshot) + '_z' + '{:.2f}'.format(header.redshift)
    print "PDF to be created:",figname
    #print "#####################################"
    #print "         CALCULATING MESH            "
    print "#####################################"
    mpart = sum(header.massarr)             ###*10**10
    data = np.zeros((dim,dim,dim))
    Nfiles = header.filenum
    Nfiles = header.filenum
    tmp = np.zeros((dim,dim,dim))
    blockfilename = totname #+ '.' + str(i)  + '.hdf5'
    pos = rs.read_block(blockfilename, "POS ")
    #npart = len(pos)
    #tmp = glib.cic(pos,mpart,boxsize,dim)
    #data = data + tmp
    
    #imageArray2Dmesh = np.mean(data**2,axis=2)
    #max=10**-1.5*imageArray2Dmesh.max()
    #min=10**-6*max
    #print "min cic array:",imageArray2Dmesh.min()
    #print "max cic array:",imageArray2Dmesh.max()
    #vmin=np.log10(min)
    #vmax=np.log10(max)
    #imageArray2Dmesh[imageArray2Dmesh<min] = min
    #imageArray2Dmesh[imageArray2Dmesh>max] = max
    #imageArray2Dmesh=np.log10(imageArray2Dmesh)
    newmap = makecolormap()
    #imsave2(fname=figname+'CIC', arr=imageArray2Dmesh, origin='lower', vmin=vmin, vmax=vmax, cmap=newmap)

    heatmap, xedges, yedges = np.histogram2d(pos[:,0], pos[:,1], bins=dim)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    meshnotzero=heatmap[heatmap != 0]
    minval = meshnotzero.min()
    heatmap[heatmap<minval] = minval

    maxval=10**0*heatmap.max()
    minval=10**-3*maxval
    print "min heatmap:",heatmap.min()
    print "max heatmap:",heatmap.max()
    vmin=np.log10(minval)
    vmax=np.log10(maxval)

    heatmap = np.rot90(np.log10(heatmap)) #rotated again cbf to work out
    ax2.imshow(heatmap, extent=extent,cmap = newmap, vmin = vmin , vmax = vmax, origin='lower',interpolation='catrom')

    #for interp in interplist:
    #    fig2 = plt.figure(frameon=False)
    #    ax3 = fig2.add_subplot(111)
    #    ax3.imshow(heatmap, extent=extent,cmap = newmap, vmin = vmin , vmax = vmax, origin='lower',interpolation=interp)
    #    ax3.set_xticks([])
    #    ax3.set_yticks([])
    #    fig2.savefig('./images/interpolation/' + interp + '_dim' + str(dim) + '_snap' + str(snapshot) + '.pdf',bbox_inches='tight', format='pdf')
    #    fig2.clf()

ax2.set_xlim([0,100])
ax2.set_ylim([0,100])
ax2.set_xlabel('x-pos [Mpc/h]',fontsize=14)
ax2.set_ylabel('y-pos [Mpc/h]',fontsize=14)

fig.savefig(figname + '.pdf',bbox_inches='tight', format='pdf')

#imsave2(fname=figname+'HIST', arr=np.log10(heatmap), origin='lower', vmin=vmin, vmax=vmax, cmap=newmap)

    #print "   Boxsize:",header.boxsize
    #print "    NFiles:",header.filenum
    #print "    Omega0:",header.omega0
    #print "    OmegaL:",header.omegaL
    #print "    Hubble:",header.hubble
    #print "  SFR Flag:",header.sfr
    #print "Cool. Flag:",header.cooling
    #print "Stel. Flag:",header.stellar_age
    #print "Metal Flag:",header.metals
    #print "Feed. Flag:",header.feedback