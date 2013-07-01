#IMPORT CORE MODULES
#import gadgetlibomp
import numpy as np
import matplotlib.pyplot as plt
#from pylab import *
import sys, os, platform
import matplotlib
#IMPORT CATERPILLAR MODULES
import readhalos.readsubf as readsubf
import readhalos.readgroup as readgroup
import readhalos.RSDataReaderv2 as RSDataReader
import readsnapshots.readsnap as rs
import readsnapshots.readsnapHDF5 as rhdf5
import mergertrees.MTCatalogue as MT
#IMPORT PERSONAL LIBRARIS
from brendanlib.grifflib import *

basepath = determinebasepath(platform.node())
#haloid = 190897
haloid = 103794
#208737
hubble = 0.6711
startsnap = 50
endsnap = 63
snaplist = range(startsnap,endsnap+1)
boxwidth = 4./hubble

# SET PATHS
gadpath = basepath + 'caterpillar/parent/512Parent/outputs'
halopath = basepath + 'caterpillar/parent/RockstarData'
basedir = basepath + 'caterpillar/contamination/halos/halo' + str(haloid)
treefile = halopath + '/trees/tree.bin'
indexfile = halopath + '/trees/treeindex.csv'

cat = MT.MTCatalogue(halopath + '/trees',numHosts=6,indexbyrsid=True)
tree = cat[haloid]
mainbranch = tree.getMainBranch(0)
scalemb = mainbranch['scale']
rvirmb = mainbranch['rvir']
posXmb = mainbranch['posX']
posYmb = mainbranch['posY']
posZmb = mainbranch['posZ']
mvirmb = mainbranch['mvir']
mallmb = mainbranch['m200c_all']
m200bmb = mainbranch['m200b']
vmaxmb = mainbranch['vmax']
vrmsmb = mainbranch['vrms']
rsmb = mainbranch['rs']
xoffmb = mainbranch['xoff']
voffmb = mainbranch['voff']
virratmb = mainbranch['T/|U|']
btamb = mainbranch['b_to_a']
ctamb = mainbranch['c_to_a']
spinmb = mainbranch['spin']
spinbmb = mainbranch['spin_bullock']
pecvxmb = mainbranch['pecVX']
pecvymb = mainbranch['pecVY']
pecvzmb = mainbranch['pecVZ']
jxmb = mainbranch['Jx']
jymb = mainbranch['Jy']
jzmb = mainbranch['Jz']
Axmb = mainbranch['A[x]']
Aymb = mainbranch['A[y]']
Azmb = mainbranch['A[z]']
haloidmb = mainbranch['origid']
#fullbranch = halotree.Trees[80609].data
# LOAD PARENT HALOS
#halodata = RSDataReader.RSDataReader(halopath,61,digits=2)
#allhalos = halodata.get_hosts()
#haloid = 254188
#print "---------------------------------------"
#print "       rockstar id: ",haloid
#print "---------------------------------------"
#print "        x-pos:",'{:.2f}'.format(float(allhalos.ix[haloid]['posX'])), "   \ [Mpc/h]"
#print "        y-pos:",'{:.2f}'.format(float(allhalos.ix[haloid]['posY'])), "   \ [Mpc/h]"
#print "        z-pos:",'{:.2f}'.format(float(allhalos.ix[haloid]['posZ'])), "   \ [Mpc/h]"
#print "  virial mass:",'{0:.2e}'.format(float(allhalos.ix[haloid]['mvir'])),"\ [Msol/h]"
#print "virial radius:",'{:.2f}'.format(float(allhalos.ix[haloid]['rvir'])),"  \ [kpc]"
#print "---------------------------------------"
#del allhalos
#sys.exit()
# FIGURE DIMENSIONS
labelsize = 14
ticksize = 12
legendfontsize = 12

figdir = './figs/posevolution/' + '/halo' + str(haloid) + '/L' + str(int(boxwidth)) + '/'

if not os.path.exists(figdir):
    os.makedirs(figdir)

mvirvec = []
rvirvec = []
vmaxvec = []
scalevec = []
vrmsvec = []
rsvec = []
xoffvec = []
voffvec = []
virratvec = []
btavec = []
ctavec = []
spinvec = []
spinbvec = []
pecvxvec = []
pecvyvec = []
pecvzvec = []
jxvec = []
jyvec = []
jzvec = []
Axvec = []
Ayvec = []
Azvec = []
npartvec = []
Tvec = []
ctbvec = []
mallvec = []
m200bvec = []

for snap in snaplist:
    index = -99
    xcen = posXmb[-1]
    ycen = posYmb[-1]
    zcen = posZmb[-1]
    
    if snap < 10:
    	digits = 1
    else:
    	digits = 2

    fig = plt.figure(figsize=(33.,20.))

    ax1 = plt.subplot2grid((4, 6), (0, 0), colspan=2, rowspan=2)
    ax2 = plt.subplot2grid((4, 6), (0, 2), colspan=2, rowspan=2)
    ax3 = plt.subplot2grid((4, 6), (0, 4), colspan=2, rowspan=2)

    ax4 = plt.subplot2grid((4, 6), (2, 0))
    ax5 = plt.subplot2grid((4, 6), (2, 1))
    ax6 = plt.subplot2grid((4, 6), (2, 2))
    ax7 = plt.subplot2grid((4, 6), (2, 3))
    ax8 = plt.subplot2grid((4, 6), (2, 4))
    ax9 = plt.subplot2grid((4, 6), (2, 5))

    ax10 = plt.subplot2grid((4, 6), (3, 0))
    ax11 = plt.subplot2grid((4, 6), (3, 1))
    ax12 = plt.subplot2grid((4, 6), (3, 2))
    ax13 = plt.subplot2grid((4, 6), (3, 3))
    ax14 = plt.subplot2grid((4, 6), (3, 4))
    ax15 = plt.subplot2grid((4, 6), (3, 5))

    font = {'size': labelsize}
    matplotlib.rc('font', **font)

    header = rs.snapshot_header(gadpath + '/snapdir_0' + str(snap) + '/snap_0' + str(snap) + '.0')
    
    halos = RSDataReader.RSDataReader(halopath,snap,digits=digits)
    hosts = halos.get_hosts()

    redshift = '{:.2f}'.format(header.redshift)
    expfact = '{:.2f}'.format(header.time)
   
    titlestr = 'snapnum = ' + str(snap) + ' | redshift = ' + redshift + ' | a = ' + expfact
    print titlestr
    
    # CENTRAL HOST
    if len(np.array(hosts['posX'])) != 0:
        for i in range(0,len(mainbranch['scale'])):
        	if '{:.2f}'.format(mainbranch['scale'][i]) == expfact:
        		index = i

        if index != -99:
            xcen = posXmb[index]
            ycen = posYmb[index]
            zcen = posZmb[index]
            haloidsnap = haloidmb[index]

    # DO HOST HALOS
    if len(np.array(hosts['posX'])) != 0: 

        xposhosts = hosts['posX']
        yposhosts = hosts['posY']
        zposhosts = hosts['posZ']
        rvirhosts = hosts['rvir']
        mvirhosts = hosts['mvir']
        cond1h = (xposhosts >= xcen - boxwidth/2.) & (xposhosts <= xcen + boxwidth/2.)
        cond2h = (yposhosts >= ycen - boxwidth/2.) & (yposhosts <= ycen + boxwidth/2.)
        cond3h = (zposhosts >= zcen - boxwidth/2.) & (zposhosts <= zcen + boxwidth/2.)
        #cond4h = (hosts['id'] != haloidsnap)
        condh = cond1h & cond2h & cond3h 
        #& cond4h
        xposhosts = xposhosts[condh]
        yposhosts = yposhosts[condh]
        zposhosts = zposhosts[condh]
        rvirhosts = rvirhosts[condh]
        mvirhosts = mvirhosts[condh]
        if len(xposhosts) != 0:
            xcirchosts,ycirchosts = drawcircle(xposhosts,yposhosts,rvirhosts/1000,vec=True)
            ax1.plot(xcirchosts,ycirchosts,'b-',linewidth=1)
            xcirchosts,zcirchosts = drawcircle(xposhosts,zposhosts,rvirhosts/1000,vec=True)
            ax2.plot(xcirchosts,zcirchosts,'b-',linewidth=1)
            ycirchosts,zcirchosts = drawcircle(yposhosts,zposhosts,rvirhosts/1000,vec=True)
            ax3.plot(ycirchosts,zcirchosts,'b-',linewidth=1)
            print '-- # hosts:',len(xposhosts)

        for hostindex in range(0,len(xposhosts)):
            xposplace = np.array(xposhosts[condh])
            yposplace = np.array(yposhosts[condh])
            zposplace = np.array(zposhosts[condh])
            mvirplace = np.array(mvirhosts[condh])
            rvirplace = np.array(rvirhosts[condh])
            if mvirplace[hostindex] > 10**11:
                mstring = '{0:.2e}'.format(mvirplace[hostindex])
                xtrans = (xposplace[hostindex] - (xcen - boxwidth/2.))/boxwidth 
                ytrans = (yposplace[hostindex] - (ycen - boxwidth/2.))/boxwidth
                ztrans = (zposplace[hostindex] - (zcen - boxwidth/2.))/boxwidth
                extra = 1.5*rvirplace[hostindex]/1000
                addon = extra/boxwidth
                placetext(ax1,xtrans+addon,ytrans,mstring,'normal',11)
                placetext(ax2,xtrans+addon,ztrans,mstring,'normal',11)
                placetext(ax3,ytrans+addon,ztrans,mstring,'normal',11)
    
    #placetext(ax1,10,78,'hello')

    # DO SUBHALOS
    subs = halos.get_subs()
    if len(np.array(subs['posX'])) != 0: 
        xpossubs = subs['posX']
        ypossubs = subs['posY']
        zpossubs = subs['posZ']
        rvirsubs = subs['rvir']
        cond1s = (xpossubs >= xcen - boxwidth/2.) & (xpossubs <= xcen + boxwidth/2.)
        cond2s = (ypossubs >= ycen - boxwidth/2.) & (ypossubs <= ycen + boxwidth/2.)
        cond3s = (zpossubs >= zcen - boxwidth/2.) & (zpossubs <= zcen + boxwidth/2.)
        conds = cond1s & cond2s & cond3s
        xpossubs = xpossubs[conds]
        ypossubs = ypossubs[conds]
        zpossubs = zpossubs[conds]
        rvirsubs = rvirsubs[conds]
        if len(xpossubs) != 0:
            xcircsubs,ycircsubs = drawcircle(xpossubs,ypossubs,rvirsubs/1000,vec=True)
            ax1.plot(xcircsubs,ycircsubs,'g-',linewidth=1)
            xcircsubs,zcircsubs = drawcircle(xpossubs,zpossubs,rvirsubs/1000,vec=True)
            ax2.plot(xcircsubs,zcircsubs,'g-',linewidth=1)
            ycircsubs,zcircsubs = drawcircle(ypossubs,zpossubs,rvirsubs/1000,vec=True)
            ax3.plot(ycircsubs,zcircsubs,'g-',linewidth=1)
            print '-- # subs:  ',len(xpossubs)

    if index != -99:
        npartch = mvirmb[index]/6.554894E6
        rvirch = rvirmb[index]
        mvirch = mvirmb[index]
        mallch = mallmb[index]
        m200bch = m200bmb[index]
        posXch = posXmb[index]
        posYch = posYmb[index]
        posZch = posZmb[index]
        vmaxch = vmaxmb[index]
        vrmsch = vrmsmb[index]
        rsch = rsmb[index]
        xoffch = xoffmb[index]
        voffch = voffmb[index]
        virratch = virratmb[index]
        btach = btamb[index]
        ctach = ctamb[index]
        ctbch = btach*(1/ctach)
        spinch = spinmb[index]
        spinbch = spinbmb[index]
        pecvxch = pecvxmb[index]
        pecvych = pecvymb[index]
        pecvzch = pecvzmb[index]
        jxch = jxmb[index]
        jych = jymb[index]
        jzch = jzmb[index]
        Axch = Axmb[index]
        Aych = Aymb[index]
        Azch = Azmb[index]
        Tch = calcT(ctbch)


        mstring = '{0:.2e}'.format(mvirch)

        circut1 = 1.4*header.time*hubble
        curcut2 = 3.*header.time*hubble
        circut3 = 4.*header.time*hubble

        xch,ych = drawcircle(posXch,posYch,rvirch/1000,vec=False)
        ax1.plot(xch,ych,'r-',linewidth=3)
        xch,zch = drawcircle(posXch,posZch,rvirch/1000,vec=False)
        ax2.plot(xch,zch,'r-',linewidth=3)
        ych,zch = drawcircle(posYch,posZch,rvirch/1000,vec=False)
        ax3.plot(ych,zch,'r-',linewidth=3)

        xdot,ydot = drawcircle(posXch,posYch,circut1,vec=False)
        ax1.plot(xdot,ydot,'c--',linewidth=2,alpha=0.5)
        xdot,zdot = drawcircle(posXch,posZch,circut1,vec=False)
        ax2.plot(xdot,zdot,'c--',linewidth=2,alpha=0.5)
        ydot,zdot = drawcircle(posYch,posZch,circut1,vec=False)
        ax3.plot(ydot,zdot,'c--',linewidth=2,alpha=0.5)

        xdot,ydot = drawcircle(posXch,posYch,curcut2,vec=False)
        ax1.plot(xdot,ydot,'c-.',linewidth=2,alpha=0.5)
        xdot,zdot = drawcircle(posXch,posZch,curcut2,vec=False)
        ax2.plot(xdot,zdot,'c-.',linewidth=2,alpha=0.5)
        ydot,zdot = drawcircle(posYch,posZch,curcut2,vec=False)
        ax3.plot(ydot,zdot,'c-.',linewidth=2,alpha=0.5)

        xdot,ydot = drawcircle(posXch,posYch,circut3,vec=False)
        ax1.plot(xdot,ydot,'c-',linewidth=2,alpha=0.5)
        xdot,zdot = drawcircle(posXch,posZch,circut3,vec=False)
        ax2.plot(xdot,zdot,'c-',linewidth=2,alpha=0.5)
        ydot,zdot = drawcircle(posYch,posZch,circut3,vec=False)
        ax3.plot(ydot,zdot,'c-',linewidth=2,alpha=0.5)
        
        scalevec.append(expfact)
        mvirvec.append(mvirch)
        mallvec.append(mallch)
        m200bvec.append(m200bch)
        vmaxvec.append(vmaxch)
        rvirvec.append(rvirch)
        vrmsvec.append(vrmsch)
        rsvec.append(rsch)
        xoffvec.append(xoffch)
        voffvec.append(voffch)
        virratvec.append(virratch)
        btavec.append(btach)
        ctavec.append(ctach)
        ctbvec.append(ctbch)
        spinvec.append(spinch)
        spinbvec.append(spinbch)
        pecvxvec.append(pecvxch)
        pecvyvec.append(pecvych)
        pecvzvec.append(pecvzch)
        jxvec.append(jxch)
        jyvec.append(jych)
        jzvec.append(jzch)
        Axvec.append(Axch)
        Ayvec.append(Aych)
        Azvec.append(Azch)
        npartvec.append(npartch)
        Tvec.append(Tch)

        ax4.plot(scalevec,np.log10(mvirvec),'r-',linewidth=2,label='virial')
        ax4.plot(scalevec,np.log10(mallvec),'g-',linewidth=2,label='+unbound')
        ax4.plot(scalevec,np.log10(m200bvec),'b-',linewidth=2,label='200*mean')
        ax5.plot(scalevec,rvirvec,'r-',linewidth=2,label='virial')
        ax5.plot(scalevec,rsvec,'b-',linewidth=2,label='scale')
        ax6.plot(scalevec,vmaxvec,'r-',linewidth=2)
        ax7.plot(scalevec,spinvec,'r-',linewidth=2,label='normal')
        ax7.plot(scalevec,spinbvec,'b-',linewidth=2,label='bullock')
        ax8.plot(scalevec,xoffvec,'r-',linewidth=2,label='position')
        ax8.plot(scalevec,voffvec,'b-',linewidth=2,label='velocity') 
        ax9.plot(scalevec,virratvec,'r-',linewidth=2)
        ax10.plot(scalevec,np.log10(npartvec),'r-',linewidth=2)
        ax11.plot(scalevec,Tvec,'r-',linewidth=2)
        ax12.plot(scalevec,btavec,'r-',linewidth=2,label='b/a')
        ax12.plot(scalevec,ctavec,'g-',linewidth=2,label='c/a') 
        ax12.plot(scalevec,ctbvec,'b-',linewidth=2,label='c/b')
        ax15.plot(scalevec,pecvxvec,'r-',linewidth=2,label='Vx') 
        ax15.plot(scalevec,pecvyvec,'g-',linewidth=2,label='Vy') 
        ax15.plot(scalevec,pecvzvec,'b-',linewidth=2,label='Vz') 
        ax14.plot(scalevec,jxvec,'r-',linewidth=2,label='Jx') 
        ax14.plot(scalevec,jyvec,'g-',linewidth=2,label='Jy') 
        ax14.plot(scalevec,jzvec,'b-',linewidth=2,label='Jz') 
        ax13.plot(scalevec,Axvec,'r-',linewidth=2,label='Ax') 
        ax13.plot(scalevec,Ayvec,'g-',linewidth=2,label='Ay') 
        ax13.plot(scalevec,Azvec,'b-',linewidth=2,label='Az') 
        
        if snap != endsnap:
            ax4.plot(expfact,np.log10(mvirch),'ro',markeredgewidth=0.0)
            ax4.plot(expfact,np.log10(mallch),'go',markeredgewidth=0.0)
            ax4.plot(expfact,np.log10(m200bch),'bo',markeredgewidth=0.0)
            ax5.plot(expfact,rvirch,'ro',markeredgewidth=0.0)
            ax5.plot(expfact,rsch,'bo',markeredgewidth=0.0)
            ax6.plot(expfact,vmaxch,'ro',markeredgewidth=0.0)
            ax7.plot(expfact,spinch,'ro',markeredgewidth=0.0)
            ax7.plot(expfact,spinbch,'bo',markeredgewidth=0.0)
            ax8.plot(expfact,xoffch,'ro',markeredgewidth=0.0)
            ax8.plot(expfact,voffch,'bo',markeredgewidth=0.0)
            ax10.plot(expfact,np.log10(npartch),'ro',markeredgewidth=0.0)
            ax9.plot(expfact,virratch,'ro',markeredgewidth=0.0)
            ax11.plot(expfact,Tch,'ro',markeredgewidth=0.0)
            ax12.plot(expfact,btach,'ro',markeredgewidth=0.0)
            ax12.plot(expfact,ctach,'go',markeredgewidth=0.0)
            ax12.plot(expfact,ctbch,'bo',markeredgewidth=0.0)
            ax15.plot(expfact,pecvxch,'ro',markeredgewidth=0.0)
            ax15.plot(expfact,pecvych,'go',markeredgewidth=0.0)
            ax15.plot(expfact,pecvzch,'bo',markeredgewidth=0.0)
            ax14.plot(expfact,jxch,'ro',markeredgewidth=0.0)
            ax14.plot(expfact,jych,'go',markeredgewidth=0.0)
            ax14.plot(expfact,jzch,'bo',markeredgewidth=0.0)
            ax13.plot(expfact,Axch,'ro',markeredgewidth=0.0)
            ax13.plot(expfact,Aych,'go',markeredgewidth=0.0)
            ax13.plot(expfact,Azch,'bo',markeredgewidth=0.0)
        
    
        ax11.text(0.5, 0.05,'oblate',
            horizontalalignment='center',
            verticalalignment='center',
            color='black',
            fontsize=12,
            transform = ax11.transAxes)

        ax11.text(0.5, 0.95,'prolate',
            horizontalalignment='center',
            verticalalignment='center',
            color='black',
            fontsize=12,
            transform = ax11.transAxes)

        xtrans = 0.5
        ytrans = 0.5
        ztrans = 0.5
        extra = 1.5*rvirch/1000
        addon = extra/boxwidth 

        placetext(ax1,xtrans+addon,ytrans,mstring,'bold',11)
        placetext(ax2,xtrans+addon,ztrans,mstring,'bold',11)
        placetext(ax3,ytrans+addon,ztrans,mstring,'bold',11)

    ax2.set_title(titlestr)

    ax1.set_xlabel('x-pos [Mpc/h]')
    ax1.set_ylabel('y-pos [Mpc/h]')
    ax2.set_xlabel('x-pos [Mpc/h]')
    ax2.set_ylabel('z-pos [Mpc/h]')
    ax3.set_xlabel('y-pos [Mpc/h]')
    ax3.set_ylabel('z-pos [Mpc/h]')


    ax4.set_ylabel('Mvir [Msol/h]')
    ax5.set_ylabel('Radius [kpc]')
    ax6.set_ylabel('Vmax [km/s]')
    ax7.set_ylabel('spin')
    ax8.set_ylabel('offset')
    ax9.set_ylabel('virial ratio')
    ax10.set_ylabel('number of particles')
    ax11.set_ylabel('triaxiality parameter')
    ax12.set_ylabel('axis ratios')
    ax15.set_ylabel('peculiar Vx, Vy, Vz [km/s]')
    ax14.set_ylabel('Jx, Jy, Jz')
    ax13.set_ylabel('ellipticity axis')

    ax4.set_xlabel('scale factor')
    ax5.set_xlabel('scale factor')
    ax6.set_xlabel('scale factor')
    ax7.set_xlabel('scale factor')
    ax8.set_xlabel('scale factor')
    ax9.set_xlabel('scale factor')
    ax10.set_xlabel('scale factor')
    ax11.set_xlabel('scale factor')
    ax12.set_xlabel('scale factor')
    ax13.set_xlabel('scale factor')
    ax14.set_xlabel('scale factor')
    ax15.set_xlabel('scale factor')

    plotlegend(ax4,legendfontsize)
    plotlegend(ax5,legendfontsize)
    plotlegend(ax7,legendfontsize)
    plotlegend(ax8,legendfontsize)
    plotlegend(ax12,legendfontsize)
    plotlegend(ax13,legendfontsize)
    plotlegend(ax14,legendfontsize)
    plotlegend(ax15,legendfontsize)

    ax1.set_xlim([xcen - boxwidth/2., xcen + boxwidth/2.])
    ax1.set_ylim([ycen - boxwidth/2., ycen + boxwidth/2.])
    ax2.set_xlim([xcen - boxwidth/2., xcen + boxwidth/2.])
    ax2.set_ylim([zcen - boxwidth/2., zcen + boxwidth/2.])
    ax3.set_xlim([ycen - boxwidth/2., ycen + boxwidth/2.])
    ax3.set_ylim([zcen - boxwidth/2., zcen + boxwidth/2.])

    ax4.set_xlim([0,1])
    ax5.set_xlim([0,1])
    ax6.set_xlim([0,1])
    ax7.set_xlim([0,1])
    ax8.set_xlim([0,1])
    ax9.set_xlim([0,1])
    ax11.set_ylim([0,1])
    ax10.set_xlim([0,1])
    ax11.set_xlim([0,1])
    ax12.set_xlim([0,1])
    ax13.set_xlim([0,1])
    ax14.set_xlim([0,1])
    ax15.set_xlim([0,1])

    #fig.savefig(figdir + 'snap' + str(snap) + '-xyzhalopos-haloparam.png',bbox_inches='tight')
    fig.savefig(figdir + 'TESTsnap' + str(snap) + '.png',bbox_inches='tight')
    plt.close(fig)


#plt.show()
