import readsnapshots.readsnapHDF5 as rsHD
import readsnapshots.readsnap as rs
import readhalos.readsubf
import readsnapshots.readids
import numpy as np
import pylab as plt
import sys
from random import randint
from matplotlib import *
from grifflib import *
import mergertrees.GregNewMTCatalogue as MT

# DEFINE PARAMS
#candidatelist = [80609,158376]
candidatelist = [208737]
#190897]
#,208737,140666,28221,147419,28188,147273,78411,131988,19910]
ticksize = 11
hubble = 0.6711

legendfontsize = 10
axislabelfontsize = 14
#Alternatives
#121761,43307

#2.9E12 Msol
#78442
#78442
#78411
#131988

#2E12 Msol
#104939
#233429
#198870
#60552
#164560

#7E11
#recent mergers
#132399
#217067
#51346

#violent history
#225903
#147515
#259960

# SET PATHS
#halopath = '/n/home01/bgriffen/data/caterpillar/parent/RockstarData'  # Odyssey
halopath = '/spacebase/data/AnnaGroup/caterpillar/parent/RockstarData' # Spacebase
#candfile = '/n/home01/bgriffen/data/caterpillar/halos/candidates.dat' #odyssey
candfile = '/spacebase/data/bgriffen/code/caterpillar/candidates.dat' #spacebase
treefile = halopath + '/trees/tree.bin'
indexfile = halopath + '/trees/treeindex.csv'

fig1 = plt.figure(figsize=(22.0,12.0))
ax1 = fig1.add_subplot(3,5,1)
ax2 = fig1.add_subplot(3,5,2)
ax3 = fig1.add_subplot(3,5,3)
ax4 = fig1.add_subplot(3,5,4)
ax5 = fig1.add_subplot(3,5,5)
ax6 = fig1.add_subplot(3,5,6)
ax7 = fig1.add_subplot(3,5,7)
ax8 = fig1.add_subplot(3,5,8)
ax9 = fig1.add_subplot(3,5,9)
ax10 = fig1.add_subplot(3,5,10)
ax11 = fig1.add_subplot(3,5,11)
ax12 = fig1.add_subplot(3,5,12)
ax13 = fig1.add_subplot(3,5,13)
ax14 = fig1.add_subplot(3,5,14)
ax15 = fig1.add_subplot(3,5,15)

plt.subplots_adjust(hspace=0.1,wspace=0.3)


ncols = 0
icand = 0
for haloid in candidatelist:
    haloid = int(haloid)
    halotree = MT.NewMTCatalogue(treefile,indexfile,hostid=haloid)
    tree = halotree.Trees[haloid]
    mainbranch = tree.getMainBranch(0)
    scale = mainbranch['scale']
    rvir = mainbranch['rvir']
    posX = mainbranch['posX']
    posY = mainbranch['posY']
    posZ = mainbranch['posZ']
    mvir = mainbranch['mvir']
    mall = mainbranch['m200c_all']
    m200b = mainbranch['m200b']
    vmax = mainbranch['vmax']
    vrms = mainbranch['vrms']
    rs = mainbranch['rs']
    xoff = mainbranch['xoff']
    voff = mainbranch['voff']
    virrat = mainbranch['T/|U|']
    bta = mainbranch['b_to_a']
    cta = mainbranch['c_to_a']
    spin = mainbranch['spin']
    spinb = mainbranch['spin_bullock']
    pecvx = mainbranch['pecVX']
    pecvy = mainbranch['pecVY']
    pecvz = mainbranch['pecVZ']
    jx = mainbranch['Jx']
    jy = mainbranch['Jy']
    jz = mainbranch['Jz']
    Ax = mainbranch['A[x]']
    Ay = mainbranch['A[y]']
    Az = mainbranch['A[z]']
    ctb = bta*(1/cta)
    Tch = calcT(ctb)

    npart = mvir/6.554894E6

    normmvir = mvir/mvir[0]
    normvmax = vmax**2/vmax[0]**2

    icand += 1
    print "---------------------------------------"
    print "Candidate:",icand,"Rockstar ID:",haloid
    print "---------------------------------------"
    print "        x-pos:",'{:.2f}'.format(posX[0]), "   \ [Mpc/h]"
    print "        y-pos:",'{:.2f}'.format(posY[0]), "   \ [Mpc/h]"
    print "        z-pos:",'{:.2f}'.format(posZ[0]), "   \ [Mpc/h]"
    print "  virial mass:",'{0:.2e}'.format(mvir[0]),"\ [Msol]"
    print "virial radius:",'{:.2f}'.format(rvir[0]),"  \ [kpc]"

    plotquantl(ax1,scale,np.log10(mvir),'virial')
    plotquantl(ax1,scale,np.log10(mall),'+unbound')
    plotquantl(ax1,scale,np.log10(m200b),'m200')
    plotquant(ax2,scale,normmvir)
    plotquantl(ax3,scale,vmax,'max')
    plotquantl(ax3,scale,vrms,'rms')
    plotquant(ax4,scale,normvmax)
    plotquantl(ax5,scale,rvir,'virial')  
    plotquantl(ax5,scale,rs,'scale')
    plotquantl(ax6,scale,spin,'normal')
    plotquantl(ax6,scale,spinb,'bullock')
    plotquantl(ax7,scale,xoff,'position')
    plotquantl(ax7,scale,voff,'velocity')
    plotquant(ax8,scale,virrat)
    plotquantl(ax9,scale,np.log10(npart),haloid)
    plotquantl(ax10,scale,bta,'b/a')
    plotquantl(ax10,scale,cta,'c/a')
    plotquant(ax11,scale,Tch)
    plotquantl(ax12,scale,pecvx,'Vx')
    plotquantl(ax12,scale,pecvy,'Vy')
    plotquantl(ax12,scale,pecvz,'Vz')
    plotquantl(ax13,scale,jx,'Jx')
    plotquantl(ax13,scale,jy,'Jy')
    plotquantl(ax13,scale,jz,'Jz')
    plotquantl(ax14,scale,Ax,'Ax')
    plotquantl(ax14,scale,Ay,'Ay')
    plotquantl(ax14,scale,Az,'Az')



    #ax1.plot(scale,np.log10(mvir),linestyle='-',linewidth=2,label=str(haloid))
    #ax2.plot(scale,np.log10(normmvir),linestyle='-',linewidth=2,label=str(haloid))
    #ax3.plot(scale,spin ,linestyle='-',linewidth=2,label=str(haloid))
    #ax4.plot(scale,rvir,linestyle='-',linewidth=2,label=str(haloid))
    #ax5.plot(scale,vmax,linestyle='-',linewidth=2,label=str(haloid))
    #ax6.plot(scale,normvmax,linestyle='-',linewidth=2,label=str(haloid))
    #ax7.plot(scale,virialratio,linestyle='-',linewidth=2)
    #ax8.plot(scale,xoff,linestyle='-',linewidth=2)
    
    
    if icand % 3 != 0:
        ncols += 1
print "---------------------------------------"

handles, labels = ax9.get_legend_handles_labels()
fig1.legend(handles, labels, 'upper center',prop={'size':10},ncol=ncols)
                
new_tick_locations = np.array([.2, .5, .9])

ax11.text(0.5, 0.05,'oblate',
    horizontalalignment='center',
    verticalalignment='center',
    color='black',
    fontsize=11,
    transform = ax11.transAxes)
ax11.text(0.5, 0.95,'prolate',
    horizontalalignment='center',
    verticalalignment='center',
    color='black',
    fontsize=11,
    transform = ax11.transAxes)

ax15.text(0.5, 0.5,'spare',
    horizontalalignment='center',
    verticalalignment='center',
    color='black',
    fontsize=11,
    transform = ax15.transAxes)

ax1.set_xlim([0,1])
ax2.set_xlim([0,1])
ax3.set_xlim([0,1])
ax4.set_xlim([0,1])
ax5.set_xlim([0,1])
ax6.set_xlim([0,1])
ax7.set_xlim([0,1])
ax8.set_xlim([0,1])
ax9.set_xlim([0,1])
ax10.set_xlim([0,1])
ax11.set_xlim([0,1])
ax12.set_xlim([0,1])
ax13.set_xlim([0,1])
ax14.set_xlim([0,1])
ax15.set_xlim([0,1])
ax11.set_ylim([0,1])

plotlegend(ax1,legendfontsize,location=4)
plotlegend(ax3,legendfontsize)
plotlegend(ax5,legendfontsize)
plotlegend(ax6,legendfontsize)
plotlegend(ax7,legendfontsize)
plotlegend(ax10,legendfontsize,location=4)
plotlegend(ax12,legendfontsize)
plotlegend(ax13,legendfontsize)
plotlegend(ax14,legendfontsize)

new_tick_locations = np.array([0.,.2, .4, 0.6, 0.8, 1.])

ax1top = ax1.twiny()
ax1top.set_xticklabels(tick_function(new_tick_locations))

ax2top = ax2.twiny()
ax2top.set_xticklabels(tick_function(new_tick_locations))

ax3top = ax3.twiny()
ax3top.set_xticklabels(tick_function(new_tick_locations))

ax4top = ax4.twiny()
ax4top.set_xticklabels(tick_function(new_tick_locations))

ax5top = ax5.twiny()
ax5top.set_xticklabels(tick_function(new_tick_locations))


ax1top.set_xlabel(r'$\mathrm{redshift}$',size=14)
ax2top.set_xlabel(r'$\mathrm{redshift}$',size=14)
ax3top.set_xlabel(r'$\mathrm{redshift}$',size=14)
ax4top.set_xlabel(r'$\mathrm{redshift}$',size=14)
ax5top.set_xlabel(r'$\mathrm{redshift}$',size=14)

ax1.set_xticks([])
ax2.set_xticks([])
ax3.set_xticks([])
ax4.set_xticks([])
ax5.set_xticks([])
ax6.set_xticks([])
ax7.set_xticks([])
ax8.set_xticks([])
ax9.set_xticks([])
ax10.set_xticks([])


ax1.set_ylabel(r'$\mathrm{log_{10}\ M(z)\ [M_\odot]}$',size=axislabelfontsize)
ax2.set_ylabel(r'$\mathrm{M_v(z)/M_v(z=0)}$',size=axislabelfontsize)
ax3.set_ylabel(r'$\mathrm{V(z)\ [km/s]}$',size=axislabelfontsize)
ax4.set_ylabel(r'$\mathrm{V_{max}(z)^2/V_{max}(z=0)^2}$',size=axislabelfontsize)
ax5.set_ylabel(r'$\mathrm{radius\ [kpc]}$',size=axislabelfontsize)
ax6.set_ylabel(r'$\mathrm{spin}$',size=axislabelfontsize)
ax7.set_ylabel(r'$\mathrm{offset}$',size=axislabelfontsize)
ax8.set_ylabel(r'$\mathrm{virial\ ratio}$',size=axislabelfontsize)
ax9.set_ylabel(r'$\mathrm{log_{10}\ number\ of\ particles}$',size=axislabelfontsize)
ax10.set_ylabel(r'$\mathrm{axis\ ratios}$',size=axislabelfontsize)
ax11.set_ylabel(r'$\mathrm{triaxiality\ parameter}$',size=axislabelfontsize)
ax12.set_ylabel(r'$\mathrm{peculiar\ V_x,\ V_y,\ V_z\ [km/s]}$',size=axislabelfontsize)
ax13.set_ylabel(r'$\mathrm{J_x, J_y, J_z}$',size=axislabelfontsize)
ax14.set_ylabel(r'$\mathrm{ellipticity\ axis}$',size=axislabelfontsize)

ax11.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
ax12.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
ax13.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
ax14.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
ax15.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)

#ax1.tick_params(axis='both', which='major', labelsize=ticksize)
#ax2.tick_params(axis='both', which='major', labelsize=ticksize)
#ax3.tick_params(axis='both', which='major', labelsize=ticksize)
#ax4.tick_params(axis='both', which='major', labelsize=ticksize)
#ax5.tick_params(axis='both', which='major', labelsize=ticksize)
#ax6.tick_params(axis='both', which='major', labelsize=ticksize)
#ax7.tick_params(axis='both', which='major', labelsize=ticksize)
#ax8.tick_params(axis='both', which='major', labelsize=ticksize)
#ax1top.tick_params(axis='both', which='major', labelsize=ticksize)
#ax2top.tick_params(axis='both', which='major', labelsize=ticksize)
#ax3top.tick_params(axis='both', which='major', labelsize=ticksize)
#ax4top.tick_params(axis='both', which='major', labelsize=ticksize)

#fig1.savefig('./mvirvmaxhistory.png',bbox_inches='tight')

plt.show()
