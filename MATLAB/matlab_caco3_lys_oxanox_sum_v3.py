# -*- coding: utf-8 -*-
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import glob
from matplotlib.ticker import MultipleLocator
"""
caco3 profiles 
"""

################ Setting up working directory ###################

#Workdir = "C:/Users/YK/Desktop/Sed_res/"
##Workdir = "C:/Users/YK/Desktop/HT_res/permgrad-20181104/"
Workdir = "/home/domhu/Documents/GitHub/caco3/MATLAB/2006_lysocline/"

org = 'oxonly'
org = 'oxanox'

# org_end = '_oxonly'
org_end = '_oxanox'

addname = '-zox_corr-dt1e8_sp12'
##addname='-turbo2'

mixing = ''
##mixing = 'turbo2'

# Workdir = 'C:/Users/YK/Desktop/Sed_res/'
Workdir = "/home/domhu/Documents/GitHub/caco3/MATLAB/2006_lysocline/"
Workdir += 'test-translabs/res/multi/'
if not mixing=='':
    Workdir += org+'-'+mixing
else:
    Workdir += org
if not addname =='':
    Workdir += addname
    
Workdir += '/'

# Or you can just specify the working directory here
#  (where result files are restored)
# Workdir += './'
Workdir = "/home/domhu/Documents/GitHub/caco3/MATLAB/2006_lysocline/"

################ Setting up working directory ###################

rr = '0_0E0'

rrlist = ['0.00','0.50','0.67','1.00','1.50']

ddd = 6 + len(mixing)+1

flx = [6,12,18,24,30,36,42,48,54,60]

title='CaCO' +r"$\mathregular{_3}$"+' rain\n'\
                + '('+r"$\mathregular{\mu}$"+'mol cm' \
                +r"$\mathregular{^{-2}}$"+' yr' \
                +r"$\mathregular{^{-1}}$"+')'



plt.rcParams['font.family'] = 'Arial' 
plt.rcParams['font.size'] = 20

linewidth = 1.5

plt.rcParams['axes.linewidth'] = linewidth

plt.rcParams['xtick.major.width'] = linewidth
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.width'] = linewidth
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.pad'] = 8

plt.rcParams['ytick.major.width'] = linewidth
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.minor.width'] = linewidth
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.major.pad'] = 8

# plt.rcParams['axes.labelpad'] = 8

plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

plt.tick_params(top=True)
plt.tick_params(right=True)

nx = 5
ny=2

fig = plt.figure(figsize=(7,10))

# file size determined by 5 rain ratios, 10 cc flxes, 25 water depths
#    and 10 and 5 data for CaCO3 wt% and bur flx, respectively
#   (thus, 10 and 5 may change depending on the number of data
#       you restore: previously these are 7 and 2)
##data_lys = np.zeros((5,10,25,10),dtype=np.float)
##data_bur = np.zeros((5,10,25,5),dtype=np.float)
data_lys = np.zeros((5,10,25,7),dtype=np.float)
data_bur = np.zeros((5,10,25,2),dtype=np.float)

cmap = plt.cm.get_cmap('jet')
##cmap = get_parula()

# First read data 
for j in range(5):
    rr = rrlist[j]
    filelist = glob.glob(Workdir+'lys_sense_cc-*_rr-'+rr+org_end+'.txt')
    filelist[0], filelist[1:] = filelist[-1], filelist[:-1]

    datatmp = np.loadtxt(filelist[0])
    data = np.zeros((len(filelist),datatmp.shape[0],datatmp.shape[1])\
                    ,dtype=np.float)

    for i in range(len(filelist)):
        data_lys[j,i,:,:]=np.sort(np.loadtxt(filelist[i]),axis=0)

    filelist = glob.glob(Workdir+'ccbur_sense_cc-*_rr-'+rr+org_end+'.txt')
    filelist[0], filelist[1:] = filelist[-1], filelist[:-1]
    
    for i in range(len(filelist)):
        data_bur[j,i,:,:]=np.sort(np.loadtxt(filelist[i]),axis=0)

# Then plot data
for j in range(5):
    ax = plt.subplot2grid((ny,nx), (0,j))

    color=cmap(np.linspace(0,1,10))

    for i in range(10):
##        if org=='ox':label = filelist[i][63+ddd:69+ddd]
##        if org=='oxanox':label = filelist[i][67+ddd:73+ddd]
        label=str(flx[i])
        plt.plot(data_lys[j,i,:,5],data_lys[j,i,:,0],'-x',c=color[i]
                 ,label=label
                 )
##    if rr == '0_0E0':plt.legend(facecolor='None',edgecolor='None',title=title)
    ax = plt.gca()
    ax.set_xlim(0,90)
    ax.set_ylim(-40,40)
    ax.set_xticks([0,45,90])
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_minor_locator(MultipleLocator(15))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    
    if j!=0:
        ax.set_xticks([45,90])
        ax.set_yticklabels([])
    fig.subplots_adjust(bottom=0.2)

##    if org=='ox':label = filelist[0][70+ddd:76+ddd]
##    if org=='oxanox':label = filelist[0][74+ddd:80+ddd]

    ##outfilename = Workdir+'lys_'+label+'.svg'
    ##plt.savefig(outfilename, transparent=True)
    ##subprocess.call('"C:\Program Files\Inkscape\inkscape.exe" -z -f ' \
    ##                + outfilename + ' --export-emf '+outfilename+\
    ##                '.emf',shell=True)
    ##plt.show()
    ##plt.clf()
    ##plt.close()

    print data_lys[j,:,:,2].max(),data_lys[j,:,:,2].min()
    print data_lys[j,:,:,4].max(),data_lys[j,:,:,4].min()
    print data_lys[j,:,:,6].max(),data_lys[j,:,:,6].min()


    ax2 = plt.subplot2grid((ny,nx), (1,j))
    color=cmap(np.linspace(0,1,10))

    for i in range(10):
##        if org=='ox':label = filelist[i][65+ddd:71+ddd]
##        if org=='oxanox':label = filelist[i][69+ddd:75+ddd]
        label=str(flx[i])
        plt.plot(data_bur[j,i,:,1],data_bur[j,i,:,0],'-x',c=color[i]
                 ,label=label
                 )
    ##plt.legend(facecolor='None',edgecolor='None')
    ax2 = plt.gca()
    ax2.set_xlim(0,40)
    ax2.set_ylim(-40,40)
    ax2.set_xticks([0,20,40])
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.xaxis.set_minor_locator(MultipleLocator(10))
    ax2.yaxis.set_minor_locator(MultipleLocator(5))
    
    if j!=0:
        ax2.set_yticklabels([])
        ax2.set_xticks([20,40])


##fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(left=0.25,bottom=0.1,wspace=0.06,hspace=0.3)

##if org=='ox':label = filelist[0][72+ddd:78+ddd]
##if org=='oxanox':label = filelist[0][76+ddd:82+ddd]

##outfilename = Workdir+'lysbur-'+org+'.svg'
outfilename = Workdir+'lysbur'+org_end+'.svg'
plt.savefig(outfilename, transparent=True)
# subprocess.call('"C:\Program Files\Inkscape\inkscape.exe" -z -f ' \
#                + outfilename + ' --export-emf '+outfilename+\
#                '.emf',shell=True)
plt.show()
plt.clf()
plt.close()
