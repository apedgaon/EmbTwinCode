import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams # The rcPArams method to customize the legend font size
rcParams['legend.fontsize'] = 20
rcParams['pdf.fonttype'] = 42  # Avoid Font type 3
rcParams['ps.fonttype'] = 42   # Avoid Font tyep 3
rcParams['axes.unicode_minus'] = False # fix the 'minus' sign display issue
plt.rc('xtick',labelsize=18)   # Change label tick font size
plt.rc('ytick',labelsize=18)
import math
import sys

from matplotlib.ticker import FormatStrFormatter

elmt_number = sys.argv[1]
elmt_file = 'elmt_'+elmt_number+'.hist'


elmt_SDV = np.loadtxt(elmt_file)
time_step = 1.0E-3

time = []

M_SDV1 = []
SB_SDV5 = []

M_SDV2 = []
SB_SDV6 = []

M_SDV3 = []
SB_SDV7 = []

M_SDV4 = []
SB_SDV8 = []

i = 0
for line in elmt_SDV:
    i = i + 1
    if line[-1] == 1.0:
        break


for line in elmt_SDV:
    time.append(int(line[0])*time_step)

    M_SDV1.append(line[1])
    SB_SDV5.append(line[5])

    M_SDV2.append(line[2])
    SB_SDV6.append(line[6])

    M_SDV3.append(line[3])
    SB_SDV7.append(line[7])

    M_SDV4.append(line[4])
    SB_SDV8.append(line[8])


fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(time, M_SDV1, linewidth = 2.0,color='blue',linestyle='--',label='Matrix')
ax.plot(time[i-1:], SB_SDV5[i-1:], linewidth = 2.0,color='red',linestyle='-',label='Shear band')

leg = ax.legend(loc='upper left',frameon=True)
for label in leg.get_lines():
	label.set_linewidth(2.0)  # the legend line width
plt.grid(True)
plt.xlim( 0, 60 )
plt.xlabel('Time ($10^{-6}$s)',fontsize=24)
plt.ylabel('Von Mises stress (MPa)',fontsize=24)
#plt.show()
ax.plot([18], [1125], 'o')
ax.annotate('Onset of shear localization', fontsize=24,  style='italic', xy=(18, 1125), xytext=(5, 600),
            arrowprops=dict(facecolor='black', shrink=0.05))
plt.tight_layout()
plt.savefig('Top_hat_elmt_'+elmt_number+'_Matrix_VS_SB_VonMisesStress_slides.eps', format='eps', dpi=1000)



fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(time, M_SDV2, linewidth = 2.0,color='blue',linestyle='--',label='Matrix')
ax.plot(time[i-1:], SB_SDV6[i-1:], linewidth = 2.0,color='red',linestyle='-',label='Shear band')

leg = ax.legend(loc='upper left',frameon=True)
for label in leg.get_lines():
	label.set_linewidth(2.0)  # the legend line width
plt.grid(True)
plt.xlim( 0, 60 )
plt.xlabel('Time ($10^{-6}$s)',fontsize=24)
plt.ylabel('Equivalent plastic strain',fontsize=24)

ax.plot([19], [0.3], 'o')
ax.annotate('Onset of shear localization', fontsize=24,  style='italic', xy=(19, 0.3), xytext=(10, 4),
            arrowprops=dict(facecolor='black', shrink=0.05))
#plt.show()
plt.tight_layout()
plt.savefig('Top_hat_elmt_'+elmt_number+'_Matrix_VS_SB_PlasticStrain_slides.eps', format='eps', dpi=1000)



fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(time, M_SDV3, linewidth = 2.0,color='blue',linestyle='--',label='Matrix')
ax.plot(time[i-1:], SB_SDV7[i-1:], linewidth = 2.0,color='red',linestyle='-',label='Shear band')

leg = ax.legend(loc='upper left',frameon=True)
for label in leg.get_lines():
	label.set_linewidth(2.0)  # the legend line width
plt.grid(True)
plt.xlim( 0, 60 )
plt.ylabel('Temperature (K)',fontsize=24)
plt.xlabel('Time ($10^{-6}$s)',fontsize=24)
#plt.show()
plt.tight_layout()
plt.savefig('Top_hat_elmt_'+elmt_number+'_Matrix_VS_SB_Temperature_slides.eps', format='eps', dpi=1000)



fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(time, M_SDV4, linewidth = 2.0,color='blue',linestyle='--',label='Matrix')
ax.plot(time[i-1:], SB_SDV8[i-1:], linewidth = 2.0,color='red',linestyle='-',label='Shear band')

leg = ax.legend(loc='upper left',frameon=True)
for label in leg.get_lines():
	label.set_linewidth(2.0)  # the legend line width
plt.grid(True)
plt.xlim( 0, 60 )
plt.xlabel('Time ($10^{-6}$s)',fontsize=24)
plt.ylim(-200,1600)
plt.ylabel('Struct. Evol. thre. stress (MPa)',fontsize=24)
#plt.show()
plt.tight_layout()
plt.savefig('Top_hat_elmt_'+elmt_number+'_Matrix_VS_SB_MTSStress_slides.eps', format='eps', dpi=1000)

