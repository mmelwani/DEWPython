"""
    profile_plots -- Authors: Seda & Emre Isik, 2024
    Plots of the depth profile of net affinity of reactions around an input species in 
    a reaction network. See Isik, S. et al. 2025 ACS Earth & Space Chemistry (in press), 
    Figures 9-13, including TCA cycle ring plot indicating relative affinity fractions
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import glob

R_earth = 6378.*1e3 # in meters
R_eur, R_enc = 1561., 252.1 # in km
R_gan, R_tit = 2634.1, 2574.73
R_planet = R_enc
planet = 'Enceladus'

path_all = 'Species_Profiles/'+planet+'-ph11/'
path_TCA = 'Species_Profiles/'+planet+'-ph11/TCA-only/'
profiles_TCA = sorted(glob.glob(os.path.join(path_TCA, "*.csv")))
profiles_all = sorted(glob.glob(os.path.join(path_all, "*.csv")))
cols_TCA = np.array(['C2','C1','C4','C9','C8','C3','C5','C6'])
cols_all = np.array(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
labels_TCA = np.array(['Citrate', 'Cis-aconitate', 'Isocitrate', 'a-Ketoglutarate', 'Succinate',
                       'Fumarate', 'Malate', 'Oxaloacetate'])
labels_all = np.array(['Acetate', 'Cis-aconitate', 'Citrate', 'Fumarate',
                       'Isocitrate', 'Malate', 'Oxaloacetate', 'Pyruvate', 'Succinate',
                       r'$\alpha$-ketoglutarate'])

print(profiles_TCA)
print(profiles_all)

plt.rcParams.update({'font.size': 14})
# Change the following according to the moon in consideration
radius_limit_sup = 0.69373117  
# Europa 0.93435835
# Titan 0.73752561 dew_max
# Ganymede 0.69373117 dew_max_radius, above which we extrapolate.
# Enceladus 0.77021467 # ??? -- no radius_limit_sup in Enceladus! 

fig, ax = plt.subplots(figsize=(12,8))
if (planet != 'LostCity'): ax2 = ax.twiny()
i0 = 0
for i in range(len(profiles_all)):
    df = pd.read_csv(profiles_all[i])
    df_dew = df[df.model == 'dew']
    df_sup = df[df.model == 'sup']
    ax.plot(df_sup.R, df_sup.aff_net/1e3, c=cols_all[i],    # /1e3 to convert to kJ
            alpha=.5, linewidth=3, label=labels_all[i])
    ax.plot(df_dew.R, df_dew.aff_net/1e3, c=cols_all[i],
            alpha=.5, linewidth=3)
    if (planet != 'LostCity' and planet != 'Enceladus'):
        ax.axvline(x=radius_limit_sup, color='gray', alpha=.3, linestyle=':')
    ax.axhline(y=0, color='gray', alpha=.3, linestyle='--')

if (planet != 'LostCity'): ax2.plot(df_sup.R*R_tit, 0*df_sup.aff_net, c=cols_all[i],
                                    alpha=0, linewidth=3)
ax.legend(bbox_to_anchor=(1.2, 1.05), loc='upper left', fontsize=10, frameon=False)
#ax.set_title(planet, loc='left')
ax.set_title('Enceladus-specific (pH = 11)', loc='left')
#ax.set_title('Enceladus (-specific activities)', loc='left')

if (planet != 'LostCity'):
    ax.set_xlabel('$r/R_{\mathrm{surface}}$', fontsize=16)
    if (planet != 'LostCity'): ax2.set_xlabel('$r~({\mathrm{km}})$')
else: ax.set_xlabel('Depth below ocean floor (m)', fontsize=16)
ax.set_ylabel(r'$A_{\mathrm{net}}~~({\mathrm{kJ/mol}})$', fontsize=16)
plt.tight_layout()
#plt.show()
#fig.savefig('Europa_profiles.pdf')

# Radius-integrated species abundances.
# First, let's integrate along the entire SUPCRT array, first as a BAR CHART:
#sums = np.zeros(10)
sums = np.zeros(len(profiles_all))
i=0
for i in range(len(profiles_all)):
    df = pd.read_csv(profiles_all[i])
    df_dew = df[df.model == 'dew']
    df_sup = df[df.model == 'sup']
    r = np.flip(df_sup.R) * R_planet
    rr = (df_sup.R*1e-3) # be careful about the consistency of this with the qty in integral!
    H = np.max(rr) - np.min(rr) # height of the Lost-City dome
    print(rr)
    print('H=', H, 'km')
    print(np.min(rr), np.max(rr))
    if (planet == 'LostCity'):
#        sums[i] = np.trapz(4.*np.pi*np.flip(df_sup.eqnet)*rr**2, x=rr)
        sums[i] = np.trapz(4.*np.pi*np.flip(df_sup.aff_net)*r**2, x=r) / np.trapz(4.*np.pi*r**2, x=r)
#        sums[i] = np.trapz(df_sup.aff_net * (H**2 - rr**2), x=rr)
    else:
#        sums[i] = np.trapz(np.flip(df_sup.eqnet), x=np.flip(df_sup.R))
        sums[i] = np.trapz(4.*np.pi*np.flip(df_sup.aff_net)*r**2, x=r) / np.trapz(4.*np.pi*r**2, x=r)
    print(i, labels_all[i], sums[i])
#    print(df_sup.aff_net)
    print('Sums of sums as an overall energetic')
#    sums[i] = np.trapz(np.flip(df_sup.aff_net), x=np.flip(df_sup.R))
    print(labels_all[i], sums[i], np.mean(df_sup.aff_net))
#    print(df_sup.aff_net)

#fig2, ax2 = plt.subplots()
axins = inset_axes(ax, width="100%", height="100%",
                   bbox_to_anchor=(1.05, .35, .35, .3),
                   bbox_transform=ax.transAxes, loc=2, borderpad=0)
axins.tick_params(left=False, right=True, labelleft=False, labelright=True, labelbottom=False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
axins.set_xticks([])
axins.bar(labels_all, sums/1e3, color=cols_all, alpha=.5) # 1e3 kJ conversion
axins.set_ylabel(r'$\bar{A}_{\mathrm{SUP}}~(\mathrm{kJ/mol})$')
axins.yaxis.set_label_position("right")
#axins.set_title(r'$A_{\mathrm{SUP}}~(\mathrm{km}^3)$')
#axins.set_title(r'$\int_{\mathrm{SUP}}\log K_{\mathrm{net}}dr$')
axins.axhline(linewidth=.5)
#axins.set_title(r'$\int_{\mathrm{SUP}}A_{\mathrm{net}}dr$')

###### RING PLOT for TCA ######
sums_TCA = np.zeros(len(labels_TCA))
size = 0.3
for i in range(len(labels_TCA)):
    df = pd.read_csv(profiles_TCA[i])
    df_dew = df[df.model == 'dew']
    df_sup = df[df.model == 'sup']
    print(i), labels_all, labels_TCA
    sums_TCA[i] = np.trapz(np.flip(df_sup.aff_net), x=np.flip(df_sup.R))

print(sums_TCA)
#fig3, ax3 = plt.subplots()
axins2 = inset_axes(ax, width="100%", height="100%",
                   bbox_to_anchor=(1.05, .001, .4, .3),
                   bbox_transform=ax.transAxes, loc=2, borderpad=0)
axins2.tick_params(left=False, right=False, labelleft=False, labelright=False, labelbottom=False)
labels_TCA_short = np.array(['Cit', 'CisA', 'I', r'$\alpha$', 'S', 'F', 'M', 'O'])
sums1 = np.array([sums[2],sums[1],sums[4],sums[9],sums[8],sums[3],sums[5],sums[6]])
#breakpoint()
#axins2.pie(3 + sums_TCA/1e4, radius=1, colors=cols_TCA, labels=labels_TCA_short,
#offset = 1 # 1 Eur - 8 Enc
offset = np.abs(1.25*np.min(sums1))
axins2.pie(offset + sums1, radius=1, colors=cols_TCA, labels=labels_TCA_short,
           counterclock=False, wedgeprops=dict(width=size, edgecolor='w', alpha=0.5))

axins2.set(aspect="equal")
#patches.wedge.set_alpha(.5)
#axins2.set_title('Abiotic TCA cycle', loc='right')
axins2.set_xlabel('TCA cycle')
plt.savefig('out/'+planet+'-ph11/'+planet+'_summary.pdf')
#plt.savefig('out/'+planet+'/'+planet+'_specific_summary.pdf')
#plt.show()