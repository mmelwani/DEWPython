"""
    affinity -- calculates chemical affinity from given activities of species entering a reaction
    AUTHORS: Seda & Emre Isik 2024 (Isik, S. et al. 2025 ACS Earth & Space Chem.)

    INPUTS: - chemical activities of the species of the reaction (read from text file)
            - equilibrium constant of the reaction = f(P,T) (calc. by DEW or SUPCRT, using DEWPython)
            - Geotherms of planetary interior from PlanetProfile - only to use PT curves on the affinity 
            plot as a function of P,T
    OUTPUTS: - csv files listing the chemical affinity of the reaction with the ambiend P and T. 
            - affinity plots with PT curves of planetary bodies.
"""

import warnings
import glob, os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def read_activities(file_path, reaction):
    df = pd.read_csv(file_path, header=None, names=['Species', 'Activity'], 
                     delim_whitespace=True, dtype={'Activity': float})

    activities = df.set_index('Species')['Activity'].to_dict()
    relevant_activities = {species: activities[species] for species in reaction if species in activities}
    return relevant_activities

def read_logk(file_path):
    df = pd.read_csv(file_path)
    temperature = df['Temperature'].to_numpy()
    pressure = df['Pressure'].to_numpy()
    logk = df['LogK'].to_numpy()
    delG = df['delG'].to_numpy()
    if np.min(pressure >= 1000):
        temperature = temperature[pressure >= 4000]
        logk = logk[pressure >= 4000]
        delG = delG[pressure >= 4000]
        pressure = pressure[pressure >= 4000]
    return temperature, pressure, logk, delG

def calculate_logq(activities, reaction):
    logq = 0
    for species, coefficient in reaction.items():
        if species in activities:
            print(coefficient, activities[species])
            logq += coefficient * np.log10(activities[species])
        else:
            warnings.warn(f"Species '{species}' not found in the activities. Skipping it in the calculation of logQ.")
    return logq

def read_planet_profiles():
    tit = pd.read_csv('Titan PT profiles/TitanProfile_MgSO4_100.0ppt_Tb255.0K_PorousRock.txt', 
                      skiprows=80, delim_whitespace=True)
    eur = pd.read_csv('Europa PT profiles/EuropaProfile_Seawater_35.2ppt_Tb268.305K.txt',
                      skiprows=80, delim_whitespace=True)
    enc = pd.read_csv('Enceladus PT profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt',
                      skiprows=80, delim_whitespace=True)
    gan = pd.read_csv('Ganymede PT profiles/GanymedeProfile_PureH2O_Tb258.86K.txt', 
                      skiprows=80, delim_whitespace=True)
    earth = pd.read_csv('Earth PT profiles/Lost_city_PT_profile_340T_U1309D.csv')
    conv = 1e-2 # from MPa to kb -- 1 megapascals = 0.01 kilobars
    titP, titT = tit["P(MPa)"]*conv, tit["T(K)"]
    eurP, eurT = eur["P(MPa)"]*conv, eur["T(K)"]
    encP, encT = enc["P(MPa)"]*conv, enc["T(K)"]
    ganP, ganT = gan["P(MPa)"]*conv, gan["T(K)"]
    earthP, earthT = earth["pressure_mean_bar"]*1e-3, earth["Temperature_C"]+273.25
    return titP, titT, eurP, eurT, encP, encT, ganP, ganT, earthP, earthT

def calculate_affinity(logk, logq, T):
    affinity = 2.303 * 8.314 * T * (logk - logq) # J/molK * K = J/mol
    return affinity

def plot_affinity(T_dew, P_dew, T_sup, P_sup, aff_dew, delG_dew, aff_sup, delG_sup, OUT_PATH, OUT_FILENAME_TITLE):
    titP, titT, eurP, eurT, encP, encT, ganP, ganT, earthP, earthT = read_planet_profiles()
    lw, wl = 6, 1
    width = 25
    height = width/5.+2
    fig = plt.figure(figsize=(width,height)) # the entire figure space

    xmin, xmax = np.log10(230), np.log10(1273)
    ymin, ymax = 1e-2, 60

    plt.rcParams.update({'font.size': 22})

    plt.subplot(1,2,1) # the first subplot in [1,2]
    plt.ticklabel_format(style='plain')
    
    aff_dew = aff_dew * 1e-3 # J/mol -> kJ/mol
    aff_sup = aff_sup * 1e-3 # J/mol -> kJ/mol

    if np.min(aff_dew)*np.max(aff_dew) < np.min(aff_sup)*np.max(aff_sup):
        if np.min(aff_dew) > 0 and np.max(aff_dew) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=0., vcenter=np.min(aff_dew), vmax=np.max(aff_dew))
        if np.min(aff_dew) < 0 and np.max(aff_dew) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(aff_dew), vcenter=0., vmax=np.max(aff_dew))
        if np.min(aff_dew) < 0 and np.max(aff_dew) < 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(aff_dew), vcenter=np.max(aff_dew), vmax=0.)
        colbarscale = aff_dew
    else:
        if np.min(aff_sup) > 0 and np.max(aff_sup) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=0., vcenter=np.min(aff_sup), vmax=np.max(aff_sup))
        if np.min(aff_sup) < 0 and np.max(aff_sup) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(aff_sup), vcenter=0., vmax=np.max(aff_sup))
        if np.min(aff_sup) < 0 and np.max(aff_sup) < 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(aff_sup), vcenter=np.max(aff_sup), vmax=0.)
        colbarscale = aff_sup
        
    cp_dew = plt.tricontourf(np.log10(T_dew), P_dew*1e-3, aff_dew, 20, cmap='Spectral', norm=divnorm)
    plt.tricontour(np.log10(T_dew), P_dew*1e-3, aff_dew, levels=[0.], linestyles='dashed', colors='black')
    cp_sup = plt.tricontourf(np.log10(T_sup), P_sup*1e-3, aff_sup, 20, cmap='Spectral', norm=divnorm)
    plt.tricontour(np.log10(T_sup), P_sup*1e-3, aff_sup, levels=[0.], linestyles='dashed', colors='black')

    plt.autoscale(False)
    plt.plot(np.log10(eurT), eurP, label='Europa', linewidth=lw)
    plt.plot(np.log10(eurT), eurP, label='Europa', linewidth=wl, color='white')
    plt.plot(np.log10(ganT), ganP, label='Ganymede', linewidth=lw)
    plt.plot(np.log10(ganT), ganP, label='Ganymede', linewidth=wl, color='white')
    plt.plot(np.log10(titT), titP, label='Titan', linewidth=lw)
    plt.plot(np.log10(titT), titP, label='Titan', linewidth=wl, color='white')
    plt.plot(np.log10(earthT), earthP, label='Lost City', linewidth=lw)
    plt.plot(np.log10(earthT), earthP, label='Lost City', linewidth=wl, color='white')
    plt.plot(np.log10(encT[1:]), encP[1:], label='Enceladus', linewidth=lw)
    plt.plot(np.log10(encT[1:]), encP[1:], label='Enceladus', linewidth=wl, color='white')
    tloc, xoff, yoff = -150, [-140,30,30,-50], [0,6,2,0.25]
    plt.text(np.log10(ganT[tloc:tloc+1]+xoff[0]),ganP[tloc:tloc+1]-yoff[0], 'Ganymede', size='small')
    plt.text(np.log10(titT[tloc:tloc+1]+xoff[1]),titP[tloc:tloc+1]-yoff[1], 'Titan', size='small')
    plt.text(np.log10(eurT[tloc:tloc+1]+xoff[2]),eurP[tloc:tloc+1]-yoff[2], 'Europa', size='small')
    plt.text(2.6,7e-1, 'Lost City', size='small')
    plt.text(2.5,1e-1, 'Enceladus', size='small')

    clb1 = plt.colorbar(cp_dew)
    clb2 = plt.colorbar(cp_sup, label='Affinity (kJ mol$^{-1}$)')

    plt.xlabel('log$_{10}$(T (K))')
    plt.ylabel('P (kbar)')
    plt.ylim([ymin, ymax])
    plt.xlim([xmin, xmax])
    plt.yscale('log')
    plt.title(OUT_FILENAME_TITLE)
    plt.tick_params(which='minor', length=6, width=1.5)
    plt.tick_params(length=10, width=1.5)

    #############################################################################################
    # Delta G plot
    plt.subplot(1,2,2) # the second subplot in [1,2]
    
    #delG3_dew = delG_dew*1e-3
    delG_dew = delG_dew * 4.184 * 1e-3 # cal to J, then to kJ
    delG_sup = delG_sup * 4.184 * 1e-3 # cal to J, then to kJ

    if np.min(delG_dew)*np.max(delG_dew) < np.min(delG_sup)*np.max(delG_sup):
        if np.min(delG_dew) > 0 and np.max(delG_dew) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=0., vcenter=np.min(delG_dew), vmax=np.max(delG_dew))
        if np.min(delG_dew) < 0 and np.max(delG_dew) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(delG_dew), vcenter=0., vmax=np.max(delG_dew))
        if np.min(delG_dew) < 0 and np.max(delG_dew) < 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(delG_dew), vcenter=np.max(delG_dew), vmax=0.)
    else:
        if np.min(delG_sup) > 0 and np.max(delG_sup) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=0., vcenter=np.min(delG_sup), vmax=np.max(delG_sup))
        if np.min(delG_sup) < 0 and np.max(delG_sup) > 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(delG_sup), vcenter=0., vmax=np.max(delG_sup))
        if np.min(delG_sup) < 0 and np.max(delG_sup) < 0:
            divnorm=colors.TwoSlopeNorm(vmin=np.min(delG_sup), vcenter=np.max(delG_sup), vmax=0.)

    cp_dew = plt.tricontourf(np.log10(T_dew), P_dew*1e-3, delG_dew, 20, cmap='Spectral', norm=divnorm)
    cp_sup = plt.tricontourf(np.log10(T_sup), P_sup*1e-3, delG_sup, 20, cmap='Spectral', norm=divnorm)
    plt.tricontour(np.log10(T_dew), P_dew*1e-3, delG_dew, levels=[0.], linestyles='dashed', colors='black')
    plt.tricontour(np.log10(T_sup), P_sup*1e-3, delG_sup, levels=[0.], linestyles='dashed', colors='black')

    plt.autoscale(False)
    plt.plot(np.log10(eurT), eurP, label='Europa', linewidth=lw)
    plt.plot(np.log10(eurT), eurP, label='Europa', linewidth=wl, color='white')
    plt.plot(np.log10(ganT), ganP, label='Ganymede', linewidth=lw)
    plt.plot(np.log10(ganT), ganP, label='Ganymede', linewidth=wl, color='white')
    plt.plot(np.log10(titT), titP, label='Titan', linewidth=lw)
    plt.plot(np.log10(titT), titP, label='Titan', linewidth=wl, color='white')
    plt.plot(np.log10(earthT), earthP, label='Lost City', linewidth=lw)
    plt.plot(np.log10(earthT), earthP, label='Lost City', linewidth=wl, color='white')
    plt.plot(np.log10(encT[1:]), encP[1:], label='Enceladus', linewidth=lw)
    plt.plot(np.log10(encT[1:]), encP[1:], label='Enceladus', linewidth=wl, color='white')

    tloc, xoff, yoff = -150, [-140,30,30,-50], [0,6,2,0.25]
    plt.text(np.log10(ganT[tloc:tloc+1]+xoff[0]),ganP[tloc:tloc+1]-yoff[0], 'Ganymede', size='small')
    plt.text(np.log10(titT[tloc:tloc+1]+xoff[1]),titP[tloc:tloc+1]-yoff[1], 'Titan', size='small')
    plt.text(np.log10(eurT[tloc:tloc+1]+xoff[2]),eurP[tloc:tloc+1]-yoff[2], 'Europa', size='small')
    plt.text(2.6,7e-1, 'Lost City', size='small')
    plt.text(2.5,1e-1, 'Enceladus', size='small')

    cbar1 = plt.colorbar(cp_dew)
    cbar1.formatter.set_powerlimits((0, 0))
    cbar1.formatter.set_useMathText(True)
    cbar2 = plt.colorbar(cp_sup, label='$\Delta G$ (kJ)')
    cbar2.formatter.set_powerlimits((0, 0))
    cbar2.formatter.set_useMathText(True)
    plt.xlabel('log$_{10}$(T (K))')
    plt.ylabel('P (kbar)')
    plt.ylim([ymin, ymax])
    plt.xlim([xmin, xmax])
    plt.yscale('log')

    #plt.title(''+r'$\alpha $'+'-ketoglutarate to Succinate')
    plt.title(OUT_FILENAME_TITLE)
    plt.tick_params(which='minor', length=6, width=1.5)
    plt.tick_params(length=10, width=1.5)

    plt.savefig(OUT_PATH+OUT_FILENAME_TITLE+'.pdf', format='pdf')
#    plt.show()

# Example usage
activities_file = 'activities-enceladus-new-ph11.txt' # EI: Activities here are from Canovas&Shock20, where they 
# take values from experimental conditions in E. Coli bacteria, etc. I only modified pH from 7 to 9,
# to represent Enceladus ocean, thought to have 9 < pH < 11. 
out_path_plots = 'out/Enceladus-ph11/Prebiotic/'
out_path_csvs = out_path_plots
suffix = ''
# The names below will be used in plot titles and filenames of PDF output plots:
titles = [ 'Glucose to Pyruvate', 'Pyruvate to Alanine', 'Pyruvate to Lactate', 'Pyruvate to Oxaloacetate', 
         'Pyruvate to Acetate', 'Acetate to Citrate']
fnames = [ 'Glucose_Pyruvate', 'Pyruvate_Alanine', 'Pyruvate_Lactate', 'Pyruvate_Oxaloacetate', 
         'Pyruvate_Acetate', 'Acetate_Citrate']
#name = [ 'Citrate formation', 'Cis-aconitate formation ', 'Isocitrate formation', 'α-ketoglutarate formation', 'Succinate formation', 'Fumarate formation',
#         'Malate formation', 'Oxaloacetate formation' ]
# titles = [ 'Citrate to Cis-aconitate', 'Cis-aconitate to Isocitrate', 'Isocitrate to α-ketoglutarate', 'α-ketoglutarate to Succinate', 
#         'Succinate to Fumarate', 'Fumarate to Malate', 'Malate to Oxaloacetate', 'Oxaloacetate to Citrate']
# fnames = [ 'Citrate_Cis-aconitate', 'Cis-aconitate_Isocitrate', 'Isocitrate_a-ketoglutarate', 'a-ketoglutarate_Succinate', 
#         'Succinate_Fumarate', 'Fumarate_Malate', 'Malate_Oxaloacetate', 'Oxaloacetate_Citrate']

# Here, you set up the reaction stoichiometry dictionary. Use one dictionary 'chapter' enclosed 
# by {...}, following the format below. The labels of all species have to be found in the 
# file 'activities.txt'. If you don't find a species in that file, you should add it in with its 
# estimated activity, consistent with the format in activity.txt, on a new line. 
# Note that the reactants have negative stoi. number, the products positive below. 
#reaction = {'CH4': 1, 'H2O': 2, 'CO2': -1, 'H2': -4} # --> Single reaction; methanogenesis (test case)

# reactions = [{'citrate-3': 1, 'H2O': 5, 'H+': 3, 'CO2': -6, 'H2': -9},
#            {'cis-aconitate-3': 1, 'H2O': 6, 'H+': 3, 'CO2': -6, 'H2': -9},
#            {'isocitrate-3': 1, 'H2O': 5, 'H+': 3, 'CO2': -6, 'H2': -9},
#            {'a-ketoglutarate-2': 1, 'H2O': 5, 'H+': 2, 'CO2': -5, 'H2': -8},
#            {'succinate-2': 1, 'H2O': 4, 'H+': 2, 'CO2': -4, 'H2': -7},
#            {'fumarate-2': 1, 'H2O': 4, 'H+': 2, 'CO2': -4, 'H2': -6},
#            {'malate-2': 1, 'H2O': 3, 'H+': 2, 'CO2': -4, 'H2': -6},
#            {'oxaloacetate-2': 1, 'H2O': 3, 'H+': 2, 'CO2': -4, 'H2': -5}]

#reactions = [{'acetate-': -1,  'oxaloacetate-2': -1, 'ATP-4': -1, 'H2O': -1, 'citrate-3': 1, 'ADP-3': 1,'PO4-3': 1}]
#reactions = [{'pyruvate-': -1, 'ADP-3': -1,'PO4-3': -1, 'acetate-': 1, 'ATP-4': 1, 'CO2' : 1}]

reactions = [{'glucose': -1, 'ADP-3': -2,'PO4-3': -2, 'pyruvate-': 2, 'H+': 2, 'ATP-4': 2, 'H2O': 2},
            {'pyruvate-': -1, 'NADH' :-1, 'H+': -1, 'NH3': -1, 'alanine': 1, 'NAD+' :1, 'OH-': 1},
            {'pyruvate-': -1, 'NADH' :-1, 'H+': -1, 'lactate-': 1, 'NAD+': 1},
            {'pyruvate-': -1, 'HCO3-': -1, 'ATP-4': -1, 'oxaloacetate-2': 1, 'ADP-3': 1, 'PO4-3': 1},
            {'pyruvate-': -1, 'ADP-3': -1,'PO4-3': -1, 'acetate-': 1, 'ATP-4': 1, 'CO2' : 1},
            {'acetate-': -1,  'oxaloacetate-2': -1, 'ATP-4': -1, 'H2O': -1, 'citrate-3': 1, 'ADP-3': 1,'PO4-3': 1}]                   
# reactions = [{'cis-aconitate-3': 1, 'H2O': 1, 'citrate-3': -1},
#              {'isocitrate-3': 1, 'cis-aconitate-3': -1, 'H2O': -1},
#              {'a-ketoglutarate-2': 1, 'NAD-2red': 1, 'CO2': 1, 'H+': 1, 'isocitrate-3': -1, 'NAD-ox': -1},
#              {'succinate-2': 1, 'NAD-2red':1,'CO2':1, 'ATP-4':1,'H+':1, 'a-ketoglutarate-2': -1, 'NAD-ox': -1, 'ADP-3': -1, 'HPO4-2': -1},
#              {'fumarate-2': 1, 'H2': 1, 'succinate-2': -1},
#              {'malate-2': 1,'fumarate-2': -1, 'H2O': -1},
#              {'oxaloacetate-2': 1, 'NAD-2red':1, 'H+': 2, 'malate-2': -1, 'NAD-ox': -1},
#              {'citrate-3': 1, 'NAD-2red': 1, 'CO2':1, 'H+': 2,'pyruvate-': -1,'NAD-ox': -1,'oxaloacetate-2': -1, 'H2O': -1}]
             
# Here, you should define where the csv source files of thermodynamic data (logK, etc.) are located 
# for DEW and SUPCRT. 
# csv_dew = sorted(glob.glob(os.path.join('LogK_csvs/TCA/dew/', "*.csv")))
# csv_sup = sorted(glob.glob(os.path.join('LogK_csvs/TCA/sup/', "*.csv")))
csv_dew = sorted(glob.glob(os.path.join('LogK_csvs/Prebiotic/dew/', "*.csv")))
csv_sup = sorted(glob.glob(os.path.join('LogK_csvs/Prebiotic/sup/', "*.csv")))

# csv_sup = ['organic_activity_estimates/acetform.csv']
# csv_dew = ['organic_activity_estimates/acetateform.csv']
# csv_sup = ['LogK_csvs/Prebiotic/dew/5-PyruvatetoAcetate-DEW.csv']
# csv_dew = ['LogK_csvs/Prebiotic/sup/5-PyruvatetoAcetate-SPCRT.csv']

#print(len(reactions))
print(csv_dew)
# The main loop over reactions: 
# Read activities, read logKs from DEW and SUPCRT, determine logQ (reac.quotient based on activities),
# calculate affinities, plot affinities along with DeltaGs. 
header_dew = "T_dew, P_dew, aff_dew, delG_dew"
header_sup = "T_sup, P_sup, aff_sup, delG_sup"
for index in range(len(reactions)):
    activities = read_activities(activities_file, reactions[index])
    print(csv_dew[index], csv_sup[index])
    T_dew, P_dew, logk_dew, delG_dew = read_logk(csv_dew[index])
    T_sup, P_sup, logk_sup, delG_sup = read_logk(csv_sup[index])
    logq = calculate_logq(activities, reactions[index])
    aff_dew = calculate_affinity(logk_dew, logq, T_dew)
    aff_sup = calculate_affinity(logk_sup, logq, T_sup)
    plot_affinity(T_dew, P_dew, T_sup, P_sup, aff_dew, delG_dew, aff_sup, delG_sup, out_path_plots, titles[index])
    np.savetxt(out_path_csvs + 'csv_dew/' + fnames[index]+suffix+'.csv', np.column_stack((T_dew, P_dew, aff_dew, delG_dew)), 
               delimiter=',', header=header_dew, comments='')
    np.savetxt(out_path_csvs + 'csv_sup/' + fnames[index]+suffix+'.csv', np.column_stack((T_sup, P_sup, aff_sup, delG_sup)), 
               delimiter=',', header=header_sup, comments='')

#breakpoint()
#print(f"Chemical Affinity: {affinity:.2f} J/mol")
