"""
    net_affinity_profiles -- Authors: Seda & Emre Isik, 2024
    Calculations of the depth profile of net affinity of reactions around an input species in 
    a reaction network. 'Depth profile' means that we read the PT profiles from PlanetProfile 
    model files, filter aqueous domains through those profiles to yield aqdataframe(...), 
    and return the net affinity profile of a given species only through the depths where water 
    can be in liquid form. See Isik, S. et al. 2025 ACS Earth & Space Chemistry (in press) for details.
    REQUIREMENTS: seafreeze package (Journaux et al. 2020)
    INPUTS: - the name of the chemical species
            - reaction affinities from csv files (output of affinity.py) provided in 
            species_affinity(...)
            - PT curves of planetary bodies (icy moons and Lost City borehole profile in our case), 
            in geotherms()
    OUTPUTS: the net affinity of reactions that share the input species, as a function of fractional 
            radius -- also indicated is the name of the background model used to calculate ∆G. 
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import glob
from scipy import interpolate
from seafreeze import seafreeze as sf
def species_affinity(species, model, TCA_only=False):
    print(f'Selected species & model: {species} - {model}')  # Press ⌘F8 to toggle the breakpoint.

    if model == 'SUPCRT':
        if (TCA_only == False):
            path = 'affinity-csvs/network/Enceladus-only-ph11/CSV_SUPCRT/' #'out/Enceladus-rev2/TCA/csv_sup/'
        else: path = 'affinity-csvs/network/Enceladus-only2/CSV_SUPCRT/TCA-only/'
    else:
        if (TCA_only == False):
            path = 'affinity-csvs/network/Enceladus-only2/CSV_DEW/' #'out/Enceladus-rev2/TCA/csv_dew/'
        else: path = 'affinity-csvs/network/Enceladus-only2/CSV_DEW/TCA-only/'
    csv_files_nom = glob.glob(os.path.join(path, "*" + species + ".csv"))
    csv_files_denom = glob.glob(os.path.join(path, species + "*.csv"))
    print(csv_files_nom)
    print(csv_files_denom)
    df = pd.read_csv(csv_files_nom[0])
    if (model == 'DEW') & (species == 'Fumarate'): df = df[df[df.keys()[1]] >= 4000]
    if (model == 'DEW') & (species != 'Fumarate'): df = df[df[df.keys()[1]] >= 3000]
    if (model == 'DEW') & (species == 'Succinate'): df = df[df[df.keys()[1]] >= 4000]
    if (model == 'DEW') & (species != 'Succinate'): df = df[df[df.keys()[1]] >= 3000]

    affinity_nom = np.zeros((len(df), len(csv_files_nom)))
    affinity_denom = np.zeros((len(df), len(csv_files_denom)))
    i = 0
    for f in csv_files_nom:  # reactions in the nominator (to be summed in log)
        df = pd.read_csv(f)
        if (model == 'DEW') & (species == 'Fumarate'): df = df[df[df.keys()[1]] >= 4000]
        if (model == 'DEW') & (species != 'Fumarate'): df = df[df[df.keys()[1]] >= 3000]
        if (model == 'DEW') & (species == 'Succinate'): df = df[df[df.keys()[1]] >= 4000]
        if (model == 'DEW') & (species != 'Succinate'): df = df[df[df.keys()[1]] >= 3000]
        print(df.shape)
        print('Nominator:', i, f)
        affinity_nom[:, i] = df[df.keys()[2]]
        i = i + 1
    i = 0
    for f in csv_files_denom:  # reactions in the denominator (to be subtracted in log)
        df = pd.read_csv(f)
        if (model == 'DEW') & (species == 'Fumarate'): df = df[df[df.keys()[1]] >= 4000]
        if (model == 'DEW') & (species != 'Fumarate'): df = df[df[df.keys()[1]] >= 3000]
        if (model == 'DEW') & (species == 'Succinate'): df = df[df[df.keys()[1]] >= 4000]
        if (model == 'DEW') & (species != 'Succinate'): df = df[df[df.keys()[1]] >= 3000]
        print('Denominator:', i, f)
        affinity_denom[:, i] = df[df.keys()[2]]
        temperature, pressure = df[df.keys()[0]], df[df.keys()[1]]
        i = i + 1
    # logKnet = log_k_nom - log_k_denom
    # Below, we first sum the branching-reaction equilibrium constants for both sides of the target species, 
    # and then calculate the log difference between the nominator and denominator of both sides, 
    # ie, by 'both sides' we mean we consider the target species as a product on the left side 
    # (nominator) and as a reactant on the right side (denominator). 

    #log_k_net = np.sum(log_k_nom, axis=1) - np.sum(log_k_denom, axis=1)
    aff_net = np.sum(affinity_nom, axis=1) - np.sum(affinity_denom, axis=1)
    return aff_net, pressure, temperature

def geotherms():
    tit = pd.read_csv('Titan PT profiles/TitanProfile_MgSO4_100.0ppt_Tb255.0K_PorousRock.txt',
                      skiprows=80, delim_whitespace=True)
    eur = pd.read_csv('Europa PT profiles/EuropaProfile_Seawater_35.2ppt_Tb268.305K.txt',
                      skiprows=80, delim_whitespace=True)
    gan = pd.read_csv('Ganymede PT profiles/GanymedeProfile_PureH2O_Tb258.86K.txt',
                      skiprows=80, delim_whitespace=True)
    enc = pd.read_csv('Enceladus PT profiles/EnceladusProfile_Seawater_10.0ppt_Tb272.4578K_PorousRock.txt',
                      skiprows=80, delim_whitespace=True)
    los = pd.read_csv('Earth PT profiles/Lost_city_PT_profile_340T_U1309D.txt')

    conv = 1e-2  # from MPa to kb -- 1 megapascals = 0.01 kilobars
    titP, titT, titD = tit["P(MPa)"] * conv, tit["T(K)"], tit["r(m)"]
    eurP, eurT, eurD = eur["P(MPa)"] * conv, eur["T(K)"], eur["r(m)"]
    ganP, ganT, ganD = gan["P(MPa)"] * conv, gan["T(K)"], gan["r(m)"]
    encP, encT, encD = enc["P(MPa)"] * conv, enc["T(K)"], enc["r(m)"]
    losP, losT, losD = (los["pressure_mean_bar"] * 1e-3, los["Temperature_C"] + 273.25,
                        los["Depth_below_seafloor_m"])
    losP, losT, losD = np.asarray(losP[::200]), np.asarray(losT[::200]), np.asarray(losD[::200])
    #eurP.name, titP.name, ganP.name, losP.name = 'P(kbar)', 'P(kbar)', 'P(kbar)', 'P(kbar)'
    titDn = np.asarray(titD / np.max(titD))
    eurDn = np.asarray(eurD / np.max(eurD))
    ganDn = np.asarray(ganD / np.max(ganD))
    encDn = np.asarray(encD / np.max(encD))
    return eurP, eurT, eurDn, titP, titT, titDn, ganP, ganT, ganDn, encP, encT, encDn, losP, losT, losD

def interpol(T, P, aff, idx, extra=1):
    # Re-index P,T,affinity arrays to show coordinates of 2D mesh-indices created above.
    Tx, Px, affx = T[idx], P[idx], aff[idx]
    Td, Pd = Tx[:, 0], Px[0, :]  # --> to be used in the interpolator
    if extra == 1:
        polate = interpolate.RegularGridInterpolator((Td, Pd), affx, method='slinear',
                                                     bounds_error=False, fill_value=None)
    elif extra == 0:
        polate = interpolate.RegularGridInterpolator((Td, Pd), affx, method='slinear',
                                                     bounds_error=False)
    if (planet == 'Titan'): prof = polate((titT, titP))
    if (planet == 'Europa'): prof = polate((eurT, eurP))
    if (planet == 'Ganymede'): prof = polate((ganT, ganP))
    if (planet == 'Enceladus'): prof = polate((encT, encP))
    if (planet == 'LostCity'): prof = polate((losT, losP))
    return prof

def water_phase():
    if (planet == 'Europa'):
        PTeur = np.empty((len(eurP),), dtype='object')
        for i in range(len(eurP)):
            PTeur[i] = (eurP[i] * 1e2, eurT[i])  # P (MPa), T (K)
        out_eur = sf.whichphase(PTeur)
        out_eur[np.isnan(out_eur)] = 7.
        out = out_eur
    if (planet == 'Titan'):
        PTtit = np.empty((len(titP),), dtype='object')
        for i in range(len(titP)):
            PTtit[i] = (titP[i] * 1e2, titT[i])  # P (MPa), T (K)
        out_tit = sf.whichphase(PTtit)
        out_tit[np.isnan(out_tit)] = 7.
        out = out_tit
    if (planet == 'Ganymede'):
        PTgan = np.empty((len(ganP),), dtype='object')
        for i in range(len(ganP)):
            PTgan[i] = (ganP[i] * 1e2, ganT[i])  # P (MPa), T (K)
        out_gan = sf.whichphase(PTgan)
        out_gan[np.isnan(out_gan)] = 7.
        out = out_gan
    if (planet == 'Enceladus'):
        PTenc = np.empty((len(encP),), dtype='object')
        for i in range(len(encP)):
            PTenc[i] = (encP[i] * 1e2, encT[i])  # P (MPa), T (K)
        out_enc = sf.whichphase(PTenc)
        out_enc[np.isnan(out_enc)] = 7.
        out = out_enc
    if (planet == 'LostCity'):
        PTlos = np.empty((len(losP),), dtype='object')
        for i in range(len(losP)):
            PTlos[i] = (losP[i] * 1e2, losT[i])  # P (MPa), T (K)
        out_los = sf.whichphase(PTlos)
        out_los[np.isnan(out_los)] = 7.
        out = out_los
    return out

def get_profiles(planet, profi):
    if (planet == 'Europa'):
        pro_ = profi
        out_ = phases
        plaDn = eurDn
        plaP = eurP
        plaT = eurT
    if (planet == 'Titan'):
        pro_ = profi
        out_ = phases
        plaDn = titDn
        plaP = titP
        plaT = titT
    if (planet == 'Ganymede'):
        pro_ = profi
        out_ = phases
        plaDn = ganDn
        plaP = ganP
        plaT = ganT
    if (planet == 'Enceladus'):
        pro_ = profi
        out_ = phases
        plaDn = encDn
        plaP = encP
        plaT = encT
    if (planet == 'LostCity'):
        pro_ = profi
        out_ = phases
        plaDn = losD
        plaP = losP
        plaT = losT
    return pro_, out_, plaDn, plaP, plaT

def aqdataframe(plaDn, plaP, plaT, phase, pro, model):
    # Send T, P, logKnet to interpolator, to return the inter/extra-polated logKnets.
    # Works for both models, sequentially.

    data = {'R': plaDn, 'P': plaP, 'Te': plaT, 'phase': phase, 'aff_net': pro,
            'model': model}

    df = pd.DataFrame(data)
    df_aq = df[df.phase == 0]
    return df_aq

if __name__ == '__main__':
    planet = 'Enceladus' # planet to plot
    TCA_only = False # If you calculate profiles exclusively based on the TCA reactions.
                    # Default is False, meaning we consider the entire network.
    species = input('Please enter the central species, e.g., Citrate\n'
                        'or type "exit" to quit:> ')
#        tmp = input('Press ENTER to continue...')
    print(f'Processing species ****{species}****')
    aff_net, P, T = species_affinity(species, 'SUPCRT', TCA_only)
    print('SUPCRT | min(aff_net)= %.2f, max(aff_net)= %.2f' % (np.min(aff_net), np.max(aff_net)))
    # Set up a dataframe for Temperature,Pressure,aff_net
    # All available from df: Temperature,Pressure,delV,delG,aff
    net_data = {"Temperature": T,
                "Pressure": P,
                "Affinity_net": aff_net}
    dfnet_sup = pd.DataFrame(net_data)
    # If method = SUPCRT, then limit temperature to 275 - 400 K and eliminate 1-bar lines.
    dfnet_sup = dfnet_sup[(dfnet_sup.Temperature > 275) & (dfnet_sup.Temperature <= 400.15)
                          & (dfnet_sup.Pressure > 1)]

    aff_net, P, T = species_affinity(species, 'DEW', TCA_only)
    print('DEW | min(aff_net)= %.2f, max(aff_net)= %.2f' % (np.min(aff_net), np.max(aff_net)))
    net_data = {"Temperature": T,
                "Pressure": P,
                "Affinity_net": aff_net}
    dfnet_dew = pd.DataFrame(net_data)

    # T, P, aff_net of SUPCRT model (final units: K, kbar, dimless)
    T_sup, P_sup, aff_net_sup = (np.array(dfnet_sup.Temperature), np.array(dfnet_sup.Pressure) / 1e3,
                                 np.array(dfnet_sup.Affinity_net))
    # T, P, aff_net of DEW model (final units: K, kbar, dimless)
    T_dew, P_dew, aff_net_dew = (np.array(dfnet_dew.Temperature), np.array(dfnet_dew.Pressure) / 1e3,
                                 np.array(dfnet_dew.Affinity_net))

    # Read geotherms in...
    eurP, eurT, eurDn, titP, titT, titDn, ganP, ganT, ganDn, encP, encT, encDn, losP, losT, losD = geotherms()

    # Define indices of 2D matrices that conform P,T arrays of sup & dew.
    # This step is independent of the planet chosen - only reshapes the PT ranges.
    idx_sup = np.lexsort((P_sup,T_sup)).reshape(126,18)
    if (species == 'Fumarate') or (species == 'Succinate'): 
        idx_dew = np.lexsort((P_dew,T_dew)).reshape(91,57) # --> For fumarate only!!!! EI 21.10.2024
    else:
        idx_dew = np.lexsort((P_dew,T_dew)).reshape(91,58) # --> For normal use! (58)

    # Determine the water/ice phases along the geotherms: out_gan, etc. defined here.
    phases = water_phase()

    # Determine LogKnet profiles for all the planets:
    profile = interpol(T_sup, P_sup, aff_net_sup, idx_sup, 1)
    # Get profiles of D,P,T,phase,LogKnet for the given planet, picking SUP-logK profile eg, proeur
    # from the line above.
    pro, phase, plaDn, plaP, plaT = get_profiles(planet, profile)
    # Filter the aqueous state only.
    #if (planet != 'LostCity'):
    df_aq_sup = aqdataframe(plaDn, plaP, plaT, phase, pro, model='sup')
    #else:
    #    df_aq_sup = aqdataframe(plaDn[::10], plaP[::10], plaT[::10], phase[::10], pro[::10], model='sup')

    # Redo the above block for DEW --> could be functionised later.
    profile = interpol(T_dew, P_dew, aff_net_dew, idx_dew, 1)
    pro, phase, plaDn, plaP, plaT = get_profiles(planet, profile)
    if (planet != 'LostCity' and planet != 'Enceladus'):
        df_aq_dew = aqdataframe(plaDn, plaP, plaT, phase, pro, model='dew')
    # Compile sup and dew results in one dataframe.
        df_aq_all = pd.concat([df_aq_sup, df_aq_dew], ignore_index=True)
    else:
        df_aq_all = df_aq_sup.copy()

    # Filter valid ranges of R and aff_net, to obtain two concatenated dataframes (one for each model).
    # By 'valid range', I mean the range that is covered by the PT range of each model.
    # So, let's first determine the boundaries of P & T in each model.
    #print('min/max(P_dew)=', np.min(P_dew), np.max(P_dew))
    # Set up the model-independent 1D interpolation tables between the pressure and the radial location,
    # ie, from the geotherms:
    RintP = interpolate.interp1d(plaP, plaDn)
    RintT = interpolate.interp1d(plaT, plaDn)
    # Determine DEW's upper boundary twice, by interpolating P and T separately:
    if (planet != 'LostCity' and planet != 'Enceladus'):
        RmaxP_dew = RintP([np.min(P_dew)])
        RmaxT_dew = RintT([np.min(T_dew)])
    # Take the mean of both interpolates, to represent the upper boundary of DEW.
        Rmax_mean_dew = np.mean([RmaxP_dew, RmaxT_dew])
        print('DEW RmaxP, RmaxT, <Rmax>=', RmaxP_dew, RmaxT_dew, Rmax_mean_dew)

    # Print out minmax values of SUPCRT PTs.
    print('min/max(P_sup)=', np.min(P_sup), np.max(P_sup))
    print('min/max(T_sup)=', np.min(T_sup), np.max(T_sup))
    # Determine SUPCRT's lower boundary twice, by interpolating P and T separately:
    if (planet != 'LostCity' and planet != 'Enceladus'):
        RminP_sup = RintP([np.max(P_sup)])
    # Take the mean of both interpolates, to represent the lower boundary of SUPCRT.
        if (planet != 'Ganymede'):
            RminT_sup = RintT([np.max(T_sup)])
            Rmin_mean_sup = np.mean([RminP_sup, RminT_sup])
        else:
            Rmin_mean_sup = np.mean([RminP_sup])
    # Determine SUPCRT's upper boundary twice, by interpolating P and T separately:
        RmaxP_sup = RintP([np.min(P_sup)])
        RmaxT_sup = RintT([np.min(T_sup)])
    # Take the mean of both interpolates, to represent the upper boundary of SUPCRT.
        Rmax_mean_sup = np.mean([RmaxP_sup, RmaxT_sup])

        print('SUP Rmin=', Rmin_mean_sup)
        print('SUP Rmax=', Rmax_mean_sup)
        print('SUP RmaxP_sup, RmaxT_sup=', RmaxP_sup, RmaxT_sup)
    #print('The lowest-temperature boundary of SUPCRT range is intersected by '
    #      'Europa geotherm at RmaxT_sup.')
    valid_sup = df_aq_all[(df_aq_all.model == 'sup')] # & (df_aq_all.R >= Rmin_mean_sup)]
    valid_dew = df_aq_all[(df_aq_all.model == 'dew')] # & (df_aq_all.R <= Rmax_mean_dew)]
    if (planet != 'LostCity' and planet != 'Enceladus'):
        valid_aq = pd.concat([valid_sup, valid_dew], ignore_index=True)
    else: valid_aq = valid_sup.copy()
    # Write out affinity_net profile to a csv file.
    if (TCA_only == False):
        valid_aq.to_csv('out/'+planet+'-ph11/'+planet + '_' + species + '.csv', columns=['R', 'aff_net', 'model'])
    else:
        valid_aq.to_csv('out/'+planet+'-ph11/'+planet + '_' + species + '_TCA-only.csv', columns=['R', 'aff_net', 'model'])

    plt.scatter(valid_sup.R, valid_sup.aff_net, marker='>', c='C0', alpha=.5, s=20)
    if (planet != 'LostCity'):
        plt.scatter(valid_dew.R, valid_dew.aff_net, marker = '<', c='C1', alpha=.5, s=20)
        plt.xlabel('$R/R_P$')
    if (planet == 'LostCity'): plt.xlabel('Depth (m)')
    if (planet == 'Enceladus'): plt.xlabel('$R/R_P$')
    plt.title(planet)
    plt.ylabel('$A_{net}$')
    plt.show()

