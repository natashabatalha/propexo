from chimera import *
import astropy.units as u 
import os
import pandas as pd

import picaso.justdoit as jdi

def chimera_chem(CtoO, Met, T, P,  directory=None, xsects=None):
    """
    Routine to extract chimera chemistry without running RT

    Parameters
    ----------

    logCtoO : float 
        Log C/O -0.26=Solar  
    logMet : float 
        Log Met relative to solar (0 is solar)
    T : array 
        Temperature profile array (T in kelvin)
    P : array 
        Pressure profie array (P in bars)
    directory : str , optional    
        Directory pointing to ABSCOEF (this will read in the chemistry array directly)
    xsects : array 
        This is the output from chimeras xsecs routine or `grab_stars_xsecs` in this black box codes program. 
        Input this if you have already read in the chemistry via CHIMERA and don't want to reread. 

    Returns
    -------
    pd.Dataframe of full chemistry input
    """

    #loading interpolatedable chemistry grid as a function of C/O, Metallicity, T, and P (for many gases)
    if isinstance(xsects, type(None)):
        hf=h5py.File(os.path.join(directory,'CHEM','chem_grid.h5'), 'r')
        logCtoO=np.array(hf['logCtoO'])  #-2.0 - 0.3
        logMet=np.array(hf['logMet']) #-2.0 - 3.0
        Tarr=np.array(hf['Tarr'])  #400 - 3400
        logParr=np.array(hf['logParr'])   #-7.0 - 2.4
        gases=np.array(hf['gases'])  ##H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He   e- h-  mmw
        loggas = np.log10(gases)
        hf.close()
    else: 
        #interpolation chem
        logCtoO, logMet, Tarr, logParr, loggas=xsects[9:]  #see xsects_HST/JWST routine...


    Ngas=loggas.shape[-2]
    gas=np.zeros((Ngas,len(P)))+1E-20

    #capping T at bounds
    TT=np.zeros(len(T))
    TT[:]=T[:]
    TT[TT>3400]=3400
    TT[TT<500]=500

    for j in range(Ngas):
        gas_to_interp=loggas[:,:,:,j,:]
        IF=RegularGridInterpolator((logCtoO, logMet, np.log10(Tarr),logParr),gas_to_interp,bounds_error=False)
        for i in range(len(P)):
            gas[j,i]=10**IF(np.array([CtoO, Met, np.log10(TT[i]), np.log10(P[i])]))

    H2Oarr, CH4arr, COarr, CO2arr, NH3arr, N2arr, HCNarr, H2Sarr,PH3arr, C2H2arr, C2H6arr, Naarr, Karr, TiOarr, VOarr, FeHarr, Harr,H2arr, Hearr,earr, Hmarr,mmw=gas
    df = pd.DataFrame({'H2O': H2Oarr, 'CH4': CH4arr, 'CO': COarr, 'CO2': CO2arr, 'NH3': NH3arr, 
                           'N2' : N2arr, 'HCN': HCNarr, 'H2S': H2Sarr, 'PH3': PH3arr, 'C2H2': C2H2arr, 
                           'C2H6' :C2H6arr, 'Na' : Naarr, 'K' : Karr, 'TiO': TiOarr, 'VO' : VOarr, 
                           'Fe': FeHarr,  'H': Harr, 'H2' : H2arr, 'He' : Hearr, 'e-':earr,'H-':Hmarr, 'temperature':T, 
                           'pressure': P})
    return df

def grab_stars_xsecs(df, xsecs_dir, stellar_dir, verbose = False, **planet_kwargs):
    for i in df.index:

        planet = df.loc[i,:].dropna().to_dict()

        stellar_file = planet.get('pl_name', planet_kwargs.get('pl_name','sum_planet')) + '_stellar_file.h5'

        temp = planet.get('st_teff', planet_kwargs.get('st_teff',np.nan))
        if np.isnan(temp) : raise Exception('Stellar temperature is not added to \
            dataframe input or to planet_kwargs through the column/key named st_teff. Please add it to one of them')

        logg = planet.get('st_logg', planet_kwargs.get('st_logg',np.nan))
        if np.isnan(logg) : raise Exception('Stellar logg is not added to \
            dataframe input or to planet_kwargs through the column/key named st_logg. Please add it to one of them')

        logmh = planet.get('st_metfe', planet_kwargs.get('st_metfe',np.nan))
        if np.isnan(logmh) : raise Exception('Stellar Fe/H is not added to \
            dataframe input or to planet_kwargs through the column/key named st_metfe. Please add it to one of them')

        stellar_db = 'phoenix'

        if logmh > 0.5: 
            if verbose: print ('Stellar M/H exceeded max value of 0.5. Value has been reset to the maximum')
            logmh = 0.5
        elif logmh < -4.0 :
            if verbose: print ('Stellar M/H exceeded min value of -4.0 . Value has been reset to the mininum')
            logmh = -4.0 

        if logg > 4.5: 
            if verbose: print ('Stellar logg exceeded max value of 4.5. Value has been reset to the maximum')
            logg = 4.5

        make_stellar(temp,logmh,logg,stellar_db,os.path.join(stellar_dir, stellar_file))

        wnomin = 750
        wnomax = 15000
        observatory='JWST'
        xsecs=xsects(wnomin, wnomax, observatory, xsecs_dir,stellar_file=os.path.join(stellar_dir, stellar_file))

        return xsecs


def chimera_transit(df, xsecs, logMH, logCtoO , logKcld=-30, logRayAmp=-30, RaySlope=4,
                wlgrid=np.arange(1, 12, 0.01), exclude_mol=[], **planet_kwargs):
    """
    Wrapper to simplify Chimera run. This really turns chimera into a black box. This was created 
    specifically Sagan School tutorial. 

    Parameters
    -----------
    df : pd.DataFrame
        See Sagan School Tutorial. This is a subset of queried Nexci API dataframe. 
    xsecs : np.array 
        Cross section array loaded from `grab_stars_xsecs`
    logMH : float , optional
        Log metallicity relative to solar (0 is solar)
    logCtoO : float , optional
        Log Carbon to Oxygen ratio (-0.26 is solar)
    logKcld : float 
        simple grey cloud parameterization (log cross section) -30 is usually no cloud. Increasing to around -25 
        introduces increasingly stronger cloud 
    logRayAmp : float
        simple Rayleigh parameterization (log cross section) -30 is usually no cloud. Increasing to around -25 
        introduces increasingly stronger rayleigh strength 
    RaySlope : float  
        Rayleigh slope (4 is earth rayleight scattering slope (Ray ~ lambda^(-4)))
    wlgrid : array 
        Wavelength grid in micron to bin chimera output to         
    exclude_mol : list 
        List of molecules to switch off contribution to spectrum without altering the scale height calculation 
    planet_kwargs : dict 
        List of parameters to supply NexSci database is values don't exist 

    Returns 
    -------
    wno, F,chemarr,y_binned

        wavenumber array, unbinned spectrum, chemistry and PT array, binned spectrum 
    """
    x = load_x()

    first = True
    for i in df.index:

        planet = df.loc[i,:].dropna().to_dict()

        stellar_file = planet.get('pl_name', planet_kwargs.get('pl_name','sum_planet')) + '_stellar_file.h5'    

        #the parameters
        #planet/star system params--typically not free parameters in retrieval
        # Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)

        Rp = planet.get('pl_radj', planet_kwargs.get('pl_radj',np.nan))
        if np.isnan(Rp) : raise Exception('Planet Radii is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_radj. J for JUPITER! \
            Please add it to one of them')


        #Stellar Radius in Solar Radii
        Rstar = planet.get('st_rad', planet_kwargs.get('st_rad',np.nan))
        if np.isnan(Rstar) : raise Exception('Stellar Radii is not added to \
            dataframe input or to planet_kwargs through the column/key named st_rad. Solar radii! \
            Please add it to one of them')

        #Mass in Jupiter Masses
        M = planet.get('pl_bmassj', planet_kwargs.get('pl_bmassj',np.nan))
        if np.isnan(M) : raise Exception('Planet Mass is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_bmassj. J for JUPITER! \
            Please add it to one of them')  

        #TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
        Tirr=planet.get('pl_eqt', planet_kwargs.get('pl_eqt',np.nan))

        if np.isnan(Tirr): 
            p =  planet.get('pl_orbper', planet_kwargs.get('pl_orbper',np.nan))
            p = p * (1*u.day).to(u.yr).value #convert to year 
            a =  (p**(2/3)*u.au).to(u.R_sun).value
            temp = planet.get('st_teff', planet_kwargs.get('st_teff',np.nan))
            Tirr = temp * np.sqrt(Rstar/(2*a))

        if np.isnan(Tirr): raise Exception('Planet Eq Temp is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_eqt. Kelvin \
            Please add it to one of them')  

        #These are the cloud ones
        logMet=logMH#x[1]#1.5742E-2 #.   #Metallicity relative to solar log--solar is 0, 10x=1, 0.1x = -1 used -1.01*log10(M)+0.6
        logCtoO=logCtoO#x[2]#-1.97  #log C-to-O ratio: log solar is -0.26
        #simple 'grey+rayleigh' parameters
        logKcld = logKcld#-40
        logRayAmp = logRayAmp#-30
        RaySlope =RaySlope#0

        
        #model parameters for TP profile 
        Tint=x['Tint']
        logKir=x['logKir']  #TP profile IR opacity controlls the "vertical" location of the gradient
        logg1=x['logg1']     #single channel Vis/IR opacity. Controls the delta T between deep T and TOA T
        #Composition parameters---assumes "chemically consistnat model" described in Kreidberg et al. 2015


        #lets ignore these for now
        logPQCarbon=x['logPQCarbon']  #CH4, CO, H2O Qunech pressure--forces CH4, CO, and H2O to constant value at quench pressure value
        logPQNitrogen=x['logPQNitrogen'] #N2, NH3 Quench pressure--forces N2 and NH3 to ""  --ad hoc for chemical kinetics--reasonable assumption
        #A&M Cloud parameters
        logKzz=9 #log Rayleigh Haze Amplitude (relative to H2)
        fsed=1.0 #haze slope--4 is Rayeigh, 0 is "gray" or flat.  
        logPbase=1.5  #gray "large particle" cloud opacity (-35 - -25)
        logCldVMR=-15.0 #cloud fraction
        xRp=x['xRp']
        

        #seting up input state vector. Must be in this order as indicies are hard wired in fx inside fm
        chimera_in=np.array([Tirr, logKir,logg1, Tint, logMet, logCtoO, logPQCarbon,logPQNitrogen, Rp*xRp, Rstar, M, logKzz, fsed,logPbase,logCldVMR, logKcld, logRayAmp, RaySlope])

        #calling forward model
        #thermochemical gas profile scaling factors
        # 0   1    2    3   4    5    6     7    8    9   10    11   12   13    14   15   16   17   18  19 20   21
        #H2O  CH4  CO  CO2 NH3  N2   HCN   H2S  PH3  C2H2 C2H6  Na    K   TiO   VO   FeH  H    H2   He   e- h-  mmw
        gas_scale=np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1., 1., 1.]) #can be made free params if desired (won't affect mmw)#can be made free params if desired (won't affect mmw)
        
        mols = avail_mols()
        for im in exclude_mol:
            gas_scale[mols.index(im)] = 0 

        y_binned,F,wno,chemarr=fx_trans(chimera_in,wlgrid,gas_scale, xsecs)#,iso=iso)  #returns model spectrum, wavenumber grid, and vertical abundance profiles from chemistry

    return wno, F,chemarr,y_binned

def picaso_load_planet(df, opacity, verbose=False,  **planet_kwargs):
    """
    Wrapper to simplify PICASO run. This really turns chimera into a black box. This was created 
    specifically Sagan School tutorial. 

    Parameters
    -----------
    df : pd.DataFrame
        See Sagan School Tutorial. This is a subset of queried Nexci API dataframe. 
    opacity : np.array 
        Opacity loaded from opannection
    chimera_dir : str 
        Pointer to chimera directory so that we have consistent chemistry
    planet_kwargs : dict 
        List of parameters to supply NexSci database is values don't exist 


    Returns 
    -------
    dict 
        dictionary of picaso full output
    """
    x = load_x()

    for i in df.index:

        planet = df.loc[i,:].dropna().to_dict()

        temp = planet.get('st_teff', planet_kwargs.get('st_teff',np.nan))
        if np.isnan(temp) : raise Exception('Stellar temperature is not added to \
            dataframe input or to planet_kwargs through the column/key named st_teff. Please add it to one of them')

        logg = planet.get('st_logg', planet_kwargs.get('st_logg',np.nan))
        if np.isnan(logg) : raise Exception('Stellar logg is not added to \
            dataframe input or to planet_kwargs through the column/key named st_logg. Please add it to one of them')

        logmh = planet.get('st_metfe', planet_kwargs.get('st_metfe',np.nan))
        if np.isnan(logmh) : raise Exception('Stellar Fe/H is not added to \
            dataframe input or to planet_kwargs through the column/key named st_metfe. Please add it to one of them')

        stellar_db = 'phoenix'

        if logmh > 0.5: 
            if verbose: print ('Stellar M/H exceeded max value of 0.5. Value has been reset to the maximum')
            logmh = 0.5
        elif logmh < -4.0 :
            if verbose: print ('Stellar M/H exceeded min value of -4.0 . Value has been reset to the mininum')
            logmh = -4.0 

        if logg > 4.5: 
            if verbose: print ('Stellar logg exceeded max value of 4.5. Value has been reset to the maximum')
            logg = 4.5   


        #the parameters
        #planet/star system params--typically not free parameters in retrieval
        # Planet radius in Jupiter Radii--this will be forced to be 10 bar radius--arbitrary (scaling to this is free par)

        Rp = planet.get('pl_radj', planet_kwargs.get('pl_radj',np.nan))
        if np.isnan(Rp) : raise Exception('Planet Radii is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_radj. J for JUPITER! \
            Please add it to one of them')


        #Stellar Radius in Solar Radii
        Rstar = planet.get('st_rad', planet_kwargs.get('st_rad',np.nan))
        if np.isnan(Rstar) : raise Exception('Stellar Radii is not added to \
            dataframe input or to planet_kwargs through the column/key named st_rad. Solar radii! \
            Please add it to one of them')

        #Mass in Jupiter Masses
        Mp = planet.get('pl_bmassj', planet_kwargs.get('pl_bmassj',np.nan))
        if np.isnan(Mp) : raise Exception('Planet Mass is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_bmassj. J for JUPITER! \
            Please add it to one of them')  

        #TP profile params (3--Guillot 2010, Parmentier & Guillot 2013--see Line et al. 2013a for implementation)
        Tirr=planet.get('pl_eqt', planet_kwargs.get('pl_eqt',np.nan))

        if np.isnan(Tirr): 
            p =  planet.get('pl_orbper', planet_kwargs.get('pl_orbper',np.nan))
            p = p * (1*u.day).to(u.yr).value #convert to year 
            a =  (p**(2/3)*u.au).to(u.R_sun).value
            temp = planet.get('st_teff', planet_kwargs.get('st_teff',np.nan))
            Tirr = temp * np.sqrt(Rstar/(2*a))

        if np.isnan(Tirr): raise Exception('Planet Eq Temp is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_eqt. Kelvin \
            Please add it to one of them') 

        p=planet.get('pl_orbper', planet_kwargs.get('pl_orbper',np.nan))

        if np.isnan(Tirr): raise Exception('Orbital Period is not added to \
            dataframe input or to planet_kwargs through the column/key named pl_orbper. Days Units') 
        else: 
            p = p * (1*u.day).to(u.yr).value #convert to year 
            a =  p**(2/3) #semi major axis in AU


        #Run picaso
        start_case = jdi.inputs()
        start_case.phase_angle(0) #radians

        #define gravity
        start_case.gravity(mass=Mp, mass_unit=u.Unit('M_jup'),
                            radius=Rp, radius_unit=u.Unit('R_jup')) #any astropy units available

        #define star
        start_case.star(opacity, temp,logmh,logg,radius=Rstar, radius_unit=u.Unit('R_sun'),
                            semi_major=a, semi_major_unit=u.Unit('au'),
                            database = stellar_db ) #opacity db, pysynphot database, temp, metallicity, logg


        start_case.guillot_pt(Tirr, x['Tint'], x['logg1'],x['logKir'],x['alpha'])


        return start_case


def load_x():
    return {"CalcType":"chemeq",
    "T":None, 
    "Rp": None,
    "Rstar":None,
    "Mp":None, 
    "logKir": -1.5, 
    "Tint":100,
    "logg1":-1, 
    "logg2":-1, 
    "alpha":0.5, 
    "logKcld":-40, 
    "logRayAmp":1, 
    "HazeSlope":4,
    "xRp":1.0, 
    "logPQCarbon":-5.0, 
    "logPQNitrogen":-5.0,
    "chemeq":{
        "logMet":0, 
        "logCtoO":-0.26
    }}

def errorbar(fig, x, y, xerr=None, yerr=None,  color='blue',
             point_kwargs={}, error_kwargs={}):

  fig.circle(x, y,color=color, **point_kwargs)

  if not isinstance(xerr,type(None)):
      x_err_x = []
      x_err_y = []
      for px, py, err in zip(x, y, xerr):
          x_err_x.append((px - err, px + err))
          x_err_y.append((py, py))
      fig.multi_line(x_err_x, x_err_y,color=color,   **error_kwargs)

  if not isinstance(yerr,type(None)):
      y_err_x = []
      y_err_y = []
      for px, py, err in zip(x, y, yerr):
          y_err_x.append((px, px))
          y_err_y.append((py - err, py + err))
      fig.multi_line(y_err_x, y_err_y, color=color,  **error_kwargs)


def write_justification(planet_name, all_ps):

    modes_to_compare =list(all_ps.keys())

    test = 0 
    chose = modes_to_compare[0]
    for i in modes_to_compare: 
        p_array = np.array(all_ps[i])
        rejected = np.array(all_ps[i])[:,2]
        nrejected = len(rejected[rejected<0.05])
        if nrejected > test: 
            test = nrejected 
            chose = i 

    return ('We choose '+chose+
              ' in order to robustly detect an atmosphere on ' +planet_name+'. '+
              'Specifically, we explored '+str(len(p_array[:,0]))
              +' atmospheric scenarios evenly spaced between log M/H='+
              str(min(p_array[:,0]))+'-'+str(max(p_array[:,0]))+'xSolar'+
              ' and log $\kappa_{cld}$= '+str(np.min(p_array[:,1]))+' - '+str(np.max(p_array[:,1]))+
              '. With '+ chose+', '+str(test)+' of these atmospherics scenarios can be ' +
              'detected with p<0.05 significance. We define atmospheric detection as the ability'+
             ' to reject, at this significance, a 2-parameter linear model. ' +
             'But really we should check with a theorist to make sure we know what we are doing. ')

def avail_mols():
    molecules = ['H2O' , 'CH4' , 'CO' , 'CO2', 'NH3' , 'N2' , 
             'HCN'  , 'H2S' , 'PH3',  'C2H2', 'C2H6' , 'Na' ,  
             'K'  , 'TiO'  , 'VO' ,  'FeH' , 'H'  ,  'H2'  , 'He'  , 
             'e-' ,'h-',  'mmw']
    return molecules
