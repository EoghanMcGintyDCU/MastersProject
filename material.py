import meep as mp

# silver (Ag)
um_scale = 0.1 # 100nm
eV_um_scale = um_scale/1.23984193

metal_range_Au = mp.FreqRange(min=um_scale/6.1992, max=um_scale/.24797)

Au_plasma_frq = 9.03*eV_um_scale
Au_f0 = 0.760
Au_frq0 = 1e-10
Au_gam0 = 0.053*eV_um_scale
Au_sig0 = Au_f0*Au_plasma_frq**2/Au_frq0**2
Au_f1 = 0.024
Au_frq1 = 0.415*eV_um_scale      # 2.988 μm
Au_gam1 = 0.241*eV_um_scale
Au_sig1 = Au_f1*Au_plasma_frq**2/Au_frq1**2
Au_f2 = 0.010
Au_frq2 = 0.830*eV_um_scale      # 1.494 μm
Au_gam2 = 0.345*eV_um_scale
Au_sig2 = Au_f2*Au_plasma_frq**2/Au_frq2**2
Au_f3 = 0.071
Au_frq3 = 2.969*eV_um_scale      # 0.418 μm
Au_gam3 = 0.870*eV_um_scale
Au_sig3 = Au_f3*Au_plasma_frq**2/Au_frq3**2
Au_f4 = 0.601
Au_frq4 = 4.304*eV_um_scale      # 0.288 μm
Au_gam4 = 2.494*eV_um_scale
Au_sig4 = Au_f4*Au_plasma_frq**2/Au_frq4**2
Au_f5 = 4.384
Au_frq5 = 13.32*eV_um_scale      # 0.093 μm
Au_gam5 = 2.214*eV_um_scale
Au_sig5 = Au_f5*Au_plasma_frq**2/Au_frq5**2

Au_susc = [mp.DrudeSusceptibility(frequency=Au_frq0, gamma=Au_gam0, sigma=Au_sig0),
           mp.LorentzianSusceptibility(frequency=Au_frq1, gamma=Au_gam1, sigma=Au_sig1),
           mp.LorentzianSusceptibility(frequency=Au_frq2, gamma=Au_gam2, sigma=Au_sig2),
           mp.LorentzianSusceptibility(frequency=Au_frq3, gamma=Au_gam3, sigma=Au_sig3),
           mp.LorentzianSusceptibility(frequency=Au_frq4, gamma=Au_gam4, sigma=Au_sig4),
           mp.LorentzianSusceptibility(frequency=Au_frq5, gamma=Au_gam5, sigma=Au_sig5)]

Au = mp.Medium(epsilon=1.0, E_susceptibilities=Au_susc, valid_freq_range=metal_range_Au)       
Au.from_um_factor = um_scale   

# gold (Au)
# fit to P.B. Johnson and R.W. Christy, Physical Review B, Vol. 6, pp. 4370-9 (1972)

Au_JC_visible_frq0 = 1/(0.139779231751333*um_scale)
Au_JC_visible_gam0 = 1/(1.12834046202759*um_scale)
Au_JC_visible_sig0 = 1

Au_JC_visible_frq1 = 1/(0.404064525036786*um_scale)
Au_JC_visible_gam1 = 1/(26.1269913352870*um_scale)
Au_JC_visible_sig1 = 2.07118534879440

Au_JC_visible_susc = [mp.DrudeSusceptibility(frequency=Au_JC_visible_frq0, gamma=Au_JC_visible_gam0, sigma=Au_JC_visible_sig0),
                      mp.LorentzianSusceptibility(frequency=Au_JC_visible_frq1, gamma=Au_JC_visible_gam1, sigma=Au_JC_visible_sig1)]

Au_JC_visible = mp.Medium(epsilon=6.1599, E_susceptibilities=Au_JC_visible_susc)


metal_range_Ag = mp.FreqRange(min=um_scale/12.398, max=um_scale/.24797)

# silver (Ag)
Ag_plasma_frq = 9.01*eV_um_scale
Ag_f0 = 0.845
Ag_frq0 = 1e-10
Ag_gam0 = 0.048*eV_um_scale
Ag_sig0 = Ag_f0*Ag_plasma_frq**2/Ag_frq0**2
Ag_f1 = 0.065
Ag_frq1 = 0.816*eV_um_scale      # 1.519 um
Ag_gam1 = 3.886*eV_um_scale
Ag_sig1 = Ag_f1*Ag_plasma_frq**2/Ag_frq1**2
Ag_f2 = 0.124
Ag_frq2 = 4.481*eV_um_scale      # 0.273 um
Ag_gam2 = 0.452*eV_um_scale
Ag_sig2 = Ag_f2*Ag_plasma_frq**2/Ag_frq2**2
Ag_f3 = 0.011
Ag_frq3 = 8.185*eV_um_scale      # 0.152 um
Ag_gam3 = 0.065*eV_um_scale
Ag_sig3 = Ag_f3*Ag_plasma_frq**2/Ag_frq3**2
Ag_f4 = 0.840
Ag_frq4 = 9.083*eV_um_scale      # 0.137 um
Ag_gam4 = 0.916*eV_um_scale
Ag_sig4 = Ag_f4*Ag_plasma_frq**2/Ag_frq4**2
Ag_f5 = 5.646
Ag_frq5 = 20.29*eV_um_scale      # 0.061 um
Ag_gam5 = 2.419*eV_um_scale
Ag_sig5 = Ag_f5*Ag_plasma_frq**2/Ag_frq5**2

Ag_susc = [mp.DrudeSusceptibility(frequency=Ag_frq0, gamma=Ag_gam0, sigma=Ag_sig0),
           mp.LorentzianSusceptibility(frequency=Ag_frq1, gamma=Ag_gam1, sigma=Ag_sig1),
           mp.LorentzianSusceptibility(frequency=Ag_frq2, gamma=Ag_gam2, sigma=Ag_sig2),
           mp.LorentzianSusceptibility(frequency=Ag_frq3, gamma=Ag_gam3, sigma=Ag_sig3),
           mp.LorentzianSusceptibility(frequency=Ag_frq4, gamma=Ag_gam4, sigma=Ag_sig4),
           mp.LorentzianSusceptibility(frequency=Ag_frq5, gamma=Ag_gam5, sigma=Ag_sig5)]

Ag = mp.Medium(epsilon=1.0, E_susceptibilities=Ag_susc, valid_freq_range=metal_range_Ag)
Ag.from_um_factor = um_scale   

SiO2_range = mp.FreqRange(min=um_scale/1.77, max=um_scale/0.25)

SiO2_frq1 = 1/(0.103320160833333*um_scale)
SiO2_gam1 = 1/(12.3984193000000*um_scale)
SiO2_sig1 = 1.12

SiO2_susc = [mp.LorentzianSusceptibility(frequency=SiO2_frq1, gamma=SiO2_gam1, sigma=SiO2_sig1)]

SiO2 = mp.Medium(epsilon=1.0, E_susceptibilities=SiO2_susc, valid_freq_range=SiO2_range)
SiO2.from_um_factor = um_scale   