import meep as mp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import PyMieScatt as ps
import sys

um_scale = 0.1 # 100nm
eV_um_scale = um_scale/1.23984193

metal_range = mp.FreqRange(min=um_scale/6.1992, max=um_scale/.24797)

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

Au = mp.Medium(epsilon=1.0, E_susceptibilities=Au_susc, valid_freq_range=metal_range)  
Au.from_um_factor = um_scale     

r = float(sys.argv[1]) # radius of sphere                                                                                                                                    

wvl_min = 2.5
wvl_max = 12

frq_min = 1/wvl_max
frq_max = 1/wvl_min
frq_cen = 0.5*(frq_min+frq_max) 
dfrq = frq_max-frq_min
nfrq = 100

## needs to be large enough to resolve skin depth of Au sphere                                                                                         
resolution = 100

dpml = 0.5*wvl_max
dair = 0.5*wvl_max

pml_layers = [mp.PML(thickness=dpml)]

s = 2*(dpml+dair+r)
cell_size = mp.Vector3(s,s,0)

# is_integrated=True necessary for any planewave source extending into PML                                                                                     
sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     center=mp.Vector3(-0.5*s+dpml),
                     size=mp.Vector3(0,s,0),
                     component=mp.Ez)]


## symmetries sometimes cause field instabilities at certain resolutions                                                                                       
symmetries = [mp.Mirror(mp.Y)]
              #,mp.Mirror(mp.Z,phase=-1)]

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    symmetries=symmetries)

dpad = 0.5*dair
box_x1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=-r-dpad),size=mp.Vector3(0,2*(r+dpad),0)))
box_x2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=+r+dpad),size=mp.Vector3(0,2*(r+dpad),0)))
box_y1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=-r-dpad),size=mp.Vector3(2*(r+dpad),0,0)))
box_y2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=+r+dpad),size=mp.Vector3(2*(r+dpad),0,0)))

sim.run(until_after_sources=10)

freqs = mp.get_flux_freqs(box_x1)
box_x1_data = sim.get_flux_data(box_x1)
box_x2_data = sim.get_flux_data(box_x2)
box_y1_data = sim.get_flux_data(box_y1)
box_y2_data = sim.get_flux_data(box_y2)

box_x1_flux0 = mp.get_fluxes(box_x1)

sim.reset_meep()

n_sphere = 2.0
geometry = [mp.Sphere(material=Au, #mp.Medium(index=n_sphere),
                      center=mp.Vector3(),
                      radius=r)]

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    symmetries=symmetries,
                    geometry=geometry)

box_x1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=-r-dpad),size=mp.Vector3(0,2*(r+dpad),0)))
box_x2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=+r+dpad),size=mp.Vector3(0,2*(r+dpad),0)))
box_y1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=-r-dpad),size=mp.Vector3(2*(r+dpad),0,0)))
box_y2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=+r+dpad),size=mp.Vector3(2*(r+dpad),0,0)))

# for normal run, load negated fields to subtract incident from refl. fields
sim.load_minus_flux_data(box_x1, box_x1_data)
sim.load_minus_flux_data(box_x2, box_x2_data)
sim.load_minus_flux_data(box_y1, box_y1_data)
sim.load_minus_flux_data(box_y2, box_y2_data)


animate = mp.Animate2D(sim,
                       fields=mp.Ez,
                       output_plane=mp.Volume(center=mp.Vector3(0,0),size=mp.Vector3(2*(r+0.1*dpad),2*(r+0.1*dpad))),
                       field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'none'},
                       eps_parameters={'interpolation':'spline36', 'cmap':'binary','alpha':0.5, 'contour':True, 'contour_linewidth':1, 'frequency':None, 'resolution':None})

sim.run(#mp.to_appended("ez-slice", mp.output_efield_z), # need to get working
        mp.at_every(1,animate),
        until_after_sources=100)

box_x1_flux = mp.get_fluxes(box_x1)
box_x2_flux = mp.get_fluxes(box_x2)
box_y1_flux = mp.get_fluxes(box_y1)
box_y2_flux = mp.get_fluxes(box_y2)

scatt_flux = np.asarray(box_x1_flux)-np.asarray(box_x2_flux)+np.asarray(box_y1_flux)-np.asarray(box_y2_flux)
intensity = np.asarray(box_x1_flux0)#/2*(r+dpad) #/(2*r)**2
#scatt_cross_section = np.divide(scatt_flux,intensity)
scatt_eff_meep = np.divide(scatt_flux,intensity)
scatt_eff_meep *= -1
#scatt_eff_meep = scatt_cross_section*-1(np.pi*r**2)

trace = lambda t: (t[0][0]+t[1][1]+t[2][2])/3
scatt_eff_theory = [ps.MieQ(np.sqrt(trace(Au.epsilon(f))),(1/um_scale)/f,2*r*(1/um_scale),asDict=True)['Qsca'] for f in freqs]

scatt_eff_theory = np.asarray(scatt_eff_theory)

#normalize
normalized_scatt_eff_meep = (scatt_eff_meep - scatt_eff_meep.min())/ (scatt_eff_meep.max() - scatt_eff_meep.min())
normalized_scatt_eff_theory = (scatt_eff_theory - scatt_eff_theory.min())/ (scatt_eff_theory.max() - scatt_eff_theory.min())

if mp.am_master():    
    animate.to_gif(fps=1E9, filename="./gifs/NanoparticleScattering2D_Au_R{}nm.gif".format(r*100))
    plt.figure(dpi=150)
    plt.plot(1/np.asarray(freqs),normalized_scatt_eff_meep,'bo-',label='Meep') #2*np.pi*r* 
    plt.plot(1/np.asarray(freqs),normalized_scatt_eff_theory,'ro-',label='theory')
    plt.grid(True,which="both",ls="-")
    plt.xlabel('wavelength, λ')
    plt.ylabel('normalized scattering efficiency, σ') # πr$^{2}$
    plt.legend(loc='upper left')
    plt.title('Mie Scattering of Au Sphere')
    plt.tight_layout()
    plt.savefig("./scattering_efficiencies/mie_scattering_Au_R{}nm.png".format(r*100))
    np.savez("./npz/scatt_eff_res_Au_R{}nm.npz".format(r*100),scatt_eff_meep=scatt_eff_meep)