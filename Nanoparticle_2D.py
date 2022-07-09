import meep as mp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import PyMieScatt as ps
import sys
from material import Au_JC_visible

um_scale = 0.1

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
geometry = [mp.Sphere(material=Au_JC_visible,
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

sim.run(#mp.to_appended("ez-slice", mp.output_efield_z),
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