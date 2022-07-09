import meep as mp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import PyMieScatt as ps
import sys, os
from material import Au, SiO2, Ag

r = float(sys.argv[1]) # radius of sphere     
d = float(sys.argv[2]) # distance between spheres                                                                                                                              

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

sx = 2*(dpml+dair+r)
sy = 2*(dpml+dair+2*r+d)
cell_size = mp.Vector3(sx,sy,0)

# is_integrated=True necessary for any planewave source extending into PML                                                                                     
sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     center=mp.Vector3(-0.5*sx+dpml),
                     size=mp.Vector3(0,sy,0),
                     component=mp.Ez)]

## symmetries sometimes cause field instabilities at certain resolutions                                                                                       
symmetries = [mp.Mirror(mp.Z,phase=-1)]

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    #symmetries=symmetries,
                    k_point=mp.Vector3())

dpad = 0.5*dair
dbox = 0.5*d+2*r+dpad
box_x1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=-r-dpad),size=mp.Vector3(0,2*(dbox),0)))
box_x2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=+r+dpad),size=mp.Vector3(0,2*(dbox),0)))
box_y1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=-dbox),size=mp.Vector3(2*(r+dpad),0,0)))
box_y2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=dbox),size=mp.Vector3(2*(r+dpad),0,0)))

sim.run(until_after_sources=10)

freqs = mp.get_flux_freqs(box_x1)
box_x1_data = sim.get_flux_data(box_x1)
box_x2_data = sim.get_flux_data(box_x2)
box_y1_data = sim.get_flux_data(box_y1)
box_y2_data = sim.get_flux_data(box_y2)

box_x1_flux0 = mp.get_fluxes(box_x1)

sim.reset_meep()

n_sphere = 2.0
geometry =      [mp.Sphere(material=mp.metal,
                        center=mp.Vector3(0,r+0.5*d,0),
                        radius=r)
                ,mp.Sphere(material=mp.metal,
                        center=mp.Vector3(0,-r-0.5*d,0),
                        radius=r)]

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    #symmetries=symmetries,
                    geometry=geometry)

box_x1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=-r-dpad),size=mp.Vector3(0,2*(dbox),0)))
box_x2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(x=+r+dpad),size=mp.Vector3(0,2*(dbox),0)))
box_y1 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=-dbox),size=mp.Vector3(2*(r+dpad),0,0)))
box_y2 = sim.add_flux(frq_cen, dfrq, nfrq, mp.FluxRegion(center=mp.Vector3(y=dbox),size=mp.Vector3(2*(r+dpad),0,0)))

# for normal run, load negated fields to subtract incident from refl. fields
sim.load_minus_flux_data(box_x1, box_x1_data)
sim.load_minus_flux_data(box_x2, box_x2_data)
sim.load_minus_flux_data(box_y1, box_y1_data)
sim.load_minus_flux_data(box_y2, box_y2_data)

vol = mp.Volume(center=mp.Vector3(0,0), size=mp.Vector3(2*(r+0.5*dpad),2*(0.5*d+2*r+0.5*dpad)))

eps_data = sim.get_array(vol=vol, component=mp.Dielectric)
eps_data = eps_data.transpose()

out_str = "R{radius}nm_D{distance}".format(radius=r*100,distance=d*100)
dir = "sim_plots/2_Nanoparticles_2D_{}".format(out_str)

if not os.path.isdir(dir):
        os.makedirs(dir)

def get_png(sim):
        ez_data = sim.get_array(vol=vol, component=mp.Ez)
        plt.figure()
        plt.imshow(eps_data, interpolation='spline36', cmap='binary')
        plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.7)
        plt.axis('off')
        plt.savefig("{dir}/{ts}.png".format(dir=dir,ts=sim.meep_time()),dpi=150,bbox_inches='tight')

sim.run(mp.at_every(1, mp.after_sources(get_png)),
        until_after_sources=100)

box_x1_flux = mp.get_fluxes(box_x1)
box_x2_flux = mp.get_fluxes(box_x2)
box_y1_flux = mp.get_fluxes(box_y1)
box_y2_flux = mp.get_fluxes(box_y2)

scatt_flux = np.asarray(box_x1_flux)-np.asarray(box_x2_flux)+np.asarray(box_y1_flux)-np.asarray(box_y2_flux)
intensity = np.asarray(box_x1_flux0)
scatt_eff_meep = np.divide(scatt_flux,intensity)
scatt_eff_meep *= -1

normalized_scatt_eff_meep = (scatt_eff_meep - scatt_eff_meep.min()) / (scatt_eff_meep.max() - scatt_eff_meep.min())

if mp.am_master():    
        plt.figure(dpi=150)
        plt.plot(1/np.asarray(freqs),normalized_scatt_eff_meep,'bo-',label='Meep')
        plt.grid(True,which="both",ls="-")
        plt.xlabel('wavelength, λ')
        plt.ylabel('normalized scattering efficiency, σ') # πr$^{2}$
        plt.legend(loc='upper left')
        plt.title('Scattering of 2 Nanospheres')
        plt.tight_layout()
        plt.savefig("scattering_efficiencies/mie_scattering_dielectric_{}.png".format(out_str))
        np.savez("npz/scatt_eff_res_{}.npz".format(out_str),scatt_eff_meep=scatt_eff_meep)