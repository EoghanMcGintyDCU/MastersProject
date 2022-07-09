import meep as mp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import PyMieScatt as ps
import sys, os
from material import Au, SiO2, Ag, Au_JC_visible

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

dpad = 0.5*dair
dbox = 0.5*d+2*r+dpad

n_sphere = 2.0
geometry =      [mp.Sphere(material=Au_JC_visible,
                        center=mp.Vector3(0,r+0.5*d,0),
                        radius=r)
                ,mp.Sphere(material=Au_JC_visible,
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

out_str = "R{radius}nm_D{distance}".format(radius=r*100,distance=d*100)
dir = "sim_plots/2_Nanoparticles_2D_Au_JC_{}".format(out_str)

if not os.path.isdir(dir):
        os.makedirs(dir)

vol = mp.Volume(center=mp.Vector3(0,0), size=mp.Vector3(2*(r+0.5*dpad),2*(0.5*d+2*r+0.5*dpad)))

def get_png(sim):
        eps_data = sim.get_array(vol=vol, component=mp.Dielectric)
        ez_data = sim.get_array(vol=vol, component=mp.Ez)
        plt.figure()
        plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
        plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.7)
        plt.axis('off')
        plt.savefig("{dir}/{ts}.png".format(dir=dir,ts=sim.meep_time()),dpi=150,bbox_inches='tight')

sim.run(mp.at_every(1, mp.after_sources(get_png)),
        until_after_sources=100)
