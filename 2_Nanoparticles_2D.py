import meep as mp
import numpy as np
import os
from material import Au, Au_JC_visible, Ag, Cu
import argparse

material_map = {"Au": Au, "Au_JC_visible": Au_JC_visible, "Ag" : Ag,  "Cu": Cu, "Eps7" : mp.Medium(epsilon=7)}

parser = argparse.ArgumentParser(
    description="Perform FDTD on 2 Nanoparticles")
parser.add_argument("material", type=str, help="Material of nanoparticles", default="Au")
parser.add_argument("savefolder", type=str, help="Folder to save results", default="./")
parser.add_argument("-r", "--radius", type=float, help="Radius of nanospheres (nm)", default=50)
parser.add_argument("-g", "--gap", type=float, help="Gap size between particles (nm)", default=5)
parser.add_argument("-s", "--surround", type=float, help="Refractive index of surrounding medium", default=1)

args = parser.parse_args()

r = args.radius / 100   # radius of sphere     
g = args.gap / 100      # distance between spheres               

save_folder = args.savefolder

if mp.am_master:
        if not os.path.isdir(save_folder):
                os.makedirs(save_folder)

material = Au
if args.material in material_map:
    material = material_map[args.material]

wvl_min = 3
wvl_max = 10

frq_min = 1/wvl_max
frq_max = 1/wvl_min
frq_cen = 0.5*(frq_min+frq_max) 
dfrq = frq_max-frq_min
nfrq = 100
                                                                                  
resolution = 50

dpml = 0.5*wvl_max
dair = 0.5*wvl_max

pml_layers = [mp.PML(thickness=dpml)]

sx = 2*(dpml+dair+r)
sy = 2*(dpml+dair+2*r+g)
cell_size = mp.Vector3(sx,sy,0)

dpad = 0.5*dair
dbox = 0.5*g+2*r+dpad
vol = mp.Volume(center=mp.Vector3(0,0), size=mp.Vector3(2*(r+0.5*dpad),2*(0.5*g+2*r+0.5*dpad)))

# is_integrated=True necessary for any planewave source extending into PML                                                                                     
sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     center=mp.Vector3(-0.5*sx+dpml),
                     size=mp.Vector3(0,sy,0),
                     component=mp.Hz)]

# overlay for images
geometry =  [mp.Sphere(material=mp.Medium(epsilon=5),
                   center=mp.Vector3(0,r+0.5*g,0),
                   radius=r),
                mp.Sphere(material=mp.Medium(epsilon=5),
                   center=mp.Vector3(0,-r-0.5*g,0),
                   radius=r)]    

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    filename_prefix="data",
                    default_material=mp.Medium(index=args.surround),
                    geometry=geometry)

sim.run(until=1)

eps_data = sim.get_array(vol=vol, component=mp.Dielectric)
eps_data = eps_data.transpose()
np.savez("dielectric.npz",eps_data=eps_data)

# normalization run i.e., without nanoparticles
sim.reset_meep()

geometry =   [] 

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    filename_prefix="",
                    default_material=mp.Medium(index=args.surround),
                    geometry=geometry)

sim.run(mp.at_every(0.5 / (frq_cen + dfrq * 0.5), mp.in_volume(vol, mp.to_appended("{}/norm".format(save_folder), mp.output_efield_x, mp.output_efield_y, mp.output_efield_z))),
        until_after_sources=100)
        
sim.reset_meep()

# now with the nanoparticles
geometry =   [mp.Sphere(material=material,
                   center=mp.Vector3(0,r+0.5*g,0),
                   radius=r),
                mp.Sphere(material=material,
                   center=mp.Vector3(0,-r-0.5*g,0),
                   radius=r)]    

sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=mp.Vector3(),
                    sources=sources,
                    filename_prefix="",
                    geometry=geometry)

sim.run(mp.at_every(0.5 / (frq_cen + dfrq * 0.5), mp.in_volume(vol, mp.to_appended("{}/data".format(save_folder), mp.output_efield_x, mp.output_efield_y, mp.output_efield_z))),
        until_after_sources=100)
        