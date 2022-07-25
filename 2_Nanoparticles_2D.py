import meep as mp
import numpy as np
import os
from material import Au, Au_JC_visible, Ag, Cu
import argparse

#scale is 100nm
wvl_min = 2
wvl_max = 10

frq_min = 1/wvl_max
frq_max = 1/wvl_min
frq_cen = 0.5*(frq_min+frq_max) 
dfrq = frq_max-frq_min
nfrq = 100
                                                                                  
resolution = 50

def GetShape(r,g,shape,material):

        shape_map =     {
                "sphere" :     [mp.Sphere(material=material,
                                        center=mp.Vector3(0,r+0.5*g,0),
                                        radius=r),
                                mp.Sphere(material=material,
                                        center=mp.Vector3(0,-r-0.5*g,0),
                                        radius=r)
                                ], 
                "rod" :         [mp.Cylinder(material=material,
                                        center=mp.Vector3(0,r+0.5*g,0),
                                        radius=r,
                                        height=4*r,
                                        axis=mp.Vector3(1,0,0)),
                                mp.Cylinder(material=material,
                                        center=mp.Vector3(0,-r-0.5*g,0),
                                        radius=r,
                                        height=4*r,
                                        axis=mp.Vector3(1,0,0))
                                ]
                }
        
        geometry =  shape_map["sphere"]
        if shape != "sphere"  and shape in shape_map:
                geometry = shape_map[shape]
                print("Shape: {}".format(args.shape))
        else:
                print('Shape: Defaulting to sphere')

        return geometry

material_map = {"Au": Au, "Au_JC_visible": Au_JC_visible, "Ag" : Ag,  "Cu": Cu, "Eps7" : mp.Medium(epsilon=7)}
component_map = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez" : mp.Ez,  "Hx": mp.Hx, "Hy" : mp.Hy, "Hz" : mp.Hz}

parser = argparse.ArgumentParser(description="Perform FDTD on 2 Nanoparticles")
parser.add_argument("-m", "--material", type=str, help="Material of nanoparticles", default="Au")
parser.add_argument("-f", "--folder", type=str, help="Folder to save results", default="./")
parser.add_argument("-r", "--radius", type=float, help="Radius of nanoparticles (nm)", default=50)
parser.add_argument("-g", "--gap", type=float, help="Gap size between particles (nm)", default=5)
parser.add_argument("-i", "--index", type=float, help="Refractive index of surrounding medium", default=1)
parser.add_argument("-c", "--component", type=str, help="Source field component i.e., Ey, Hz etc.", default="Ey")
parser.add_argument("-s", "--shape", type=str, help="Shape of nanoparticle (sphere, rod)", default="sphere")

args = parser.parse_args()

material = Au
if args.material != "Au" and args.material in material_map:
        material = material_map[args.material]
        print("Material: {}".format(args.material))
else:
        print('Material: Defaulting to Au')

save_folder = args.folder

if not os.path.isdir(save_folder):
        print('save folder does not exist')
        exit() 

r = args.radius / 100   # radius of shape     
g = args.gap / 100      # distance between shape               

index = args.index

component = mp.Ey
if args.component != "Ey" and args.component in component_map:
        component = component_map[args.component]
        print("Component: {}".format(args.component))
else:
        print('Component: Defaulting to Ey')

geometry = GetShape(r,g,args.shape,material)

dpml = 0.5*wvl_max
dair = 0.5*wvl_max

pml_layers = [mp.PML(thickness=dpml)]

sx = 2*(dpml+dair+r)
sy = 2*(dpml+dair+2*r+g)
cell_size = mp.Vector3(sx,sy,0)

dpad = 0.5*dair
vol = mp.Volume(center=mp.Vector3(0,0), size=mp.Vector3(2*(r+0.5*dpad),2*(0.5*g+2*r+0.5*dpad)))

# Generating overlay for images
geometry_overlay = GetShape(r,g,args.shape,mp.Medium(epsilon=5))

sim = mp.Simulation(resolution=resolution,
                cell_size=cell_size,
                boundary_layers=pml_layers,
                geometry=geometry_overlay)

sim.run(until=1)

overlay_data = sim.get_array(vol=vol, component=mp.Dielectric)
overlay_data = overlay_data.transpose()
np.savez("{}/overlay.npz".format(save_folder),overlay_data=overlay_data)

sim.reset_meep()

# is_integrated=True necessary for any planewave source extending into PML                                                                                     
sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     center=mp.Vector3(-0.5*sx+dpml),
                     size=mp.Vector3(0,sy,0),
                     component=component)]

empty_geometry =   [] 

## normalization run i.e., without nanoparticles
sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    filename_prefix="",
                    default_material=mp.Medium(index=index),
                    geometry=empty_geometry)

sim.run(mp.at_every(1 / (2*frq_max), mp.in_volume(vol, mp.to_appended("{}/norm".format(save_folder), mp.output_efield_x, mp.output_efield_y, mp.output_efield_z))),
        until_after_sources=100)
        
sim.reset_meep()

## now with the nanoparticles
sim = mp.Simulation(split_chunks_evenly=False,
                    resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=mp.Vector3(),
                    sources=sources,
                    filename_prefix="",
                    default_material=mp.Medium(index=index),
                    geometry=geometry)

sim.run(mp.at_every(1 / (2*frq_max), mp.in_volume(vol, mp.to_appended("{}/data".format(save_folder), mp.output_efield_x, mp.output_efield_y, mp.output_efield_z))),
        until_after_sources=100)
        