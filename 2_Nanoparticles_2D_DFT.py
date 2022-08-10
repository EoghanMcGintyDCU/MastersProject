import meep as mp
import numpy as np
import os
from material import Au, Ag, Cu
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#scale is 100nm
wvl_min = 3
wvl_max = 10

frq_min = 1/wvl_max
frq_max = 1/wvl_min
frq_cen = 0.5*(frq_min+frq_max) 
dfrq = frq_max-frq_min
nfrq = 100
                                                                                  
resolution = 80

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
                                ],
                "cone" :       [mp.Cone(material=material,
                                        center=mp.Vector3(0,r+0.5*g,0),
                                        radius=r,
                                        height=2*r,
                                        axis=mp.Vector3(0,-1,0)),
                                mp.Cone(material=material,
                                        center=mp.Vector3(0,-r-0.5*g,0),
                                        radius=r,
                                        height=2*r,
                                        axis=mp.Vector3(0,1,0))
                                ]
                }
        
        geometry =  shape_map["sphere"]
        if shape != "sphere"  and shape in shape_map:
                geometry = shape_map[shape]
                print("Shape: {}".format(parse_args.shape))
        else:
                print('Shape: Defaulting to sphere')

        return geometry

material_map = {"Au": Au, "Ag" : Ag,  "Cu": Cu, "Eps7" : mp.Medium(epsilon=7)}
component_map = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez" : mp.Ez,  "Hx": mp.Hx, "Hy" : mp.Hy, "Hz" : mp.Hz}

parser = argparse.ArgumentParser(description="Perform FDTD on 2 Nanoparticles")
parser.add_argument("-m", "--material", type=str, help="Material of nanoparticles", default="Au")
parser.add_argument("-f", "--folder", type=str, help="Folder to save results", default="./")
parser.add_argument("-r", "--radius", type=float, help="Radius of nanoparticles (nm)", default=50)
parser.add_argument("-g", "--gap", type=float, help="Gap size between particles (nm)", default=5)
parser.add_argument("-i", "--index", type=float, help="Refractive index of surrounding medium", default=1)
parser.add_argument("-c", "--component", type=str, help="Source field component i.e., Ex, Hz etc.", default="Hz")
parser.add_argument("-s", "--shape", type=str, help="Shape of nanoparticle (sphere, rod)", default="sphere")

parse_args = parser.parse_args()

material = Au
if parse_args.material != "Au" and parse_args.material in material_map:
        material = material_map[parse_args.material]
        print("Material: {}".format(parse_args.material))
else:
        print('Material: Defaulting to Au')

save_folder = parse_args.folder

if not os.path.isdir(save_folder):
        print('save folder does not exist')
        exit() 

r = parse_args.radius / 100   # radius of shape     
g = parse_args.gap / 100      # distance between shape               

index = parse_args.index

component = mp.Hz
if parse_args.component != "Hz" and parse_args.component in component_map:
        component = component_map[parse_args.component]
        print("Component: {}".format(parse_args.component))
else:
        print('Component: Defaulting to Hz')

geometry = GetShape(r,g,parse_args.shape,material)

dpml = 0.5*wvl_max
dair = 0.5*wvl_max

pml_layers = [mp.PML(thickness=dpml)]

sx = 2*(dpml+dair+r)
sy = 2*(dpml+dair+2*r+g)
cell_size = mp.Vector3(sx,sy,0)

dpadx = 1.1*r
dpady = r
vol = mp.Volume(center=mp.Vector3(0,0), size=mp.Vector3(2*(r+dpadx),2*(0.5*g+2*r+dpady)))

symmetries = [mp.Mirror(mp.Y,phase=-1)]

# is_integrated=True necessary for any planewave source extending into PML                                                                                     
sources = [mp.Source(mp.GaussianSource(frq_cen,fwidth=dfrq,is_integrated=True),
                     center=mp.Vector3(-0.5*sx+dpml),
                     size=mp.Vector3(0,sy,0),
                     component=component)]

# Generating overlay for images and computational cell makeup
geometry_overlay = GetShape(r,g,parse_args.shape,mp.Medium(epsilon=5))

sim = mp.Simulation(    resolution=resolution,
                        sources=sources,
                        cell_size=cell_size,
                        boundary_layers=pml_layers,
                        geometry=geometry_overlay)

sim.add_dft_fields([mp.Ex,mp.Ey,mp.Ez], frq_cen, dfrq, nfrq, where=vol)

sim.run(until=1)

sim.plot2D()
plt.savefig('{}/cell.png'.format(save_folder))
plt.close()

overlay_data = sim.get_array(vol=vol, component=mp.Dielectric)
overlay_data = overlay_data.transpose()

sim.reset_meep()

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
                    symmetries=symmetries,         
                    geometry=empty_geometry)

dft_norm = sim.add_dft_fields([mp.Ex,mp.Ey,mp.Ez], frq_cen, dfrq, nfrq, where=vol)

sim.run(until_after_sources=100)

dft_shape = sim.get_dft_array(dft_norm, mp.Ey, 0).shape
dft_e_fields_norm = np.zeros((dft_shape[0], dft_shape[1], nfrq),dtype=np.complex128)
for i in range(nfrq):
        dft_e_fields_norm[:,:,i] += sim.get_dft_array(dft_norm, mp.Ex, i)        
        dft_e_fields_norm[:,:,i] += sim.get_dft_array(dft_norm, mp.Ey, i)
        dft_e_fields_norm[:,:,i] += sim.get_dft_array(dft_norm, mp.Ez, i)
        
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
                    symmetries=symmetries,
                    geometry=geometry)

dft = sim.add_dft_fields([mp.Ex,mp.Ey,mp.Ez], frq_cen, dfrq, nfrq, where=vol)

sim.run(until_after_sources=100)

dft_e_fields = np.zeros((dft_shape[0], dft_shape[1], nfrq),dtype=np.complex128)
for i in range(nfrq):
        dft_e_fields[:,:,i] += sim.get_dft_array(dft, mp.Ex, i)
        dft_e_fields[:,:,i] += sim.get_dft_array(dft, mp.Ey, i)
        dft_e_fields[:,:,i] += sim.get_dft_array(dft, mp.Ez, i)

dft_e_fields_norm = np.abs(dft_e_fields_norm)
norm_power = np.square(dft_e_fields_norm)

dft_e_fields = np.abs(dft_e_fields)
power = np.square(dft_e_fields)

enhancement = power / norm_power

freqs = np.linspace(frq_min, frq_max, nfrq) 

hotspot_enhancement = np.zeros_like(freqs)

for i in range(nfrq):
    hotspot_enhancement[i] = np.max(enhancement[:,:,i])

plt.figure()
plt.plot(100 / freqs[:], hotspot_enhancement[:])
plt.xlabel("Wavelength, nm")
plt.ylabel("$|\\vec{E}$/$\\vec{E_0}|^2$")
plt.grid()
plt.savefig("{}/EnhancementSpectra.png".format(save_folder),dpi=150,bbox_inches='tight')
plt.close()

np.savez("{}/hotspot_enhancement.npz".format(save_folder),hotspot_enhancement=hotspot_enhancement)

max_enhancement = np.argmax(hotspot_enhancement[:])
max_f = freqs[max_enhancement]
max_wv = 100/max_f

plt.figure()
plt.title("Enhancement {}nm".format(max_wv))
plt.imshow(overlay_data, interpolation='spline36', cmap='binary')
plt.imshow(enhancement[:,:,max_enhancement].transpose(), interpolation='spline36', cmap='plasma', alpha=0.9)
plt.axis('off')
plt.colorbar()
plt.savefig("{}/MaxEnhancement.png".format(save_folder),dpi=150,bbox_inches='tight')
plt.close()
