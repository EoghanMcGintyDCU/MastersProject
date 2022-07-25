import h5py
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse

wvl_min = 2
wvl_max = 10

frq_min = 1/wvl_max
frq_max = 1/wvl_min
nfrq = 100

parser = argparse.ArgumentParser(description="Perform Enhancement Analysis on 2 Nanoparticles FDTD output")
parser.add_argument("-f", "--folder", type=str, help="Folder to read results from", default="./")

args = parser.parse_args()

save_folder = args.folder

norm_data = h5py.File("{}/norm.h5".format(save_folder), "r")
keys = list(norm_data.keys())
norm_data = [norm_data[key] for key in keys]  # read all datasets in the hf5 file
print("Components: {}".format(keys))
print("Data {}".format(norm_data[0].shape))

l = norm_data[0].shape[2]           # length of temporal dimension
flen = l // 2                       # length of data in frequency
print("Temporal length of data {}, length of freq domain {}".format(l, flen))

print("Freq max {}".format(frq_max))
freqs = np.linspace(frq_min, frq_max, nfrq)     # change linspace to arange 0...f_nyq

data = h5py.File("{}/data.h5".format(save_folder), "r")
data = [data[key] for key in keys]              # read all datasets in the hf5 file
print("Components: {}".format(keys))
print("Data {}".format(data[0].shape))

def GetPower(time_data):
    # convert chunk of data from time domain to
    # squared frequency domain (power spectrum)    
    fdata = np.fft.fft(time_data, flen * 2, axis=-1) / l
    fdata = 2 * np.abs(fdata[:, :, 0:flen])
    power = np.square(fdata)
    return power

sh = data[0][:].shape
power = np.zeros((sh[0], sh[1], flen))
# X, Y, Z
for j in data:
    print(j.name)
    power += GetPower(j[:])

ref_power = np.zeros((sh[0], sh[1], flen))
# X, Y, Z
for j in norm_data:
    print(j.name)
    ref_power += GetPower(j[:])

enhancement = power[:, :, :] / ref_power[:, :, :]

hotspot_enhancement = np.zeros_like(freqs)

overlay_data = np.load("overlay.npz")
overlay_data = overlay_data['overlay_data']

# skip the nonphysical low-frequncy part of spectra
skip_freq = 5

for i in range(flen):
    hotspot_enhancement[i] = np.percentile(enhancement[:,:,i], 99)

plt.plot(100 / freqs[skip_freq:], hotspot_enhancement[skip_freq:])
plt.xlabel("Wavelength, nm")
plt.ylabel("$|\\vec{E}$/$\\vec{E_0}|^2$")
plt.grid()
plt.savefig("{}/EnhancementSpectra.png".format(save_folder),dpi=150,bbox_inches='tight')

max_enhancement = np.argmax(hotspot_enhancement[skip_freq:]) + skip_freq
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
