import os
import imageio

png_dir ="../../Figures/Arctic/Predictions/"
images = []
for subdir, dirs, files in os.walk(png_dir):
    for file in files:
        file_path = os.path.join(subdir, file)
        if file_path.endswith("multi.png"):
            images.append(imageio.imread(file_path))
kargs = { 'duration': 0.3, 'loop':1 }
imageio.mimsave(png_dir+'/animatedForecasts.gif', images, loop=2)
imageio.mimsave(png_dir+'/animatedForecasts.mp4', images)