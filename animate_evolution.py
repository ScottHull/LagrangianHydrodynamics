import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import moviepy.editor as mpy

from src import plot

plt.style.use('seaborn-colorblind')

min_iteration = 0
max_iteration = 3000
output_directory = "jet_test_outputs"
to_dir = "animate"
ic_filename = output_directory.split("/")[-1] + ".txt"

ic_df = pd.read_csv(ic_filename, header=None, index_col=0)

if os.path.exists(to_dir):
    shutil.rmtree(to_dir)
os.mkdir(to_dir)

r_0 = float(ic_df[1]['r_0'])  # initial planet radius
T_0 = float(ic_df[1]['T_0'])
P_0 = float(ic_df[1]['P_0'])
rho_0 = float(ic_df[1]['rho_0'])
c_s_0 = float(ic_df[1]['c_s_0'])
vesc = float(ic_df[1]['v_esc'])

fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
dir_content = os.listdir(output_directory)
for index, i in enumerate(range(min_iteration, max_iteration + 1)):
    if "{}.csv".format(i) in dir_content:
        plot.plot_time(
            output_path=output_directory,
            iteration=i,
            fig=fig,
            ax_density=ax_density,
            ax_pressure=ax_pressure,
            ax_velocity=ax_velocity,
            ax_temperature=ax_temperature,
            r_0=r_0,
            rho_0=rho_0,
            P_0=P_0,
            vesc=vesc,
            T_0=T_0,
            fig_path="",
            color=colors[index],
            # min_x=0.8,
            # max_x=1.005
        )
        plt.savefig(to_dir + f"/{i}.png", format='png')

def animate(start_time, end_time, interval, path, filename="animation.mp4", fps=30, reverse=False):
    dir_contents = os.listdir(path)
    frames = [path + "/{}.png".format(time) for time in np.arange(start_time, end_time + interval, interval)
              if "{}.png".format(time) in dir_contents]
    if reverse:
        frames = list(reversed(frames))
    animation = mpy.ImageSequenceClip(frames, fps=fps, load_images=True)
    animation.write_videofile(filename, fps=fps)


animate(min_iteration, max_iteration, 1, to_dir, filename="animation.mp4", fps=30, reverse=False)

