import sys
import os
from random import random

import matplotlib.pyplot as plt
import numpy as np


def main():
    files = os.listdir("groups")

    data = dict()
    global_groups = []
    global_epochs = []

    for file in files:
        file_name = "groups/" + file
        epoch = file.split(".")[0]
        group, channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt(file_name, unpack=True)
        groups = list(set(group))
        global_groups += groups
        global_epochs.append(epoch)

        for g in groups:
            if g not in data.keys():
                data[g] = dict()
                data[g]["vel"] = []
                data[g]["ra"] = []
                data[g]["dec"] = []
                data[g]["epoch"] = []

            vel = []
            inten = []
            ra_ = []
            dec_ = []
            epoch_ = []
            for ch in range(0, len(channel)):
                if group[ch] == g:
                    vel.append(velocity[ch])
                    inten.append(intensity[ch])
                    ra_.append(ra[ch])
                    dec_.append(dec[ch])
                    epoch_.append(epoch)

            data[g]["vel"] += vel
            data[g]["ra"] += ra_
            data[g]["dec"] += dec_
            data[g]["epoch"] += epoch_

    global_groups = list(set(global_groups))
    global_epochs = list(set(global_epochs))
    epoch_colors = []
    legends = []
    for ge in global_epochs:
        color = (random(), random(), random())
        epoch_colors.append(color)
        legends.append(ge + "\n" + str(color))

    for g in global_groups:
        vel = data[g]["vel"]
        ra_ = data[g]["ra"]
        dec_ = data[g]["dec"]
        epoch_ = data[g]["epoch"]

        colors = []

        for e in epoch_:
            colors.append(epoch_colors[global_epochs.index(e)])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
        ax.set_zlabel("Velocity")
        scatter = ax.scatter(ra_, dec_, vel, color=colors)
        ax.set_title("Group index is " + str(int(g)))
        ax.legend(handles=[scatter], loc="lower left", title=legends)

    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
