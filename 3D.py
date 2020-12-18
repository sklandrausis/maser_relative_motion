import sys
import os
import matplotlib.pyplot as plt
import numpy as np


def main():
    files = os.listdir("groups")

    data = dict()
    global_groups = []

    for file in files:
        file_name = "groups/" + file
        group, channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt(file_name, unpack=True)
        groups = list(set(group))
        global_groups += groups

        for g in groups:
            if g not in data.keys():
                data[g] = dict()
                data[g]["vel"] = []
                data[g]["ra"] = []
                data[g]["dec"] = []

            vel = []
            inten = []
            ra_ = []
            dec_ = []
            for ch in range(0, len(channel)):
                if group[ch] == g:
                    vel.append(velocity[ch])
                    inten.append(intensity[ch])
                    ra_.append(ra[ch])
                    dec_.append(dec[ch])

            data[g]["vel"] += vel
            data[g]["ra"] += ra_
            data[g]["dec"] += dec_

    global_groups = list(set(global_groups))
    for g in global_groups:
        vel = data[g]["vel"]
        ra_ = data[g]["ra"]
        dec_ = data[g]["dec"]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
        ax.set_zlabel("Velocity")
        ax.scatter(ra_, dec_, vel)

    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
