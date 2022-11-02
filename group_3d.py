import sys
import os

import matplotlib.pyplot as plt
import numpy as np

from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """
    :param section: configuration file section
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main():
    groups_file_path = get_configs("paths", "groups")
    files = os.listdir(groups_file_path)

    data = dict()
    global_groups = []
    global_epochs = []

    for file in files:
        file_name = groups_file_path + file
        epoch = file.split(".")[0]
        group, channel, velocity, ra, dec = np.loadtxt(file_name, unpack=True, usecols=(0, 1, 2, 5, 6))
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
            ra_ = []
            dec_ = []
            epoch_ = []
            for ch in range(0, len(channel)):
                if group[ch] == g:
                    vel.append(velocity[ch])
                    ra_.append(ra[ch])
                    dec_.append(dec[ch])
                    epoch_.append(epoch)

            data[g]["vel"] += vel
            data[g]["ra"] += ra_
            data[g]["dec"] += dec_
            data[g]["epoch"] += epoch_

    global_groups = list(set(global_groups))
    global_epochs = list(set(global_epochs))

    for g in global_groups:
        vel = data[g]["vel"]
        ra_ = data[g]["ra"]
        dec_ = data[g]["dec"]
        epoch_ = data[g]["epoch"]

        vel_tmp = []
        ra_tmp = []
        dec_tmp = []

        for ge in global_epochs:
            vel_tmp.append([])
            ra_tmp.append([])
            dec_tmp.append([])

        for epoch_index in range(0, len(epoch_)):
            index = global_epochs.index(epoch_[epoch_index])
            vel_tmp[index].append(vel[epoch_index])
            ra_tmp[index].append(ra_[epoch_index])
            dec_tmp[index].append(dec_[epoch_index])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
        ax.set_zlabel("Velocity")
        for ep in global_epochs:
            index = global_epochs.index(ep)
            ax.scatter(ra_tmp[index], dec_tmp[index], vel_tmp[index], label=ep)

        ax.set_title("Group index is " + str(int(g)))
        ax.legend()

    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
