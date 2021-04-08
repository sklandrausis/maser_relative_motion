import sys
import os
import argparse

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt

from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """
    :param section: configuration file secti
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main(group_number):
    cloudlet_sub_files = [file for file in os.listdir("./") if file.startswith("cloudlet_sub")
                         if file.split("_")[-2] == str(group_number) + "."]

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]

    cloudlet_sub_files_tmp = [0] * len(cloudlet_sub_files)
    for file in file_order:
        epoch = file.split(".")[0]
        for c_file in cloudlet_sub_files:
            c_epoch = c_file.split("_")[3]
            if epoch == c_epoch:
                cloudlet_sub_files_tmp[file_order.index(file)] = c_file

    cloudlet_sub_files = cloudlet_sub_files_tmp


    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=120)
    for file in cloudlet_sub_files:
        epoch = file.split("_")[3]
        date = dates[epoch]

        cloudlet_sub_data = ascii.read(file)
        vel_fit = np.array(cloudlet_sub_data["vel_fit"], dtype=np.float)
        fit_amp = np.array(cloudlet_sub_data["fit_amp"], dtype=np.float)
        ax.scatter([date]*len(vel_fit), vel_fit, s=1000*fit_amp)


    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='vel_fit vs time')
    parser.add_argument('group_number', type=int, help='group number')
    args = parser.parse_args()
    main(args.group_number)
    sys.exit(0)