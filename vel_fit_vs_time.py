import sys
import os
import argparse
from random import random

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


def get_cloudlet_sub_files_for_group(cloudlet_sub_files, group):
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    cloudlet_sub_files_for_group = [f for f in cloudlet_sub_files if int(f.split("_")[4].replace(".", "")) == group]
    cloudlet_sub_files_for_group_tmp = [0] * len(file_order)
    for file in file_order:
        epoch = file.split(".")[0]
        for c_file in cloudlet_sub_files_for_group:
            c_epoch = c_file.split("_")[3]
            if epoch == c_epoch:
                cloudlet_sub_files_for_group_tmp[file_order.index(file)] = c_file

    cloudlet_sub_files_for_group = cloudlet_sub_files_for_group_tmp
    return cloudlet_sub_files_for_group


def main():
    cloudlet_sub_files = [file for file in os.listdir("./") if file.startswith("cloudlet_sub")]
    groups = sorted(list(set([int(f.split("_")[4].replace(".", "")) for f in cloudlet_sub_files])))
    cloudlet_sub_files_for_all_groups = {g: get_cloudlet_sub_files_for_group(cloudlet_sub_files, g) for g in groups}

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=120)

    for group in groups:
        color = [random(), random(), random()]
        cloudlet_sub_files_for_group = cloudlet_sub_files_for_all_groups[group]
        for file in cloudlet_sub_files_for_group:
            if file != 0:
                epoch = file.split("_")[3]
                date = dates[epoch]
                cloudlet_sub_data = ascii.read(file)
                vel_fit = np.array(cloudlet_sub_data["vel_fit"], dtype=np.float)
                fit_amp = np.array(cloudlet_sub_data["fit_amp"], dtype=np.float)
                main_index = np.where(fit_amp == max(fit_amp))[0][0]
                main_vel_fit = vel_fit[main_index]
                main_fit_amp = fit_amp[main_index]
                ax.scatter(date, main_vel_fit, s=10 * main_fit_amp, color=color)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='vel_fit vs time')
    args = parser.parse_args()
    main()
    sys.exit(0)