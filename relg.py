import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator

from parsers.configparser_ import ConfigParser


def check_if_group_is_in_file(file, group):
    input_file = "groups/" + "/" + file
    group_nr = np.loadtxt(input_file, unpack=True, usecols=0)

    if group not in group_nr:
        return False
    else:
        return True


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
    minorLocatorx = MultipleLocator(20)
    minorLocatory = MultipleLocator(20)
    minorLocatorvel = MultipleLocator(1)

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    bad_files = []
    for file in input_files:
        if not check_if_group_is_in_file(file.split(".")[0] + ".groups", group_number):
            bad_files.append(file)
            del dates[file.split(".")[0]]
    input_files = [file for file in input_files if file not in bad_files]

    data = dict()
    max_intensity = []
    ch_for_max_intensity = None
    for index in range(0, len(input_files)):
        epoch = input_files[index].split(".")[0]
        data[epoch] = dict()
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        intensity = np.empty(0)
        channels = np.empty(0)
        ra = np.empty(0)
        dec = np.empty(0)
        velocity = np.empty(0)
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, ra_tmp, dec_tmp = \
            np.loadtxt(input_file, unpack=True, usecols=(0, 1, 2, 3, 5, 6))
        for i in range(0, len(channel_tmp)):
            if group_tmp[i] == int(group_number):
                intensity = np.append(intensity, intensity_tmp[i])
                channels = np.append(channels, channel_tmp[i])
                ra = np.append(ra, ra_tmp[i])
                dec = np.append(dec, dec_tmp[i])
                velocity = np.append(velocity, velocity_tmp[i])

        max_intensity.append(max(intensity))
        data[epoch]["index_for_max_intensity"] = np.where(intensity == max(intensity))[0][0]
        data[epoch]["intensity"] = intensity
        data[epoch]["channels"] = channels
        data[epoch]["velocity"] = velocity
        data[epoch]["ra"] = ra
        data[epoch]["dec"] = dec

    max_max_intensity = max(max_intensity)
    epoch_index_for_max_intensity = max_intensity.index(max_max_intensity)
    epochs = list(data.keys())
    epoch_with_max_intensity = epochs[epoch_index_for_max_intensity]
    print("epoch with max intensity", epoch_with_max_intensity)
    reference_index = data[epoch_with_max_intensity]["index_for_max_intensity"]
    references_ra = data[epoch_with_max_intensity]["ra"][reference_index]
    references_dec = data[epoch_with_max_intensity]["dec"][reference_index]
    references_velocity = data[epoch_with_max_intensity]["velocity"][reference_index]

    print("references ra", references_ra, "references dec", references_dec, "references velocity", references_velocity)

    fig, ax = plt.subplots(nrows=2, ncols=len(input_files), figsize=(16, 16), dpi=90)
    for epoch in data:
        data[epoch]["ra"] -= references_ra
        data[epoch]["dec"] -= references_dec

        if epoch != epoch_with_max_intensity:
            data[epoch]["velocity"] -= references_velocity

        velocity = data[epoch]["velocity"]
        v_min = min(velocity)
        v_max = max(velocity)
        intensity = data[epoch]["intensity"]

        index = list(data.keys()).index(epoch)
        ra = data[epoch]["ra"]
        dec = data[epoch]["dec"]
        title = epoch.upper() + "-" + dates[epoch]
        for i in range(len(velocity) - 1):
            if velocity[i] < v_min or velocity[i] > v_max:
                c = (0, 0, 0)
            else:
                c = cm.turbo((velocity[i] - v_min) / (v_max - v_min), 1)

            ax[0][index].scatter((velocity[i], velocity[i + 1]), (intensity[i], intensity[i + 1]), color=c, lw=2)
            ax[0][index].set_xlim(min(velocity) - 0.5, max(velocity) + 0.5)
            ax[0][index].xaxis.set_minor_locator(minorLocatorvel)
            ax[0][index].set_title(title)
            ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)')

            el = Circle((ra[i], dec[i]), radius=0.1 * np.sqrt(intensity[i]), angle=0, lw=2)
            el.set_facecolor(c)
            ax[1][index].add_artist(el)

        coord_range = max(max(ra) - min(ra), max(dec) - min(dec))
        ax[0][index].set_ylim((min(intensity)) - 0.1, (max(intensity) + 0.1))
        ax[1][index].set_aspect("equal", adjustable='box')
        ax[1][index].set_xlim(-10, 10)
        ax[1][index].set_ylim(-10, 10)
        ax[1][index].set_xlabel('$\\Delta$ RA (mas)')
        ax[1][index].xaxis.set_minor_locator(minorLocatorx)
        ax[1][index].yaxis.set_minor_locator(minorLocatory)
        ax[1][index].invert_xaxis()

    ax[0][0].set_ylabel('Flux density (Jy)')
    ax[1][0].set_ylabel('$\\Delta$ Dec (mas)')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    args = parser.parse_args()
    main(args.group_number)
    sys.exit(0)