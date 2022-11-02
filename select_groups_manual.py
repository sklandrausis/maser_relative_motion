import argparse
import sys
from random import random

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
import mplcursors

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


def main(input_file):
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    minor_locator_x = MultipleLocator(20)
    minor_locator_y = MultipleLocator(20)
    minor_locator_vel = MultipleLocator(1)
    data_file_path = get_configs("paths", "dataFiles")
    file = data_file_path + input_file
    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}
    title = input_file.split(".")[0].upper() + "-" + dates[input_file.split(".")[0]]
    channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt(file, unpack=True)
    velocity = velocity/1000

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16, 16))
    ax[0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)

    scatter1 = ax[0].scatter(ra, dec, picker=True)
    ax[0].annotate('1 Jy beam$^{-1}$', [275, -200], fontsize=12)
    ax[0].set_aspect("equal", adjustable='box')
    ax[0].set_xlim(-200, 200)
    ax[0].set_ylim(-50, 345)
    ax[0].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
    ax[0].xaxis.set_minor_locator(minor_locator_x)
    ax[0].yaxis.set_minor_locator(minor_locator_y)
    ax[0].invert_xaxis()
    ax[0].xaxis.set_minor_locator(minor_locator_vel)
    ax[0].set_title(title, size=12)

    def create_labels1(velocity, intensity):
        labels1 = []
        for i in range(0, len(velocity)):
            labels1.append("velocity is " + str(velocity[i]) + "\nintensity is " + str(intensity[i]))

        return labels1

    def create_labels2(ra, dec):
        labels2 = []
        for i in range(0, len(ra)):
            labels2.append("RA is " + str(ra[i]) + "\nDEC is " + str(dec[i]))

        return labels2

    cursor1 = mplcursors.cursor(scatter1, hover=True, highlight=True)
    labels1 = create_labels1(velocity, intensity)
    cursor1.connect("add", lambda sel: sel.annotation.set_text(labels1[sel.target.index]))

    groups = [[]]
    global group_index
    group_index = 0
    group_indexes = [group_index]
    scatter2 = ax[1].scatter(velocity, intensity, picker=True)
    ax[1].set_xlabel('$V_{\\rm LSR}$ [km s$^{-1}$]', fontsize=12)
    ax[1].xaxis.set_minor_locator(minor_locator_vel)
    ax[1].set_ylabel('Flux density [Jy]', fontsize=12)
    cursor2 = mplcursors.cursor(scatter2, hover=True, highlight=True)
    labels2 = create_labels2(ra, dec)
    cursor2.connect("add", lambda sel: sel.annotation.set_text(labels2[sel.target.index]))

    colors = [(random(), random(), random())]
    selected_points = []

    def onpick1(event):
        global group_index
        ind = event.ind[0]
        if [group_index, channel[ind], velocity[ind], intensity[ind], integral_intensity[ind], ra[ind], dec[ind]] \
                not in selected_points:
            selected_points.append([group_index, channel[ind], velocity[ind], intensity[ind], integral_intensity[ind],
                                    ra[ind], dec[ind]])
            ax[1].plot(velocity[ind], intensity[ind],"x", markersize=10, c=colors[group_index])
            ax[0].plot(ra[ind], dec[ind], "x", markersize=10, c=colors[group_index])
            groups[group_index].append([group_index, channel[ind], velocity[ind], intensity[ind],
                                        integral_intensity[ind], ra[ind], dec[ind]])
            event.canvas.draw()

    def press(event):
        global group_index
        if event.key.isdigit():
            if int(event.key) in group_indexes:
                print("group_index changed to ", event.key)
                group_index = int(event.key)
            elif int(event.key) == group_indexes[-1] + 1:
                print("group_index changed to ", event.key)
                group_index = int(event.key)
                groups.append([])
                group_indexes.append(int(event.key))
                colors.append((random(), random(), random()))
            else:
                print("Wrong number")

        else:
            if event.key == "shift":
                print("group_index changed to ", group_index + 1)
                group_index += 1
                if group_index not in group_indexes:
                    group_indexes.append(group_index)
                    groups.append([])
                    colors.append((random(), random(), random()))
            elif event.key == "alt":
                print("group_index changed to ", group_index - 1)
                group_index -= 1
                if group_index not in group_indexes:
                    group_indexes.append(group_index)
                    groups.append([])
                    colors.append((random(), random(), random()))
            else:
                print("Not digit")

    fig.canvas.mpl_connect('key_press_event', press)
    fig.canvas.mpl_connect('pick_event', onpick1)

    plt.tight_layout()

    plt.show()
    groups_file_path = get_configs("paths", "groups")
    print(groups_file_path + input_file.split(".")[0] + ".groups")
    with open(groups_file_path + input_file.split(".")[0] + ".groups", "w") as output_file:
        for group in groups:
            for chann in group:
                for i in chann:
                    output_file.write(str(i) + " ")
                output_file.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='create group file')
    parser.add_argument('input_file', type=str, help='input file')
    args = parser.parse_args()
    main(args.input_file)
    sys.exit(0)
