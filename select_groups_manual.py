import sys
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


def main():
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    minorLocatorx = MultipleLocator(20)
    minorLocatory = MultipleLocator(20)
    minorLocatorvel = MultipleLocator(1)
    file = "data_files3/ea063.out"
    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}
    title = file.split("/")[1].split(".")[0].upper() + "-" + dates[file.split("/")[1].split(".")[0]]
    channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt(file, unpack=True)
    velocity = velocity/1000

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16, 16))
    ax[0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)

    scatter1 = ax[0].scatter(ra, dec, picker=5)
    ax[0].annotate('1 Jy beam$^{-1}$', [275, -200], fontsize=12)
    ax[0].set_aspect("equal", adjustable='box')
    ax[0].set_xlim(-200, 200)
    ax[0].set_ylim(-50, 345)
    ax[0].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
    ax[0].xaxis.set_minor_locator(minorLocatorx)
    ax[0].yaxis.set_minor_locator(minorLocatory)
    ax[0].invert_xaxis()
    ax[0].xaxis.set_minor_locator(minorLocatorvel)
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
    scatter2 = ax[1].scatter(velocity, intensity, picker=5)
    ax[1].set_xlabel('$V_{\\rm LSR}$ [km s$^{-1}$]', fontsize=12)
    ax[1].xaxis.set_minor_locator(minorLocatorvel)
    ax[1].set_ylabel('Flux density [Jy]', fontsize=12)
    cursor2 = mplcursors.cursor(scatter2, hover=True, highlight=True)
    labels2 = create_labels2(ra, dec)
    cursor2.connect("add", lambda sel: sel.annotation.set_text(labels2[sel.target.index]))

    colors = ["b", "g", "r", "c", "m", "y", "k"]
    selected_points = []

    def onpick1(event):
        global group_index
        ind = event.ind[0]
        if [group_index, channel[ind], velocity[ind], intensity[ind], integral_intensity[ind], ra[ind], dec[ind]] not in selected_points:
            print("qqqq")
            selected_points.append([group_index, channel[ind], velocity[ind], intensity[ind], integral_intensity[ind], ra[ind], dec[ind]])
            ax[1].plot(velocity[ind], intensity[ind], colors[group_index] + "x", markersize=10)
            ax[0].plot(ra[ind], dec[ind], colors[group_index] + "x", markersize=10)
            groups[group_index].append([group_index, channel[ind], velocity[ind], intensity[ind], integral_intensity[ind], ra[ind], dec[ind]])
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
            else:
                print("Wrong number")

        else:
            print("Not digit")

    fig.canvas.mpl_connect('key_press_event', press)
    fig.canvas.mpl_connect('pick_event', onpick1)

    plt.tight_layout()

    plt.show()
    print("Output file ", file.split("/")[1].split(".")[0] + ".groups")
    with open(file.split("/")[1].split(".")[0] + ".groups", "w") as output_file:
        for group in groups:
            for chann in group:
                for i in chann:
                    output_file.write(str(i) + " ")
                output_file.write("\n")


if __name__ == "__main__":
    main()
    sys.exit(0)
