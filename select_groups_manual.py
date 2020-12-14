import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator

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
    dv = (velocity.max() - velocity.min())
    vm = velocity.min()

    rel = []
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16, 16))
    ax[0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)

    ax[0].scatter(ra, dec)
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

    groups = [[]]
    global group_index
    group_index = 0
    group_indexes = [group_index]
    ax[1].plot(velocity, intensity, ls="", marker="o", markersize=5, markerfacecolor="white", markeredgewidth=1, picker=5)
    ax[1].set_xlabel('$V_{\\rm LSR}$ [km s$^{-1}$]', fontsize=12)
    ax[1].xaxis.set_minor_locator(minorLocatorvel)
    ax[1].set_ylabel('Flux density [Jy]', fontsize=12)

    def onpick1(event):
        global group_index
        ind = event.ind[0]
        this_line = event.artist
        xdata = this_line.get_xdata()[ind]
        ydata = this_line.get_ydata()[ind]
        ax[1].plot(xdata, ydata, "rx", markersize=10)
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
    #plt.subplots_adjust(top=0.97, bottom=0, wspace=0.18, hspace=0, left=0.05, right=0.99)

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
