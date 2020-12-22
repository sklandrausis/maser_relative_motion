import sys

from matplotlib import rc, cm
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np

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

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    group_number = 5

    fig, ax = plt.subplots(nrows=2, ncols=len(input_files), figsize=(16, 16))

    for index in range(0, len(input_files)):
        title = input_files[index].split(".")[0].upper() + "-" + dates[input_files[index].split(".")[0]]
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"

        velocity = np.empty(0)
        intensity = np.empty(0)
        ra = np.empty(0)
        dec = np.empty(0)
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, integral_intensity_tmp, ra_tmp, dec_tmp = np.loadtxt(input_file, unpack=True)
        for i in range(0, len(channel_tmp)):
            if group_tmp[i] == group_number:
                velocity = np.append(velocity, velocity_tmp[i])
                intensity = np.append(intensity, intensity_tmp[i])
                ra = np.append(ra, ra_tmp[i])
                dec = np.append(dec, dec_tmp[i])

        dv = (velocity.max() - velocity.min())
        vm = velocity.min()
        vx = velocity.max()

        ax[0][0].set_ylabel('Flux density (Jy)', fontsize=12)

        for i in range(len(velocity) - 1):
            if velocity[i] < vm or velocity[i] > vx:
                c = (0, 0, 0)
            else:
                c = cm.jet((velocity[i] - vm) / dv, 1)

            ax[0][index].scatter((velocity[i], velocity[i + 1]), (intensity[i], intensity[i + 1]), c=c, lw=2)
            ax[0][index].set_xlim(min(velocity) - 0.5, max(velocity) + 0.5)
            ax[0][index].xaxis.set_minor_locator(minorLocatorvel)
            ax[0][index].set_title(title, size=12)
            ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)', fontsize=12)

            rel = []
            ax[1][0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)
            for i in range(len(ra)):
                el = Circle((ra[i], dec[i]), radius=0.1 * np.sqrt(intensity[i]), angle=0, lw=2)
                ax[1][index].add_artist(el)
                c = cm.jet((velocity[i] - vm) / dv, 1)
                el.set_facecolor(c)
                rel.append([ra[i], dec[i], velocity[i]])
            #ax[1][index].add_artist(
                #Circle((285, -200), radius=10, angle=0, edgecolor='black', facecolor='white', alpha=0.9))
            #ax[1][index].annotate('1 Jy beam$^{-1}$', [275, -200], fontsize=12)
            ax[1][index].set_aspect("equal", adjustable='box')
            ax[1][index].set_xlim(-162, -150)
            ax[1][index].set_ylim(98, 110)
            ax[1][index].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
            ax[1][index].xaxis.set_minor_locator(minorLocatorx)
            ax[1][index].yaxis.set_minor_locator(minorLocatory)
            ax[1][index].invert_xaxis()

    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0, wspace=0.18, hspace=0, left=0.05, right=0.99)
    plt.show()



if __name__ == "__main__":
    main()
    sys.exit(0)
