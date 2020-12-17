import os
import sys
from random import random

import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from parsers.configparser_ import ConfigParser


def gauss(x, *p):
    a, b, c = p
    return a*np.exp(-(x-b)**2*np.log(2)/(c**2))


def DGauss(x, a1, b1, c1, a2, b2, c2):
    return a1*np.exp(-(x-b1)**2*np.log(2)/c1**2) + \
           a2*np.exp(-(x-b2)**2*np.log(2)/c2**2)


def TGauss(x, a1, b1, c1, a2, b2, c2, a3, b3, c3):
    return a1*np.exp(-(x-b1)**2*np.log(2)/c1**2) + \
          a2*np.exp(-(x-b2)**2*np.log(2)/c2**2) + \
          a3*np.exp(-(x-b3)**2*np.log(2)/c3**2)


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
    dpi = 300
    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    files = os.listdir("groups")
    files_in_order = []

    for fo in file_order:
        for f in files:
            if fo.split(".")[0] == f.split(".")[0]:
                files_in_order.append(f)

    minorLocatorvel = MultipleLocator(0.5)
    fig, ax = plt.subplots(nrows=len(file_order), ncols=1, figsize=(11.7, 8.3), dpi=dpi, sharex="all")

    colors = []
    group_index = []
    max_vel = []
    min_vel = []
    print("\hline")
    for file in files_in_order:
        index = files_in_order.index(file)
        title = file.split(".")[0].upper() + "-" + dates[file.split(".")[0]]
        group, channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt("groups/" + file, unpack=True)
        groups = list(set(group))
        max_vel.append(max(velocity))
        min_vel.append(min(velocity))

        data = dict()
        for g in groups:
            if int(g) not in group_index:
                group_index.append(int(g))
                colors.append((random(), random(), random()))

            data[g] = dict()
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

            data[g]["vel"] = vel
            data[g]["inten"] = inten
            data[g]["ra"] = ra_
            data[g]["dec"] = dec_

        for g in groups:

            vel = data[g]["vel"]
            inten = data[g]["inten"]
            ra_ = data[g]["ra"]
            dec_ = data[g]["dec"]

            group_len = len(inten)
            p0 = [max(inten), min(vel) + 0.5 * (max(vel) - min(vel)), 0.2]
            color = colors[int(g)]

            try:

                if group_len >= 3:
                    coeff, var_matrix = curve_fit(gauss, vel, inten, p0=p0)
                    for k in range(0, len(vel)):
                        print("{:3} & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f}\\\\".format(int(g), ra_[k], dec_[k], vel[k], coeff[0], inten[k], coeff[1]))
                    q = np.linspace(min(vel), max(vel), 1000)
                    hist_fit = gauss(q, *coeff)

                    ax[index].plot(q, hist_fit, 'k')

                ax[index].scatter(vel, inten, c=np.array([color]))

            except:
                pass

        ax[index].set_ylabel('Flux density [Jy]', fontsize=12)
        ax[index].text(-5.5, 5, title, size=12)
        ax[index].xaxis.set_minor_locator(minorLocatorvel)

    for file in files_in_order:
        index = files_in_order.index(file)
        ax[index].set_xlim(min(min_vel) - 0.1, max(max_vel) + 0.1)
    ax[-1].set_xlabel('$V_{\\rm LSR}$ [km s$^{-1}$]', fontsize=12)

    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0.06, wspace=0, hspace=0.05, left=0.05, right=0.99)
    #plt.show()
    plt.savefig("gauss.eps", papertype='a4', orientation='portrait', format='eps', dpi=dpi)
    print("\hline")


if __name__ == "__main__":
    main()
    sys.exit(0)
