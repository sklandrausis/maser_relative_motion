import sys
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import rc
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection
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
    relative_motion_path = get_configs("paths", "relative_motion")
    input_file = relative_motion_path + "linearity_errors_fitted_cm.dat"
    input_data = ascii.read(input_file)
    v1 = input_data["vel"]
    f = input_data["f"]
    x1 = input_data["x1"]
    y1 = input_data["y1"]
    x2 = input_data["x2"]
    y2 = input_data["y2"]
    xlong2 = input_data["xlong2"]
    ylong2 = input_data["ylong2"]
    errxminlong = input_data["errxminlong"]
    erryminlong = input_data["erryminlong"]
    errxmaxlong = input_data["errxmaxlong"]
    errymaxlong = input_data["errymaxlong"]

    dv = (v1.max() - v1.min())
    v1mi = v1.min()

    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    plt.figure(figsize=(7, 5), dpi=100)
    ax1 = plt.subplot( 111, aspect='equal' )

    ls = []
    lsreal = []
    for i in range(len(x1)):
        leng = np.sqrt((x1[i] - xlong2[i]) ** 2 + (y1[i] - ylong2[i]) ** 2)
        lengreal = np.sqrt((x1[i] - x2[i]) ** 2 + (y1[i] - y2[i]) ** 2)
        lsreal.append([lengreal])
        ls.append([leng])
        ax1.annotate("", xy=(x1[i], y1[i]), xycoords='data', xytext=(errxmaxlong[i], errymaxlong[i]), textcoords='data',
                     arrowprops=dict( arrowstyle="-", color="grey", connectionstyle="arc3", alpha=0.5))
        ax1.annotate("", xy=(x1[i], y1[i]), xycoords='data', xytext=(errxminlong[i], erryminlong[i]), textcoords='data',
                     arrowprops=dict( arrowstyle="-", color="grey", connectionstyle="arc3", alpha=0.5))
        ax1.annotate("", xy=(errxminlong[i], erryminlong[i]), xycoords='data', xytext=(errxmaxlong[i], errymaxlong[i]),
                     textcoords='data', arrowprops=dict(arrowstyle="-", color="grey", connectionstyle="arc3",
                                                        alpha=0.5))
        ax1.annotate("", xy=(x1[i], y1[i]), xycoords='data', xytext=(xlong2[i], ylong2[i]), textcoords='data', size=7,
                     arrowprops=dict(arrowstyle="<-", color="black", connectionstyle="arc3"))
        el = Circle((x1[i], y1[i]), radius=3 * np.log10(f[i] * 1000.), angle=0, lw=0.5)
        ax1.add_artist(el)
        c = cm.jet((v1[i] - v1mi) / dv, 1)
        el.set_facecolor(c)

    patches = []
    colors = v1
    p = PatchCollection(patches, cmap=cm.jet)
    p.set_array(np.array(colors))
    ax1.add_collection(p)
    plt.colorbar(p)
    als = np.array(ls)

    ax1.annotate("", xy=(50, -210), xycoords='data', xytext=(50 + als.max() / 2, -210), textcoords='data', size=7,
                 arrowprops=dict(arrowstyle="<-", connectionstyle="arc3"))
    plt.text(110, -195, "0.3 mas yr$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='k')
    plt.text(100, -230, "6 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='k')
    ax1.annotate("", xy=(150, 210), xycoords='data', xytext=(50, 210), textcoords='data', size=7,
                 arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))

    #6 mas x 4.18 kpc = 25.08 AU = 376.2e7
    plt.text(125, 220, "418 AU", size=8, rotation=0.0, ha="left", va="center", color='k')
    plt.plot(0, 0, marker='+', c='k')

    ax1.set_xlim(200, -300)
    ax1.set_ylim(-250, 250)
    plt.xlabel('$\\Delta$RA (mas)', fontsize=12)
    plt.ylabel('$\\Delta$Dec (mas)', fontsize=12)
    plt.text(-275, 270, "$V_{LSR}$ (km s$^{-1}$)", size=12, rotation=0.0, ha="left", va="center", color='k')
    plt.title("G78", size=12)

    alsreal = np.array(lsreal)
    print("max lenght is ", als.max())
    print("max lenght is ", alsreal.max())
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
