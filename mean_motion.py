import sys
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection

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
    output_file = "output/output.dat"
    output_data = ascii.read(output_file)
    output_data_headers = output_data.keys()
    velocity = output_data["vel"]
    number_of_points = len(get_configs("parameters", "fileOrder").split(","))
    velocity_range = max(velocity) - min(velocity)
    ras = []
    decs = []
    fluxs = []

    tmp = 1
    for index in range(1, number_of_points * 3 + 1):
        tmp_data = output_data[output_data_headers[index]]

        if tmp == 1:
            ras.append(tmp_data)
            tmp = 2
        elif tmp == 2:
            decs.append(tmp_data)
            tmp = 3
        else:
            fluxs.append(tmp_data)
            tmp = 1

    reference_ras = ras[0]
    reference_decs = decs[0]

    groups_indexies = []
    group_tmp = []

    for group_index in range(0, len(reference_ras) - 1):
        ra_diff = reference_ras[group_index] - reference_ras[group_index + 1]
        dec_diff = reference_decs[group_index] - reference_decs[group_index + 1]

        if abs(ra_diff) <= 1.27 and abs(dec_diff) <= 1.27:
            group_tmp.append(group_index)
            group_tmp.append(group_index + 1)

        else:
            if sorted(group_tmp) not in groups_indexies:
                groups_indexies.append(sorted(group_tmp))
            group_tmp = []

    for group_index in range(len(reference_ras)-1, 0, -1):
        ra_diff = reference_ras[group_index] - reference_ras[group_index - 1]
        dec_diff = reference_decs[group_index] - reference_decs[group_index - 1]

        if abs(ra_diff) <= 1.27 and abs(dec_diff) <= 1.27:
            group_tmp.append(group_index)
            group_tmp.append(group_index - 1)

        else:
            if sorted(group_tmp) not in groups_indexies:
                groups_indexies.append(sorted(group_tmp))
            group_tmp = []

    groups_indexies = [set(g) for g in groups_indexies if len(g) > 0]
    print("Groups ", groups_indexies)

    plt.figure()
    ax1 = plt.subplot(111, aspect='equal')
    spots_parameters = []
    vectors_parameters = []

    for group in groups_indexies:
        number_of_elements_in_group = len(group)
        sum_of_ra_diffs = []
        sum_of_dec_diffs = []
        sum_of_ras = []
        sum_of_decs = []
        lengths = []
        ls_tmp = []
        vel_for_group = [velocity[gi] for gi in group]
        sum_of_vel = sum(vel_for_group)
        flux_for_group = []

        for index in range(0, len(ras)):
            if index != 0:
                ra_diff = [ras[index][gi] - ras[0][gi] for gi in group]
                dec_diff = [decs[index][gi] - decs[0][gi] for gi in group]
                sum_ra_diff = sum(ra_diff)
                sum_dec_diff = sum(dec_diff)
                sum_of_ra_diffs.append(sum_ra_diff/number_of_elements_in_group)
                sum_of_dec_diffs.append(sum_dec_diff/number_of_elements_in_group)
                length = np.sqrt(sum_ra_diff ** 2 + sum_dec_diff ** 2)

            ra_for_group = [ras[index][gi] for gi in group]
            sum_of_ra = sum(ra_for_group)
            dec_for_group = [decs[index][gi] for gi in group]
            sum_of_dec = sum(dec_for_group)
            sum_of_ras.append(sum_of_ra)
            sum_of_decs.append(sum_of_dec)
            flux_for_group_tmp = [fluxs[index][gi] for gi in group]
            flux_for_group.extend(flux_for_group_tmp)

        coords = (sum_of_ras[0]/number_of_elements_in_group, sum_of_decs[0]/number_of_elements_in_group)
        spot_parameters = {"coords": coords, "vel": sum_of_vel / number_of_elements_in_group,
                           "radius": 3 * np.log10(max(flux_for_group) * 1000.)}
        spots_parameters.append(spot_parameters)
        vector_parameters = {"sum_of_ra_diffs": sum_of_ra_diffs, "sum_of_dec_diffs": sum_of_dec_diffs}
        vectors_parameters.append(vector_parameters)

    vector_colors = ["black", "grey", "blue", "yellow", "red"]
    vector_color_index = 0
    vector_count = len(vectors_parameters[0]["sum_of_ra_diffs"])
    for spt_index in range(0, len(spots_parameters)):
        spt = spots_parameters[spt_index]
        spot = Circle(spt["coords"], angle=0, lw=0.5, radius=spt["radius"])
        ax1.add_artist(spot)
        spot_color = cm.jet((spt["vel"] - min(velocity)) / velocity_range, 1)
        spot.set_facecolor(spot_color)

        vect = vectors_parameters[spt_index]
        vect_ra = vect["sum_of_ra_diffs"]
        vect_dec = vect["sum_of_dec_diffs"]

        for vec in range(vector_count):
            plt.annotate("", xy=spt["coords"], xycoords='data',
                         xytext=(spt["coords"][0] + (20 * vect_ra[vec]), spt["coords"][1] + (20 * vect_dec[vec])),
                         textcoords='data',
                         arrowprops=dict(arrowstyle="<-", color=vector_colors[vector_color_index], connectionstyle="arc3"))
            vector_color_index += 1

            if vector_color_index == len(vector_colors):
                vector_color_index = 0

    # vector legend
    plt.annotate("", xy=(50, -150), xycoords='data', xytext=(50 + (20 * 3), -150), textcoords='data',
                  arrowprops=dict( arrowstyle="<-", connectionstyle="arc3"))
    plt.text(105, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='k')
    plt.text(110, -170, "5.8 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='k')

    plt.annotate("", xy=(-60, -150), xycoords='data', xytext=(-60 + (20 * 3), -150), textcoords='data',
                  arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3"))
    plt.text(-5, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='grey')
    plt.text(0, -170, "7.2 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='grey')

    plt.text(100, 220, "31/Oct/2019", size=12, rotation=0.0, ha="left", va="center", color='k')
    plt.text(100, 250, "31/Oct/2011", size=12, rotation=0.0, ha="left", va="center", color='grey')
    plt.text(100, 280, "11/Mar/2009", size=12, rotation=0.0, ha="left", va="center", color='k')
    patches = []
    colors = velocity
    p = PatchCollection(patches, cmap=cm.jet)
    p.set_array(np.array(colors))
    ax1.add_collection(p)
    plt.colorbar(p)
    all_ra = [ellipse["coords"][0] for ellipse in spots_parameters]
    all_dec = [ellipse["coords"][1] for ellipse in spots_parameters]
    ax1.set_xlim(max(all_ra) + 50, min(all_ra) - 50)
    ax1.set_ylim(min(all_dec) - 50, max(all_dec) + 50)
    plt.xlabel('$\\Delta$ RA [mas]', fontsize=12)
    plt.ylabel('$\\Delta$ Dec [mas]', fontsize=12)
    plt.title("G78: SW excluded", size=12)
    plt.grid(True)
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
