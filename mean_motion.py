import sys
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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

    ls = []
    lstex = []
    linearity = []

    plt.figure()
    ax1 = plt.subplot(111, aspect='equal')
    #print(ras[2], ras[0],  ras[2] - ras[0])
    for group in groups_indexies:
        sum_of_ra_diff = []
        sum_of_dec_diff = []
        sum_of_ra = []
        sum_of_dec = []
        sum_of_vel = []
        length = []
        ls_tmp = []

        for index in range(0, len(ras)):
            #print(group)
            #print(ras[index])
            if index != 0:
                ra_diff = [ras[index][gi] - ras[0][gi] for gi in group]
                dec_diff = [decs[index][gi] - decs[0][gi] for gi in group]
                sum_ra_diff = sum(ra_diff)
                sum_dec_diff = sum(dec_diff)

            ra_for_group = [ras[index][gi] for gi in group]
            sum_of_ra = sum(ra_for_group)
            dec_for_group = [decs[index][gi] for gi in group]
            sum_of_dec = sum(dec_for_group)
            print(sum_of_ra)

            #print(sum(ra_diff), len(ra_diff), ra_diff)
        print("\n\n")











        '''


            

            dec_for_group = [decs[index][gi] for gi in group]
            sum_of_dec.append(sum(dec_for_group) / len(group))

            ls_tmp.append(sum(ra_for_group) / len(group))
            ls_tmp.append(sum(dec_for_group) / len(group))

            flux_for_group = [fluxs[index][gi] for gi in group]
            ls_tmp.append(max(flux_for_group))



            sum_of_ra_diff.append( sum( ra_diff_tmp ) )


            sum_of_dec_diff.append( sum( dec_diff_tmp ) )

            length.append( np.sqrt( np.array( sum_of_ra_diff ) ** 2 + np.array( sum_of_dec_diff ) ** 2 ) )

            ls_tmp.append( sum( ra_diff_tmp ) / len( group ) )
            ls_tmp.append( sum( dec_diff_tmp ) / len( group ) )

        ls_tmp.append(length)
        ls_tmp.append(np.array(length, dtype=object) /len(group))

        for i in range(0, len(sum_of_ra_diff)):
            ax1.annotate("", xy=(sum_of_ra[i] / len( group ), sum_of_dec[i] / len( group )), xycoords='data',
                          xytext=(sum_of_ra[i] / len( group ) + (20 * sum_of_ra_diff[i] / len( group )),
                                  sum_of_dec[i] / len( group ) + (20 * sum_of_dec_diff[i] / len( group ))),
                          textcoords='data', arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3"))

        #plt.annotate( "", xy=(wek_x13 / k13, wek_y13 / k13), xycoords='data',
                  #xytext=(wek_x13 / k13 + (20 * xx13 / k13), wek_y13 / k13 + (20 * yy13 / k13)),
                  #textcoords='data', arrowprops=dict(arrowstyle="<-", connectionstyle="arc3" ) )

        ls.append(ls_tmp)
    for group_index in range(0, len(groups_indexies)):
        #print(ls[group_index][9])
        ellipse = Ellipse((ls[group_index][1], ls[group_index][2]), angle=0, lw=0.5, width=3*np.log10(ls[group_index][9]*1000.), height=3*np.log10(ls[group_index][9]*1000.))
        ax1.add_artist(ellipse)
        collor = cm.jet((ls[group_index][0] - min(velocity)) / velocity_range, 1)
        ellipse.set_facecolor(collor)

    # vector legend
    plt.annotate("", xy=(50, -150), xycoords='data', xytext=(50 + (20 * 3), -150), textcoords='data',
              arrowprops=dict( arrowstyle="<-", connectionstyle="arc3"))
    plt.text(105, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='k')
    plt.text(110, -170, "5.8 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='k')

    plt.annotate( "", xy=(-60, -150), xycoords='data', xytext=(-60 + (20 * 3), -150), textcoords='data',
              arrowprops=dict( arrowstyle="<-", color="grey", connectionstyle="arc3"))
    plt.text(-5, -135, "3 mas", size=8, rotation=0.0, ha="left", va="center", color='grey')
    plt.text(0, -170, "7.2 km s$^{-1}$", size=8, rotation=0.0, ha="left", va="center", color='grey')

    patches = []
    colors = velocity
    p = PatchCollection( patches, cmap=cm.jet)
    p.set_array(np.array(colors))
    ax1.add_collection(p)
    plt.colorbar(p)
    plt.text(100, 220, "31/Oct/2019", size=12, rotation=0.0, ha="left", va="center", color='k')
    plt.text(100, 250, "31/Oct/2011", size=12, rotation=0.0, ha="left", va="center", color='grey')
    plt.text(100, 280, "11/Mar/2009", size=12, rotation=0.0, ha="left", va="center", color='k')
    ax1.set_xlim(150, -350)
    ax1.set_ylim(-200, 300)
    plt.xlabel('$\\Delta$ RA [mas]', fontsize=12)
    plt.ylabel('$\\Delta$ Dec [mas]', fontsize=12)
    plt.title("G78: SW excluded", size=12)
    plt.grid(True)
    plt.show()
    '''
    sys.exit(0)


if __name__ == "__main__":
    main()
