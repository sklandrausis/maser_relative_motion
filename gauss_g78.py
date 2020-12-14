import sys
import argparse
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


def distance(ra1, ra2, dec1, dec2):
    return np.sqrt((ra1-ra2)**2 + (dec1-dec2)**2)


def get_all_channels_from_group(group):
    channels = []
    for ch in group:
        channels.append(ch[0])
    return channels


def get_velocity_from_group(group):
    velocities = []
    for ch in group:
        velocities.append(ch[1])
    return velocities


def get_configs(section, key):
    """

    :param section: configuration file secti
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main(input_files_dir):
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    minorLocatorvel = MultipleLocator(0.5)
    fig, ax = plt.subplots(nrows=1, ncols=len(file_order), figsize=(16, 16))

    for file in file_order:
        index = file_order.index(file)
        title = file.split(".")[0].upper() + "-" + dates[file.split(".")[0]]
        input_file = input_files_dir + "/" + file
        channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt(input_file, unpack=True)
        velocity = velocity/1000

        groups = []
        group = [[channel[0], velocity[0], intensity[0], integral_intensity[0], ra[0], dec[0]]]
        groups.append(group)


        '''
        for j in range(1, len(channel)):
            chan = [channel[j], velocity[j], intensity[j], integral_intensity[j], ra[j], dec[j]]
            for group in groups:
                tmp = group[-1]
                vel_tmp = chan[1]
                distance_tmp = np.abs(distance(chan[4], tmp[4], chan[5], tmp[5]))
                velocities_for_group = get_velocity_from_group(group)
                max_velocities_for_group = max(velocities_for_group)
                min_velocities_for_group = min(velocities_for_group)
                if distance_tmp >= 12 and np.abs(tmp[1] - vel_tmp) <= 0.65:
                    group.append(chan)
                    break
                else:
                    new_group = [chan]
                    groups.append(new_group)
                    break

        '''

        all_chans = []
        for j in range(0, len(channel)):
            chan = [channel[j], velocity[j], intensity[j], integral_intensity[j], ra[j], dec[j]]
            group = []
            for i in range(0, len(channel)):
                if channel[i] not in get_all_channels_from_group(group):
                    if np.abs(velocity[i] - chan[1]) <= 0.9727272727272727 and \
                            np.abs(distance(ra[i], chan[4], dec[i], chan[5])) <= 9.72072072072072:
                        if [channel[i], velocity[i], intensity[i], integral_intensity[i], ra[i], dec[i]] not in all_chans:
                            all_chans.append([channel[i], velocity[i], intensity[i], integral_intensity[i], ra[i], dec[i]])
                            group.append([channel[i], velocity[i], intensity[i], integral_intensity[i], ra[i], dec[i]])

            if len(group) >= 1:
                groups.append(group)
        #'''

        '''
        test = []
        for g in groups:
            test.append(get_all_channels_from_group(g))
    
        test = [item for sublist in test for item in sublist]
        test.sort()
    
        for a in range(0, len(channel) - len(test)):
            test.append(0)
    
        tmp = list(channel - np.array(test)).index(-1)
        print(tmp)
        sys.exit(0)
       
        '''

        print("Numbers of groups", len(groups), "Number of grouped chans", len(all_chans), "Total channals ", len(channel))

        for group in groups:
            vel = []
            inten = []
            group_len = len(group)
            for chan in group:
                vel.append(chan[1])
                inten.append(chan[2])

            p0 = [max(inten), min(vel) + 0.5 * (max(vel) - min(vel)), 0.2]

            try:
                if group_len >= 3:
                    coeff, var_matrix = curve_fit(gauss, vel, inten, p0=p0)
                    q = np.linspace(min(vel) - 0.2, max(vel) + 0.2, 1000)
                    hist_fit = gauss(q, *coeff)
                    ax[index].plot(q, hist_fit, 'k')

                size = []
                for j in range(0, len(vel)):
                    for k in range(j + 1, len(vel)):
                        dist = np.sqrt(((ra[j] - ra[k]) * 11281) ** 2 + ((dec[j] - dec[k]) * 1000) ** 2)
                        size.append(dist)

                ax[index].plot(vel, inten, ls="", marker="o", markersize=5, markerfacecolor="white", markeredgewidth=1, label=str(groups.index(group)))

            except:
                pass

        ax[index].set_xlabel('$V_{\\rm LSR}$ [km s$^{-1}$]', fontsize=12)
        ax[index].set_title(title, size=12)
        ax[index].xaxis.set_minor_locator(minorLocatorvel)
    ax[0].set_ylabel('Flux density [Jy]', fontsize=12)
    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0.08, wspace=0.18, hspace=0, left=0.05, right=0.99)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot gauss')
    parser.add_argument('input_file_dir', type=str, help='Input files directory')
    args = parser.parse_args()
    main(args.input_file_dir)
    sys.exit(0)
