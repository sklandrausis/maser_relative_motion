import sys
import argparse
import random
from multiprocessing import Pool
import numpy as np

from parsers.configparser_ import ConfigParser

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


def create_groups(test):
    print(len(test))
    channel, velocity, intensity, integral_intensity, ra, dec = test

    groups = []
    param1 = random.uniform(0.1, 1)
    param2 = random.uniform(1, 100)

    #print(param1, param2)

    all_chans = []
    for j in range(0, len(channel)):
        chan = [channel[j], velocity[j], intensity[j], integral_intensity[j], ra[j], dec[j]]
        group = []
        for i in range(0, len(channel)):
            if channel[i] not in get_all_channels_from_group(group):
                if np.abs(velocity[i] - chan[1]) <= param1 and \
                        np.abs(distance(ra[i], chan[4], dec[i], chan[5])) <= param2:
                    if [channel[i], velocity[i], intensity[i], integral_intensity[i], ra[i], dec[i]] not in all_chans:
                        all_chans.append([channel[i], velocity[i], intensity[i], integral_intensity[i], ra[i], dec[i]])
                        group.append([channel[i], velocity[i], intensity[i], integral_intensity[i], ra[i], dec[i]])

        if len(group) >= 1:
            groups.append(group)

    return groups


def main(input_files_dir):
    p = Pool(7)
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]

    data = []
    for file in file_order:
        input_file = input_files_dir + "/" + file
        channel, velocity, intensity, integral_intensity, ra, dec = np.loadtxt(input_file, unpack=True)
        velocity = velocity / 1000
        data.append([channel, velocity, intensity, integral_intensity, ra, dec])

    result = dict()
    for i in range(0, 1):
        print("i", i)
        result[str(i)] = []
        groups = p.map(create_groups, data)
            #result[str(i)].append(len(groups))

    '''
    total_results = dict()
    keys = result.keys()
    for key in keys:
        total_results[key] = sum(result[key])

    print("Min numbers of groups", min(total_results.values()))
    tmp = list(total_results.values()).index(min(total_results.values()))
    print("params", list(total_results.keys())[tmp])
    print(result)
    '''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot gauss')
    parser.add_argument('input_file_dir', type=str, help='Input files directory')
    args = parser.parse_args()
    main(args.input_file_dir)
    sys.exit(0)
