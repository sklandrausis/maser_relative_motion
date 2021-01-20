import sys
import argparse

import numpy as np

from parsers.configparser_ import ConfigParser


def check_if_group_is_in_file(file, group):
    input_file = "groups/" + "/" + file
    group_nr = np.loadtxt(input_file, unpack=True, usecols=0)

    if group not in group_nr:
        return False
    else:
        return True


def get_configs(section, key):
    """
    :param section: configuration file secti
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main(group_number):
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    bad_files = []
    for file in input_files:
        if not check_if_group_is_in_file(file.split(".")[0] + ".groups", group_number):
            bad_files.append(file)
            del dates[file.split(".")[0]]
    input_files = [file for file in input_files if file not in bad_files]

    data = dict()
    for index in range(0, len(input_files)):
        data[input_files[index].split(".")[0]] = dict()
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        intensity = np.empty(0)
        channels = np.empty(0)
        ra = np.empty(0)
        dec = np.empty(0)
        group_tmp, channel_tmp, intensity_tmp, ra_tmp, dec_tmp = \
            np.loadtxt(input_file, unpack=True, usecols=(0, 1, 3, 5, 6))
        for i in range(0, len(channel_tmp)):
            if group_tmp[i] == int(group_number):
                intensity = np.append(intensity, intensity_tmp[i])
                channels = np.append(channels, channel_tmp[i])
                ra = np.append(ra, ra_tmp[i])
                dec = np.append(dec, dec_tmp[i])

        index_for_new_zero = np.where(intensity == max(intensity))[0][0]
        data[input_files[index].split(".")[0]]["max_intensity"] = max(intensity)
        data[input_files[index].split(".")[0]]["ch_for_max_intensity"] = channels[index_for_new_zero]
        data[input_files[index].split(".")[0]]["ra"] = ra - ra[index_for_new_zero]
        data[input_files[index].split(".")[0]]["dec"] = dec - dec[index_for_new_zero]

    print(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    args = parser.parse_args()
    main(args.group_number)
    sys.exit(0)
