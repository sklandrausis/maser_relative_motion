import sys
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

from parsers.configparser_ import ConfigParser


def check_if_group_is_in_file(file, group):
    input_file = file
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


def get_configs_items():
    """
    :return: None
    """
    config_file_path = "config/plot.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items("main")


def compute_distance(ra1, ra2, dec1, dec2):
    return np.sqrt((ra1-ra2)**2 + (dec1-dec2)**2)


def main():
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    groups_file_path = get_configs("paths", "groups")
    warnings.filterwarnings("ignore")

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    groups = [1, 2, 3, 4, 5, 6, 7]
    dtype = [('group_nr', int), ('velocity', float), ('intensity', np.float128), ("ra", float), ("dec", float)]

    distances = dict()
    epochs = []
    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}
    dates_ = []

    for index in range(0, len(input_files)):
        epoch = input_files[index].split(".")[0]
        epochs.append(epoch)
        dates_.append(dates[epoch])
        input_file = groups_file_path + epoch + ".groups"
        group_tmp, velocity_tmp, intensity_tmp, ra_tmp, dec_tmp = \
            np.loadtxt(input_file, unpack=True, usecols=(0, 2, 3, 5, 6))
        values = [(group_tmp[ch], velocity_tmp[ch], intensity_tmp[ch], ra_tmp[ch], dec_tmp[ch])
                  for ch in range(0, len(group_tmp))]
        data = np.array(values, dtype=dtype)
        data = np.sort(data, order=['group_nr', 'velocity'])

        data_for_groups = {str(group_number): data[data["group_nr"] == group_number] for group_number in groups}

        del_groups = []
        for group in data_for_groups:
            if len(data_for_groups[group]) == 0:
                del_groups.append(group)

        for del_group_index in del_groups:
            del data_for_groups[del_group_index]

        max_intensity_for_group = []
        for group in data_for_groups:
            max_intensity = max(data_for_groups[group]["intensity"])
            reference_index = np.where(data_for_groups[group]["intensity"] == max_intensity)[0][0]
            max_intensity_for_group.append(data_for_groups[group][reference_index])

        group_index = [g[0] for g in max_intensity_for_group]
        distances_index = []

        for g in group_index:
            for g_ in group_index:
                if str(g) + "_" + str(g_) and str(g_) + "_" + str(g) not in distances_index:
                    if str(g) != str(g_):
                        distances_index.append(str(g) + "_" + str(g_))
                        if str(g) + "_" + str(g_) not in distances.keys():
                            distances[str(g) + "_" + str(g_)] = []

        group_distances = {key: None for key in distances_index}

        for g in group_distances.keys():
            g_ = g.split("_")

            g1 = max_intensity_for_group[int(g_[0])-1]
            g2 = max_intensity_for_group[int(g_[1])-1]
            ra1 = g1[3]
            ra2 = g2[3]
            dec1 = g1[4]
            dec2 = g2[4]
            distances[g].append({epoch:compute_distance(ra1, ra2, dec1, dec2)})

    for g in distances.keys():
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=100)

        for e in distances[g]:
            ax.scatter(dates[list(e.keys())[0]], e.values(), s=100)

        ax.set_title(g)

    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
