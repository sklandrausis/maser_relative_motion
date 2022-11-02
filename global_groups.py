import sys
import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file key
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def compute_weight_of_edge(node1, node2, data):
    ra1 = data[node1]["ra"]
    ra2 = data[node2]["ra"]
    dec1 = data[node1]["dec"]
    dec2 = data[node2]["dec"]
    return np.sqrt((ra1 - ra2) ** 2 + (dec1 - dec2) ** 2)


def compere_velocities(node1, node2, data):
    vel1 = data[node1]["velocity"]
    vel2 = data[node2]["velocity"]
    a1 = vel1 - 0.02
    b1 = vel1 + 0.02
    a2 = vel2 - 0.02
    b2 = vel2 + 0.02

    if vel1 == vel2:
        return True
    else:
        if a1 <= vel2 <= b1:
            return True

        if a2 <= vel1 <= b2:
            return True

    return False


def get_data():
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    data_file_path = get_configs("paths", "dataFiles")
    data_files = os.listdir(data_file_path)
    order = [data_files.index(ind) for ind in file_order]
    data_file_tmp = []

    for ind in order:
        data_file_tmp.append(data_files[ind])
    data_files = data_file_tmp
    data = []

    for file in data_files:
        data_tmp = np.loadtxt(data_file_path + file, usecols=(0, 1, 2, 3, 4, 5), unpack=True, dtype=np.float)
        rows = data_tmp.shape[1]
        for r in range(rows):
            row = data_tmp[:, r]
            data.append({"ch": row[0], "velocity": row[1], "flux1": row[2], "flux2": row[3], "ra": row[4],
                         "dec": row[5], "file": file})

    return data, len(data_files)


def print_group(group, graph):
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    velocities = [nx.get_node_attributes(graph, 'velocity')[g] for g in group]
    files = [nx.get_node_attributes(graph, 'file')[g] for g in group]

    if len(set(files)) == len(file_order):
        order = [files.index(ind) for ind in file_order]

        velocities_tmp = []
        files_tmp = []

        for ind in order:
            velocities_tmp.append(velocities[ind])
            files_tmp.append(files[ind])

        velocities = velocities_tmp
        files = files_tmp

        print(velocities, files)


def create_output(groups, graph):
    index = 0
    header = ['vel']
    data = []
    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]

    for group in groups:
        files = [nx.get_node_attributes(graph, 'file')[g] for g in group]
        ras = [nx.get_node_attributes(graph, 'ra')[g] for g in group]
        decs = [nx.get_node_attributes(graph, 'dec')[g] for g in group]
        flux1s = [nx.get_node_attributes(graph, 'flux1')[g] for g in group]
        velocities = [nx.get_node_attributes(graph, 'velocity')[g] for g in group]

        if len(set(files)) == len(file_order):
            order = [files.index(ind) for ind in file_order]

            ras_tmp = []
            dec_tmp = []
            files_tmp = []
            flux1s_tmp = []
            velocities_tmp = []

            for ind in order:
                files_tmp.append(files[ind])
                ras_tmp.append(ras[ind])
                dec_tmp.append(decs[ind])
                flux1s_tmp.append(flux1s[ind])
                velocities_tmp.append(velocities[ind])

            index += 1
            velocities = velocities_tmp
            files = files_tmp
            ras = ras_tmp
            decs = dec_tmp
            flux1s = flux1s_tmp

            if index == 1:
                for file in files:
                    header.append("ra" + "_" + file)
                    header.append("dec" + "_" + file)
                    header.append("flux1" + "_" + file)

            data_tmp = [velocities[-1]]
            for file_index in range(0, len(files)):
                data_tmp.extend([ras[file_index], decs[file_index], flux1s[file_index]])

            data.append(data_tmp)

    relative_motion_path = get_configs("paths", "relative_motion")
    np.savetxt(relative_motion_path + 'relative_motion_group.dat', np.array(data), delimiter=",",
               header=",".join(header))


def main():
    data, file_count, = get_data()
    graph = nx.Graph()
    number_of_points = len(data)
    nodes = [node for node in range(number_of_points)]
    pos = {node: (data[node]["ra"], data[node]["dec"]) for node in nodes}
    labels = {node: str(data[node]["velocity"]) + "_" + data[node]["file"] for node in nodes}
    radius = [data[node]["flux1"] * 25 for node in nodes]

    for node in nodes:
        graph.add_node(node, velocity=data[node]["velocity"], label=data[node]["velocity"],
                       file=data[node]["file"], ra=data[node]["ra"],
                       dec=data[node]["dec"], flux1=data[node]["flux1"])

    edges = []
    for i in range(0, len(nodes)):
        for j in range(0, len(nodes)):
            if i != j and compute_weight_of_edge(nodes[i], nodes[j], data) < 10.0 and \
                    compere_velocities(nodes[i], nodes[j], data) and data[nodes[i]]["file"] != data[nodes[j]]["file"]:
                edges.append((nodes[i], nodes[j], compute_weight_of_edge(nodes[i], nodes[j], data)))

    graph.add_weighted_edges_from(edges)

    groups = list(nx.connected_components(graph))
    number_of_connected_components = nx.number_connected_components(graph)
    total_group_count = 0
    group_count_that_have_all_files = 0
    full_groups = []
    for group in groups:
        if len(group) > 1:
            total_group_count += 1
        if len(group) == file_count:
            group_count_that_have_all_files += 1
            full_groups.append(group)

    create_output(groups, graph)

    print("Total Group count is ", total_group_count)
    print("Group count that have all files is ", group_count_that_have_all_files)
    print("Single maser count ", number_of_connected_components - total_group_count)

    fig, ax = plt.subplots()
    nx.draw(graph, with_labels=False, pos=pos, cmap="jet", node_color=[vel["velocity"] for vel in data], ax=ax,
            node_size=radius)
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=6, ax=ax)
    plt.axis('on')
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    vmax = min([vel["velocity"] for vel in data])
    vmin = max([vel["velocity"] for vel in data])
    sm = plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    ax.grid(True)
    plt.colorbar(sm, shrink=0.5, ax=ax)
    fig.tight_layout()
    plt.tight_layout()
    plt.ylabel("DEC")
    plt.xlabel("RA")
    plt.subplots_adjust(top=0.986, bottom=0.053, left=0.037, right=0.992, hspace=0.3, wspace=0.3)
    plt.show()
    plt.draw()

    sys.exit()


if __name__ == "__main__":
    main()
