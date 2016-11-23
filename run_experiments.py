import argparse
import sys
import os
import subprocess


def execute_one_experiment(executable, otn_topology_file, ip_topology_file,
                           ip_node_mapping_file, ip_link_mapping_file,
                           ip_port_info_file, vn_topology_file,
                           location_constraint_file):
    process = subprocess.Popen([executable, "--otn_topology_file=" + otn_topology_file,
            "--ip_topology_file=" + ip_topology_file, 
            "--ip_node_mapping_file=" + ip_node_mapping_file,
            "--ip_link_mapping_file=" + ip_link_mapping_file,
            "--ip_port_info_file=" + ip_port_info_file,
            "--vn_topology_file=" + vn_topology_file,
            "--vn_location_file=" + location_constraint_file],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    with open(vn_topology_file + ".status") as f:
        print vn_topology_file + ": " + f.readline()

def main():
    parser = argparse.ArgumentParser(
        description="Script for automating Multilayer VNE experiments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--testcase_root',
        help='Root directory for test cases',
        required=True)
    parser.add_argument(
        '--executable',
        help='Name of the executable file to run',
        required=True)
    args = parser.parse_args()
    root = args.testcase_root
    executable = './' + args.executable
    subprocess.Popen(['make'], shell=False, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    otn_topology_file = args.testcase_root + "/otn.topo"
    ip_topology_file = args.testcase_root + "/ip.topo"
    ip_node_mapping_file = args.testcase_root + "/ip.nmap"
    ip_link_mapping_file = args.testcase_root + "/ip.emap"
    ip_port_info_file = args.testcase_root + "/ip-port"
    i = 0
    while True:
        vn_topology_file = os.path.join(args.testcase_root, "vn" + str(i) + ".topo")
        location_constraint_file = os.path.join(args.testcase_root, "vn" + str(i) + "loc")
        if not os.path.isfile(vn_topology_file):
            break
        execute_one_experiment(executable, otn_topology_file, ip_topology_file,
                               ip_node_mapping_file, ip_link_mapping_file,
                               ip_port_info_file, vn_topology_file,
                               location_constraint_file)
        i = i + 1
if __name__ == "__main__":
    main()
