# Given text files with the following format:
#
# name,h,m,s,ms,success,notes
# variant_graph,0,0,0,0,1,NA
# optimize_d_plus_n,0,0,0,33,1,Optimal;n_reads=152;n_paths=43;n_edges=470
#
# isolate only the optimize_d_plus_n rows, compute the time in seconds, extract n_edges, and plot time vs n_edges as a line plot.
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import re
import os

def extract_optimize_d_plus_n_data(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("optimize_d_plus_n"):
                parts = line.strip().split(',')
                h, m, s, ms = map(int, parts[1:5])
                time_in_seconds = h * 3600 + m * 60 + s + ms / 1000
                match = re.search(r'n_edges=(\d+)', parts[6])
                if match:
                    return time_in_seconds, int(match.group(1))
    return None

def collect_optimize_d_plus_n_data(input_dir):
    times = []
    n_edges = []
    log_files = []

    for entry in sorted(os.scandir(input_dir), key=lambda item: item.name):
        if not entry.is_dir():
            continue

        log_file = os.path.join(entry.path, 'log.csv')
        if not os.path.isfile(log_file):
            continue

        log_files.append(log_file)
        record = extract_optimize_d_plus_n_data(log_file)
        if record is None:
            continue

        time_in_seconds, edge_count = record
        times.append(time_in_seconds)
        n_edges.append(edge_count)

    return times, n_edges, log_files

def plot_time_vs_n_edges_heatmap(times, n_edges):
    max_n_edges = 16000
    max_time = 1

    # Filter data to only include points within specified ranges
    filtered_times = []
    filtered_n_edges = []
    for t, e in zip(times, n_edges):
        if 0 <= e <= max_n_edges and 0 <= t <= max_time:
            filtered_times.append(t)
            filtered_n_edges.append(e)
    
    plt.figure(figsize=(10, 10))
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_under('white')
    _, _, _, image = plt.hist2d(
        filtered_n_edges,
        filtered_times,
        bins=(250, 250),
        cmap=cmap,
        cmin=1,
        norm=LogNorm(),
    )

    plt.title('Runtime of optimize_d_plus_n on 1024 samples')
    plt.xlabel('Number of edges')
    plt.ylabel('Seconds')
    plt.xlim(0, max_n_edges)
    plt.ylim(0, max_time)
    plt.grid()
    plt.colorbar(image, label='Number of windows')
    plt.show()

if __name__ == "__main__":
    input_dir = '/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/runtime_logs/n1024/tmp.w2HIMS/tmp.d2QV32P002/output/run'
    times, n_edges, log_files = collect_optimize_d_plus_n_data(input_dir)

    if not log_files:
        raise ValueError(f'No log.csv files found in subdirectories of: {input_dir}')
    if not times:
        raise ValueError('No optimize_d_plus_n records found in discovered log.csv files')

    plot_time_vs_n_edges_heatmap(times, n_edges)
    
