# Given a text file with format:
#   [Tue Apr  8 02:24:59 UTC 2025]
# * CPU usage: 0.0%
# * Memory usage: 1.97 GiB 1.6%
# * Disk usage: 15.000 GiB 57%
# * Read/Write IO: 0.000 MiB/s 0.002 MiB/s
# Extract only memory usage rows, collect their valuev before "GiB", and plot it as a line plot.
import matplotlib.pyplot as plt
import re
import os

def extract_memory_usage(file_path):
    memory_usages = []
    with open(file_path, 'r') as file:
        for line in file:
            if "Memory usage" in line:
                match = re.search(r'(\d+\.\d+) GiB', line)
                if match:
                    memory_usages.append(float(match.group(1)))
    return memory_usages

def plot_memory_usage(memory_usages_by_file):
    plt.figure(figsize=(10, 5))
    for file_path, memory_usages in memory_usages_by_file.items():
        label = os.path.basename(file_path)
        plt.plot(memory_usages, marker='', label=label)

    plt.title('RAM usage over time')
    plt.xlabel('Time (every 10 seconds)')
    plt.ylabel('RAM (GiB)')
    plt.grid()
    plt.legend()
    plt.show()

def plot_max_memory_histogram(max_memory_usages):
    plt.figure(figsize=(8, 5))
    plt.hist(max_memory_usages, bins='auto', edgecolor='black')
    plt.title('Histogram of max RAM usage per log file')
    plt.xlabel('Max RAM (GiB)')
    plt.ylabel('Number of files')
    plt.grid(axis='y')
    plt.show()

if __name__ == "__main__":
    input_dir = '/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/ram_logs/n0512'
    #input_dir = '/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/ram_logs/n1024'
    file_paths = sorted(
        os.path.join(input_dir, name)
        for name in os.listdir(input_dir)
        if name.endswith('.log')
    )
    if not file_paths:
        raise ValueError(f'No .log files found in directory: {input_dir}')

    memory_usages_by_file = {
        file_path: extract_memory_usage(file_path)
        for file_path in file_paths
    }

    max_memory_usages = [
        max(memory_usages)
        for memory_usages in memory_usages_by_file.values()
        if memory_usages
    ]
    if not max_memory_usages:
        raise ValueError('No memory usage values found in the discovered .log files')

    plot_memory_usage(memory_usages_by_file)
    plot_max_memory_histogram(max_memory_usages)
    
