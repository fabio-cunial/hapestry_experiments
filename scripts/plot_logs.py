from cmath import log
from collections import defaultdict
from math import log10
import math
import os
import matplotlib
from matplotlib import pyplot
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, LogLocator
import numpy as np


'''
Remark: when `graphaligner` fails, `window_total` is marked as failed, but 
when `optimize_d_plus_n` fails, `window_total` is not marked as failed.
Examples:

name,h,m,s,ms,success,notes
variant_graph,0,0,0,0,1,NA
graphaligner,0,6,40,7,0,NA
window_total,0,6,40,8,0,NA

name,h,m,s,ms,success,notes
variant_graph,0,0,0,0,1,NA
graphaligner,0,0,10,772,1,NA
align_reads_to_paths,0,1,4,310,1,NA
rescale_weights,0,0,17,891,1,NA
feasibility_construct,0,0,11,124,1,NA
feasibility_init,0,0,20,457,1,NA
feasibility,0,0,58,93,1,Optimal;n_read_hap_vars=2803944;r_in=8151;r_out=8151
compress_transmap_initial,0,0,3,141,1,NA
compress_transmap_d,0,0,14,527,1,NA
optimize_d_prune_construct,0,0,1,58,1,NA
optimize_d_prune_init,0,0,1,712,1,NA
optimize_d_prune,0,0,5,979,1,Optimal;n_reads=7051;n_paths=344;n_edges=336661
optimize_d_prune_parse,0,0,0,459,1,NA
compress_transmap_d_plus_n,0,0,7,340,1,NA
optimize_d_plus_n_construct,0,0,3,381,1,NA
optimize_d_plus_n_init,0,0,6,47,1,NA
optimize_d_plus_n,0,15,0,895,0,Feasible;n_reads=6933;n_paths=148;n_edges=985134
optimize_d_plus_n_parse,0,0,3,801,1,NA
window_total,0,19,11,312,1,NA
'''


def aggregate_logs(parent_dirs, times, common_windows, name):
    '''
    Remark: `times[name][task]` contains every window with `task` in its log,
    even if that window failed on that task. `common_windows` contains only 
    windows that have every task in their log.

    Remark: `times` are in minutes.
    '''
    success_counter = [set(),set()]
    fail_reasons = defaultdict(int)
    windows = defaultdict(set)
    for directory in parent_dirs:
        for root, dirs, files in sorted(os.walk(directory)):
            for file in files:
                file_path = os.path.join(root, file)
                if "log.csv" in file_path:
                    window = root.strip('/').split("/")[-1]
                    with open(file_path, 'r') as file:
                        for l,line in enumerate(file):
                            if l == 0:
                                continue
                            task, h, m, s, ms, success, notes = line.strip().split(',')
                            windows[task].add(window)
                            if not bool(int(success)):
                                if task != "window_total":
                                    fail_reasons[task] += 1
                                    print(success, window, task, notes)
                                # continue
                            success_counter[int(success)].add(window)
                            s_total = float(h)*60*60 + float(m)*60 + float(s) + float(ms)/1000.0
                            m_total = s_total / 60.0
                            if m_total >= 100:
                                print("Window runtime >=100m ?! ", window, task, m_total)
                            times[name][task].append((window,m_total))
    print(name, len(success_counter[1]), len(success_counter[0]), ','.join(["%s:%d"%(k,v) for k,v in fail_reasons.items()]))
    for task in windows:
        if len(common_windows[task]) > 0:
            common_windows[task] = common_windows[task].intersection(windows[task])
            # just for convenience of debugging or general edification
            # uncommon_windows = windows[task].difference(common_windows[task])
            # print(task, uncommon_windows)
        else:
            common_windows[task] = windows[task]
    return times, common_windows




def plot_scaling(times, task_names, log_scale=True, ax=None, show=True):
    '''
    The runtime of the few slowest tasks as a function of the number of samples, 
    represented as a log2-log2 line plot.
    '''
    x_ticks = list()
    colormap = matplotlib.colormaps["nipy_spectral"]
    if ax is None:
        fig, ax = pyplot.subplots()
    x_per_task = defaultdict(list)
    y_per_task = defaultdict(list)
    polynomial_descriptions = list()
    x = list()
    ordering = lambda x: int(x[1:])  # Assume runs are named n0001,...,n1024.
    i=0
    for r,run_name in enumerate(sorted(times, key=ordering)):
        i+=1
        #x.append(i)
        x.append(ordering(run_name))
        #x_labels.append(ordering(run_name))
        for t,task_name in enumerate(task_names):
            if task_name not in times[run_name]:
                continue
            #y = sum([t for w,t in times[run_name][task_name] if w in common_windows[task_name]])
            x_per_task[task_name].append(ordering(run_name))
            y = sum([t for w,t in times[run_name][task_name]])
            y_per_task[task_name].append(y)
    for i,(task_name,y) in enumerate(y_per_task.items()):
        x_task = x_per_task[task_name]
        color_index = (float(i) + 1) /float(len(task_names))
        color = colormap(color_index)
        if "total" in task_name:
            color = "gray"
        ax.plot(x_task, y, color=color, label=task_name, marker='o')
        x_ticks.append(x)

        # Compute the interpolating polynomial in log2-log2 space so it is
        # compatible with the optional log-scaled axes and stays positive.
        valid = [(xv, yv) for xv, yv in zip(x_task, y) if xv > 0 and yv > 0 and xv >= 64]  # Only fit on the larger n_samples, which are more indicative of the asymptotic behavior.
        if len(valid) >= 2:
            x_fit = np.array([p[0] for p in valid], dtype=float)
            y_fit = np.array([p[1] for p in valid], dtype=float)
            degree = 1
            coeffs = np.polyfit(np.log2(x_fit), np.log2(y_fit), degree)
            polynomial = np.poly1d(coeffs)

            x_eval = np.logspace(6, 10, base=2.0)
            y_eval = 2**polynomial(np.log2(x_eval))
            ax.plot(x_eval, y_eval, color=color, linestyle='--', linewidth=1.0, alpha=0.8)

            coeff_text = np.array2string(coeffs, precision=3, separator=',')
            polynomial_descriptions.append(f"{task_name}: coeffs={coeff_text}")
    ax.set_xticks(x)
    ax.set_xlabel("Number of samples")
    if log_scale:
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.set_ylabel("Minutes")
    else:
        ax.set_ylabel("Minutes")
    ax.set_box_aspect(1)
    # ax.set_adjustable("box")
    # ax.xaxis.set_major_locator(LogLocator(base=2))
    # ax.yaxis.set_major_locator(LogLocator(base=2, subs=(1.0,), numticks=1000))
    # ax.xaxis.set_major_formatter(LogFormatterMathtext(base=2))
    # ax.yaxis.set_major_formatter(LogFormatterMathtext(base=2))
    # ax.set_aspect(1, adjustable="box")
    # ax.grid(True, axis="both", which="both", color='lightgray', linestyle='-', linewidth=0.5)
    # ax.legend(loc='lower left', bbox_to_anchor=(1.02, 0.0), borderaxespad=0.0)

    ax.grid(True, axis="both", which="both", color='lightgray', linestyle='-', linewidth=0.5)
    ax.set_title("Total sequential runtime")
    ax.legend(loc='lower right')
    if len(polynomial_descriptions) > 0:
        print("Interpolating polynomials:\n" + "\n".join(polynomial_descriptions))
    if show:
        pyplot.tight_layout()
        pyplot.show()
        pyplot.close()




def plot_barchart(times, task_names, ax=None, show=True):
    '''
    For every number of samples, the fraction of total runtime for each task,
    represented as a horizontal barchart.
    '''
    y_ticks = list()
    y_labels = list()
    observed = set()
    colormap = matplotlib.colormaps["nipy_spectral"]
    if ax is None:
        fig, ax = pyplot.subplots()
    bar_height = 0.92
    for r,run_name in enumerate(times):
        left = 0
        y_pos = float(r)
        total_time = 0
        for t,task_name in enumerate(task_names):
            if task_name not in times[run_name]:
                continue
            if "total" in task_name:
                total_time += sum([t for w,t in times[run_name][task_name]])
        if total_time <= 0:
            continue
        for t,task_name in enumerate(task_names):
            if task_name not in times[run_name]:
                continue
            # y = sum([t for w,t in times[run_name][task_name] if w in common_windows[task_name]])
            y = sum([t for w,t in times[run_name][task_name]]) / total_time
            print('\t'.join(list(map(str,[run_name, task_name, y]))))
            l = str(t) + "_" + task_name
            color_index = (float(t) + 1) /float(len(task_names))
            color = colormap(color_index)
            x_left = left
            z = 1
            h = bar_height
            if "total" in task_name:
                x_left = 0
                z = 0
                color = "gray"
                h = 0.96
            if task_name not in observed:
                p = ax.barh(y_pos, y, h, label=task_name, left=x_left, color=color, zorder=z)
                observed.add(task_name)
            else:
                p = ax.barh(y_pos, y, h, left=x_left, color=color, zorder=z)
            if "total" not in task_name:
                left += y
        y_ticks.append(y_pos)
        y_labels.append(int(run_name[1:]))
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Fraction")
    ax.set_ylabel("Number of samples")
    ax.set_xlim(0, 1)
    n_bars = len(y_ticks)
    ax.set_ylim(-0.5, n_bars - 0.5)
    ax.legend(loc='lower left')
    ax.set_title("Composition of total sequential runtime")
    if show:
        pyplot.tight_layout()
        pyplot.show()
        pyplot.close()


def plot_barchart_original(times, task_names):
    x_ticks = list()
    x_labels = list()

    observed = set()

    colormap = matplotlib.colormaps["nipy_spectral"]

    fig = pyplot.figure()
    ax = pyplot.axes()

    for r,run_name in enumerate(times):
        bottom = 0

        for t,task_name in enumerate(task_names):
            if task_name not in times[run_name]:
                continue

            # y = sum([t for w,t in times[run_name][task_name] if w in common_windows[task_name]])
            y = sum([t for w,t in times[run_name][task_name]])

            print('\t'.join(list(map(str,[run_name, task_name, y]))))

            x = r + 0.5

            l = str(t) + "_" + task_name

            color_index = (float(t) + 1) /float(len(task_names))

            color = colormap(color_index)

            b = bottom
            z = 1
            w = 0.8
            if "total" in task_name:
                b = 0
                z = 0
                color = "gray"
                w = 0.82

            if task_name not in observed:
                p = ax.bar(x, y, w, label=l, bottom=b, color=color, zorder=z)
                observed.add(task_name)
            else:
                p = ax.bar(x, y, w, bottom=b, color=color, zorder=z)

            if "total" not in task_name:
                x_t = (float(x) - w/2.0) + (color_index * w/1.0)
                y_t = (b + b + y) / 2
                ax.text(x_t, y_t, "-" + str(t) + "-",va='center', ha='center')

                x_ticks.append(x)
                x_labels.append(run_name)

                bottom += y

    pyplot.xticks(x_ticks, x_labels, rotation=45, ha="right")

    ax.set_xlabel("Solver")
    ax.set_ylabel("Time (m)")

    pyplot.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    pyplot.tight_layout()

    pyplot.show()
    pyplot.close()




def plot_window_runtime_distribution_histogram(times, task="window_total", bins=100, normalized=True, ax=None, show=True):
    '''
    For a given task, the histogram of window runtimes.
    '''
    if ax is None:
        fig, ax = pyplot.subplots(figsize=(6, 6))
    run_values = dict()
    for run_name in sorted(times):
        if task not in times[run_name]:
            continue
        y = [t for _, t in times[run_name][task]]
        if len(y) == 0:
            continue
        run_values[run_name] = y
    for run_name, y in run_values.items():
        weights = [1.0 / float(len(y))] * len(y) if normalized else None
        counts, edges = np.histogram(y, bins=bins, weights=weights)
        centers = (edges[:-1] + edges[1:]) / 2.0
        ax.plot(centers, counts, marker='o', markersize=3, linewidth=0.0, label=task)
    ax.set_xlabel("Minutes")
    if normalized:
        ax.set_ylabel("Fraction of all windows")
    else:
        ax.set_ylabel("Number of windows")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_box_aspect(1)
    ax.set_title(f"Window runtimes in {int(run_name[1:])} samples")
    ax.grid(True, color='lightgray', linestyle='-', linewidth=0.5)
    ax.legend(loc='upper right')
    if show:
        pyplot.tight_layout()
        pyplot.show()
        pyplot.close()




def plot_graphaligner_vs_optimize(parent_dirs, task_x, task_y, bins=60):
    '''
    Tries to show two window runtimes on the same plot: each point is a window, 
    with x and y coordinates corresponding to the runtime of the two tasks.
    '''
    for directory in parent_dirs:
        task_times_by_window = defaultdict(dict)

        for root, dirs, files in sorted(os.walk(directory)):
            for file in files:
                file_path = os.path.join(root, file)
                if "log.csv" not in file_path:
                    continue

                window = root.strip('/').split("/")[-1]
                with open(file_path, 'r') as file_handle:
                    for l, line in enumerate(file_handle):
                        if l == 0:
                            continue

                        task, h, m, s, ms, success, notes = line.strip().split(',')
                        if task not in (task_x, task_y):
                            continue
                        if not bool(int(success)):
                            continue

                        s_total = float(h) * 60 * 60 + float(m) * 60 + float(s) + float(ms) / 1000.0
                        m_total = s_total / 60.0
                        task_times_by_window[window][task] = m_total

        x = list()
        y = list()
        for window, task_times in task_times_by_window.items():
            if task_x in task_times and task_y in task_times:
                x.append(task_times[task_x])
                y.append(task_times[task_y])

        if len(x) == 0:
            print("No successful windows with both tasks in", directory)
            continue

        fig = pyplot.figure()
        ax = pyplot.axes()
        cmap = matplotlib.colormaps["viridis"].copy()
        cmap.set_bad("white")
        h = ax.hist2d(x, y, bins=bins, cmap=cmap, norm=LogNorm(), cmin=1)

        min_v = min(min(x), min(y))
        max_v = max(max(x), max(y))
        if min_v == max_v:
            min_v -= 0.5
            max_v += 0.5

        ax.set_xlim(min_v, max_v)
        ax.set_ylim(min_v, max_v)
        ax.set_aspect("equal", adjustable="box")
        ax.plot([min_v, max_v], [min_v, max_v], color="lightgray", linewidth=1.0, zorder=3)

        pyplot.colorbar(h[3], ax=ax, label="Window count")
        ax.set_xlabel(task_x + " runtime (m)")
        ax.set_ylabel(task_y + " runtime (m)")
        ax.set_title(os.path.basename(directory) + ": " + task_x + " vs " + task_y)
        ax.grid(True, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.4)
        pyplot.tight_layout()
        pyplot.show()
        pyplot.close()




def get_failed_windows(parent_dirs, n_samples, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other):
    '''
    Remark: when `graphaligner` fails, `window_total` is marked as failed, but 
    when `optimize_d_plus_n` fails, `window_total` is not marked as failed.
    '''
    for directory in parent_dirs:
        for root, dirs, files in sorted(os.walk(directory)):
            for file in files:
                file_path = os.path.join(root, file)
                if "log.csv" in file_path:
                    out_total[n_samples] += 1
                    graphaligner_failed = False
                    optimize_failed = False
                    other_failed = False
                    with open(file_path, 'r') as file:    
                        for l,line in enumerate(file):
                            if l == 0:
                                continue
                            task, h, m, s, ms, success, notes = line.strip().split(',')
                            if int(success) == 0:
                                if task == "graphaligner":
                                    graphaligner_failed = True
                                elif task == "optimize_d_plus_n":
                                    optimize_failed = True
                                elif task != "window_total":
                                    other_failed = True
                    if graphaligner_failed:
                        out_failed_graphaligner[n_samples] += 1
                    elif optimize_failed:
                        out_failed_optimize_d_plus_n[n_samples] += 1
                    elif other_failed:
                        out_failed_other[n_samples] += 1




def plot_failed_windows(total, failed_graphaligner, failed_optimize_d_plus_n, failed_other, ax=None, show=True):
    '''
    A simple line plot showing, for each number of samples, the percentage of 
    failed windows.
    '''
    if ax is None:
        fig, ax = pyplot.subplots(figsize=(6, 6))
    x = sorted(set(total.keys()))

    # Graphaligner
    y = list()
    for n_samples in x:
        y.append(100.0*float(failed_graphaligner.get(n_samples,0)) / float(total.get(n_samples,0))) 
    ax.plot(x, y, marker='o', label="graphaligner")

    # Optimize_d_plus_n
    y = list()
    for n_samples in x:
        y.append(100.0*float(failed_optimize_d_plus_n.get(n_samples,0)) / float(total.get(n_samples,0))) 
    ax.plot(x, y, marker='o', label="optimize_d_plus_n")

    # # Other. No plotted since it is always zero in practice.
    # y = list()
    # for n_samples in x:
    #     y.append(100.0*float(failed_other.get(n_samples,0)) / float(total.get(n_samples,0))) 
    # ax.plot(x, y, marker='o', label="other")

    ax.set_xlabel("Number of samples")
    ax.set_ylabel("% of all windows")
    ax.set_title("Windows that exceed the maximum runtime")
    ax.set_box_aspect(1)
    ax.grid(True, color='lightgray', linestyle='-', linewidth=0.5)
    ax.legend()
    if show:
        pyplot.tight_layout()
        pyplot.show()
        pyplot.close()




def main():
    # Input directories
    n0001 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0001/tmp.4yhJQ1/tmp.LI4dRERe8f/output" ]
    n0002 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0002/tmp.lzEWtV/tmp.RvGaXcDD9p/output" ]
    n0004 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0004/tmp.2OGV9S/tmp.FxI82dHKTs/output" ] 
    n0008 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0008/tmp.DXG5WY/tmp.colQIV7d31/output" ]
    n0016 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0016/tmp.mWYem8/tmp.7U9bHvVnxo/output" ]
    n0032 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0032/tmp.pyTad9/tmp.k4r789DL4D/output" ]
    n0064 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0064/tmp.iUCVjm/tmp.PbpK9Shgtk/output" ]
    n0128 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0128/tmp.d20J8O/tmp.dGKvM8WJR1/output" ]
    n0256 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0256/tmp.7eVuZ0/tmp.BMOrl3n6dm/output" ]
    n0512 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n0512/tmp.62Jaof/tmp.s6v8nb8xU4/output" ]
    n1024 = [ "/Users/fcunial/Downloads/HAPESTRY/experimental_section_3/logs/n1024/tmp.w2HIMS/tmp.d2QV32P002/output" ]

    # Loading data used in Plots 1 and 2.
    # Remark: `times[name][task]` contains every window with `task` in its log,
    # including windows that failed on that task.
    times = defaultdict(lambda: defaultdict(list))
    common_windows = defaultdict(set)
    times, common_windows = aggregate_logs(parent_dirs=n0002, times=times, common_windows=common_windows, name="n0002")
    times, common_windows = aggregate_logs(parent_dirs=n0004, times=times, common_windows=common_windows, name="n0004")
    times, common_windows = aggregate_logs(parent_dirs=n0008, times=times, common_windows=common_windows, name="n0008")
    times, common_windows = aggregate_logs(parent_dirs=n0016, times=times, common_windows=common_windows, name="n0016")
    times, common_windows = aggregate_logs(parent_dirs=n0032, times=times, common_windows=common_windows, name="n0032")
    times, common_windows = aggregate_logs(parent_dirs=n0064, times=times, common_windows=common_windows, name="n0064")
    times, common_windows = aggregate_logs(parent_dirs=n0128, times=times, common_windows=common_windows, name="n0128")
    times, common_windows = aggregate_logs(parent_dirs=n0256, times=times, common_windows=common_windows, name="n0256")
    times, common_windows = aggregate_logs(parent_dirs=n0512, times=times, common_windows=common_windows, name="n0512")
    times, common_windows = aggregate_logs(parent_dirs=n1024, times=times, common_windows=common_windows, name="n1024")

    fig, axes = pyplot.subplots(2, 2, figsize=(12, 12))

    # 1. Runtime scaling plot.
    # Remark: we only plot the tasks with largest runtimes here. We show more 
    # tasks in Plot 2.
    task_names_scaling = [
        # "variant_graph",
        "graphaligner",
        "align_reads_to_paths",
        # "feasibility_construct",
        # "feasibility_init",
        #"feasibility",
        # "compress_transmap",
        # "optimize_d_prune_construct",
        # "optimize_d_prune_init",
        #"optimize_d_prune",
        # "optimize_d_construct",
        # "optimize_d_init",
        # "optimize_d",
        # "optimize_d_prune_parse",
        # "optimize_d_plus_n_construct",
        # "optimize_d_plus_n_init",
        "optimize_d_plus_n",
        # "optimize_d_plus_n_parse",
        "window_total",
    ]
    plot_scaling(times, task_names_scaling, True, ax=axes[0, 0], show=False)

    # 2. Fraction of total runtime plot
    # Remark: we show more tasks here, but not all for simplicity.
    task_names_barchart = [
        # "variant_graph",
        "graphaligner",
        "align_reads_to_paths",
        # "feasibility_construct",
        # "feasibility_init",
        "feasibility",
        # "compress_transmap",
        # "optimize_d_prune_construct",
        # "optimize_d_prune_init",
        "optimize_d_prune",
        # "optimize_d_construct",
        # "optimize_d_init",
        # "optimize_d",
        # "optimize_d_prune_parse",
        # "optimize_d_plus_n_construct",
        # "optimize_d_plus_n_init",
        "optimize_d_plus_n",
        # "optimize_d_plus_n_parse",
        "window_total",
    ]
    plot_barchart(times, task_names_barchart, ax=axes[0, 1], show=False)

    # 3. Window runtime distribution
    times = defaultdict(lambda: defaultdict(list))
    common_windows = defaultdict(set)
    # times, common_windows = aggregate_logs(parent_dirs=n0064, times=times, common_windows=common_windows, name="n0064")
    times, common_windows = aggregate_logs(parent_dirs=n1024, times=times, common_windows=common_windows, name="n1024")
    # plot_window_runtime_distribution_histogram(times,"graphaligner")
    plot_window_runtime_distribution_histogram(times, "optimize_d_plus_n", ax=axes[1, 0], show=False)
    #plot_window_runtime_distribution_histogram(times,"align_reads_to_paths")
    
    # 4. Failed windows
    out_total = defaultdict(int)
    out_failed_graphaligner = defaultdict(int)
    out_failed_optimize_d_plus_n = defaultdict(int)
    out_failed_other = defaultdict(int)
    get_failed_windows(n0001, 1, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0002, 2, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0004, 4, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0008, 8, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0016, 16, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0032, 32, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0064, 64, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0128, 128, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0256, 256, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n0512, 512, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    get_failed_windows(n1024, 1024, out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other)
    plot_failed_windows(out_total, out_failed_graphaligner, out_failed_optimize_d_plus_n, out_failed_other, ax=axes[1, 1], show=False)

    fig.tight_layout()
    pyplot.show()
    pyplot.close(fig)




if __name__ == "__main__":
    main()
