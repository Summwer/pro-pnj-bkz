# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time
import math
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from typing import Optional, Sequence


from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma

# Display notation used in the original figures.  The numerical value is
# computed stably from the log-GS profiles as exp(real_log - sim_log).
RATIO_LATEX = (
    r'$\mathrm{Real}(\Vert\mathbf{b}_{i}^{*}\Vert)/'
    r'\mathrm{Sim}(\Vert\mathbf{b}_{i}^{*}\Vert)$'
)

# Read GS lengths.
def read_gs_lengths(filename):
    """Read positive GS lengths and return both logarithmic and raw profiles."""
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError("GS-length file was not found: %s" % path)

    rows = []
    with path.open('r', encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            row = [float(value) for value in fields]
            if any(value <= 0.0 for value in row):
                raise ValueError(
                    "All GS lengths must be positive; invalid value in %s at line %d"
                    % (path, line_number)
                )
            rows.append(row)

    if not rows:
        raise ValueError("No GS-length data were found in %s" % path)

    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("Inconsistent row lengths in %s" % path)

    log_rows = [[math.log(value) for value in row] for row in rows]
    return log_rows, rows


def read_blocksizes(file_name):
    """Read a whitespace-separated block-size list."""
    path = Path(file_name)
    if not path.is_file():
        raise FileNotFoundError("Block-size file was not found: %s" % path)
    with path.open('r', encoding='utf-8') as handle:
        return [int(value) for value in handle.read().split()]


def read_file(n, alpha_, jump, Pumpdown, Blocksize, Tours, test_number, d):
    directory = Path(
        'gs-lengths-simulator'
    ) / (
        'n=%d,alpha=%s,jump=%d,Pumpdown=%s,Blocksize=%d,Tours=%dd=%d,'
        % (n, alpha_, jump, Pumpdown, Blocksize, Tours, d)
    )
    filename = directory / ('rr_set,test_number=%d.txt' % test_number)
    log_GS_lengths, GS_lengths = read_gs_lengths(filename)
    return log_GS_lengths, GS_lengths, str(directory)


#When jump=1, PnJBKZ simulator degenerates to CN11 simulator
def pnjBKZ_simulator(log_rr0,beta,N,d,jump):
    l = []
    for i in range(d):
        temp = log_rr0[i]
        l.append(temp)
    #print(len(l))
    l_ = []
    J = jump
    randomSquaredAverages = [2.95747208 , 2.905419231 , 2.808209731 , 2.704122291 , 2.589140986 , 2.473617892 , 2.363752978 , 2.254288677 , 2.160395606 , 2.056417068 , 1.962571394 , 1.86329581 , 1.781471788 , 1.691206302 , 1.607855519 , 1.526258286 , 1.45048576 , 1.377169955 , 1.304800173 , 1.238714637 , 1.172164783 , 1.11055388 , 1.04872671 , 0.9933346386 , 0.9385826564 , 0.8858378414 , 0.8346947761 , 0.7869433841 , 0.743196743 , 0.6999110693 , 0.6612589966 , 0.620890347 , 0.5840212199 , 0.5494461008 , 0.5153120689 , 0.4851862344 , 0.4542255962 , 0.426615261 , 0.4011004324 , 0.3761683243 , 0.3521434119 , 0.3305054658 , 0.3109962405 , 0.2976205498 , 0.2897254705]
    #print(len(randomSquaredAverages))
    rk = []
    for i in range(len(randomSquaredAverages)):
        RK = log(math.sqrt(randomSquaredAverages[i]))
        rk.append(RK)
    #print(rk)
    cd = []
    if jump/beta > 0.95:
        for i in range(12):
            cd.append((lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2)*0.5)
        for i in range(12,beta):
            cd.append(lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2)
    
    elif jump/beta > 0.6:
        for i in range(12):
            cd.append((lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2)*0.5)
        for i in range(12,beta):
            cd.append(lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2)
        
    #print(cd)
    else:
        for i in range(beta):
            cd.append(lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2)
    print("PnJ-BKZ Simulator")
    for tours in range(N):
        flag = True
        for k in range(d - beta):
            beta_ = min(beta, d - k)
            if k % J == 0:
                f = min(k + beta , d)
                #之前的 f = min(k + beta, d)
                sumf = 0
                sumk = 0
                for i in range(f):
                    sumf = sumf + l[i]
                for i in range(k):
                    sumk = sumk + l_[i]
                logV = sumf - sumk
                if flag == True:
                    if logV / beta_ + cd[beta_-1] < l[k]:
                        l_k = logV / beta_ + cd[beta_-1]
                        if tours == 0:
                            l_.append(l_k)
                        else:
                            l_[k] = copy.deepcopy(l_k)
                        flag = False
                        # print("\n%d",k)
                        # print(l_[k] - l_[k-1])
                    else:
                        l_k = copy.deepcopy(l[k])
                        if tours == 0:
                            l_.append(l_k)
                        else:
                            l_[k] = copy.deepcopy(l_k)
                else:
                    l_k = logV / beta_ + cd[beta_-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = copy.deepcopy(l_k)
                if jump/beta > 0.12:
                    l_[k] = copy.deepcopy((l_[k-1] + l_k)/2)
                if k > 1 and jump/beta > 0.12:
                    l_[k] = copy.deepcopy((l_[k-1] + l_k)/2)
            else:
                f = min(k - (k % J) + beta, d)
                sumf = 0
                sumk = 0
                for i in range(f):
                    sumf = sumf + l[i]
                for i in range(k):
                    sumk = sumk + l_[i]
                logV = sumf - sumk
                if flag == True:
                    if logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1] < l[k]:
                        l_k = logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1]
                        if tours == 0:
                            l_.append(l_k)
                        else:
                            l_[k] = copy.deepcopy(l_k)
                            flag = False
                    else:
                        l_k = l[k]
                        if tours == 0:
                            l_.append(l_k)
                        else:
                            l_[k] = copy.deepcopy(l_k)
                else:
                    l_k = logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = copy.deepcopy(l_k)
                if jump/beta >= 0.5:
                    if k % J < 6:
                        l_[k] = copy.deepcopy((l_[k-1] + l_k)/2)
                        #l_[k] = copy.deepcopy((l_[k-2] + l_[k-1] + l_k)/3)
                if l_[k-1] - l_[k] > 0.09 and jump/beta > 0.12:
                    l_[k] = copy.deepcopy((l_[k-1] + l_k)/2)
        for k in range(d - beta,d - 45):
            beta_ = d - k
            f = d
            sumf = 0
            sumk = 0
            for i in range(f):
                sumf = sumf + l[i]
            for i in range(k):
                sumk = sumk + l_[i]
            logV = sumf - sumk
            if flag == True:
                if logV / beta_ + cd[beta_-1] < l[k]:
                    l_k = logV / beta_ + cd[beta_-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = copy.deepcopy(l_k)
                    flag = False
                else:
                    l_k = l[k]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = copy.deepcopy(l_k)
            else:
                l_k = logV / beta_ + cd[beta_-1]
                if tours == 0:
                    l_.append(l_k)
                else:
                    l_[k] = copy.deepcopy(l_k)
            if k < d - beta + 6 and jump/beta > 0.12:
                l_[k] = copy.deepcopy((l_[k-1] + l_k)/2)
            #print(l_)
        k_ = min (len(randomSquaredAverages), beta)
        #last 45 norms
        sumf=0;
        sumk=0;
        for i in range(d):
            sumf = sumf + l[i]
        for i in range(d - 45):
            sumk = sumk + l_[i]
        logV = sumf - sumk
        for k in range(d - k_, d):
            l_k = logV/k_ + rk[k + 45 - d]
            if tours == 0:
                l_.append(l_k)
            else:
                l_[k] = copy.deepcopy(l_k)
                #Set l[i] as l_[i]
            # if l_[k] > l_[k-1]:
            #     l_[k] = copy.deepcopy(l_[k-1] - 0.143841036)
        for i in range(d):
            l[i] = copy.deepcopy(l_[i])
    return l_


def CN11_simulator(log_rr0,beta,N,d):
    l = []
    for i in range(d):
        temp = log_rr0[i]
        l.append(temp)
    #print(len(l))
    l_ = []
    
    randomSquaredAverages = [2.95747208 , 2.905419231 , 2.808209731 , 2.704122291 , 2.589140986 , 2.473617892 , 2.363752978 , 2.254288677 , 2.160395606 , 2.056417068 , 1.962571394 , 1.86329581 , 1.781471788 , 1.691206302 , 1.607855519 , 1.526258286 , 1.45048576 , 1.377169955 , 1.304800173 , 1.238714637 , 1.172164783 , 1.11055388 , 1.04872671 , 0.9933346386 , 0.9385826564 , 0.8858378414 , 0.8346947761 , 0.7869433841 , 0.743196743 , 0.6999110693 , 0.6612589966 , 0.620890347 , 0.5840212199 , 0.5494461008 , 0.5153120689 , 0.4851862344 , 0.4542255962 , 0.426615261 , 0.4011004324 , 0.3761683243 , 0.3521434119 , 0.3305054658 , 0.3109962405 , 0.2976205498 , 0.2897254705]
    #print(len(randomSquaredAverages))
    rk = []
    for i in range(len(randomSquaredAverages)):
        RK = log(math.sqrt(randomSquaredAverages[i]))
        rk.append(RK)
    #print(rk)
    k_ = min (len(randomSquaredAverages), beta)
    cd = []
    for i in range(beta):
        CD = lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2
        cd.append(CD)
    #print(cd)
    print("CN11")
    for j in range(N):
        flag = True
        for k in range(d):
            beta_ = min(beta, d - k)
            f = min(k + beta, d)
            sumf = 0
            sumk = 0
            for i in range(f):
                sumf = sumf + l[i]
            for i in range(k):
                sumk = sumk + l_[i]
            logV = sumf - sumk
            if flag == True:
                if logV / beta_ + cd[beta_-1] < l[k]:
                    l_k = logV / beta_ + cd[beta_-1]
                    if j == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = copy.deepcopy(l_k)
                    flag = False
                else:
                    l_k = copy.deepcopy(l[k])
                    if j == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = copy.deepcopy(l_k)
            else:
                l_k = logV / beta_ + cd[beta_-1]
                if j == 0:
                    l_.append(l_k)
                else:
                    l_[k] = copy.deepcopy(l_k)
            #print(l_)
        #Set l[i] as l_[i]
        for i in range(d):
            l[i] = copy.deepcopy(l_[i])
    return l_

def _ensure_same_profile_length(observed_log_gs, simulated_log_gs):
    observed = np.asarray(observed_log_gs, dtype=float)
    simulated = np.asarray(simulated_log_gs, dtype=float)
    if observed.ndim != 1 or simulated.ndim != 1:
        raise ValueError("GS profiles must be one-dimensional")
    if observed.shape != simulated.shape:
        raise ValueError(
            "Observed and simulated profiles have different lengths: %d and %d"
            % (observed.size, simulated.size)
        )
    return observed, simulated


def real_to_simulated_norm_ratio(real_log_gs, simulated_log_gs):
    """Return Real(||b_i^*||) / Sim(||b_i^*||).

    Both inputs are log-GS profiles.  The ratio shown in the figures is the
    direct GS-norm ratio, while the numerically stable implementation is

        exp(real_log_gs - simulated_log_gs).
    """
    real_profile, simulated_profile = _ensure_same_profile_length(
        real_log_gs, simulated_log_gs
    )
    return np.exp(real_profile - simulated_profile)


# Backward-compatible name retained for any external code importing it.
def observed_to_simulated_norm_ratio(observed_log_gs, simulated_log_gs):
    return real_to_simulated_norm_ratio(observed_log_gs, simulated_log_gs)


def average_ratio(ratio):
    """Arithmetic mean of the Real/Sim GS-norm ratios."""
    ratio_array = np.asarray(ratio, dtype=float)
    if ratio_array.size == 0:
        raise ValueError("Cannot compute an average ratio from an empty array")
    return float(np.mean(ratio_array))


def mean_absolute_relative_error(ratio):
    """Mean absolute relative error from the ideal ratio 1."""
    ratio_array = np.asarray(ratio, dtype=float)
    if ratio_array.size == 0:
        raise ValueError("Cannot compute MARE from an empty ratio array")
    return float(np.mean(np.abs(ratio_array - 1.0)))


def ratio_axis_limits(ratio):
    """Include every ratio and both accuracy reference bands, with padding."""
    values = np.asarray(ratio, dtype=float).ravel()
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Ratio axis limits require finite, nonempty data")
    lower = min(0.90, float(values.min()))
    upper = max(1.10, float(values.max()))
    padding = max(0.01, 0.05 * (upper - lower))
    return lower - padding, upper + padding


def ratio_file_path(directory, n, alpha_, dimension, blocksize, jump, current_tour):
    return Path(directory) / (
        'ratio, n=%d, alpha = %s, dimension = %d, Blocksize = %d, '
        'Jump = %d, Current Tours= %d.txt'
        % (n, alpha_, dimension, blocksize, jump, current_tour)
    )


def show_gs_slope_figure(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error,
    blocksize,
    jump,
    tour_index,
    test_number,
):
    """Plot the observed and simulated log-GS profiles using the original style."""
    observed, simulated = _ensure_same_profile_length(
        log_gs_length, sim_log_gs_lengths
    )
    t = 0
    observed_x = [_ + 1 for _ in range(t, observed.size)]
    simulated_x = [_ + 1 for _ in range(t, simulated.size)]

    plt.figure(figsize=(15, 10), dpi=100)
    plt.scatter(
        observed_x,
        observed[t:],
        marker='o',
        c='none',
        edgecolors='b',
    )
    plt.scatter(
        simulated_x,
        simulated[t:],
        marker='x',
        c='r',
    )
    plt.rcParams.update({'font.size': 20})
    plt.legend(['Experimental', 'jump simulator'])

    output = Path(directory) / (
        'n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, '
        'Current Tours= %d.png'
        % (n, alpha_, dimension, blocksize, jump, tour_index + 1)
    )
    plt.savefig(output)
    plt.close()
    return output


def show_gs_slope_figure2(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    cn11_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error1,
    square_error2,
    blocksize,
    jump,
    tour_index,
    test_number,
):
    """Compare observed, PnJBKZ-simulated, and CN11 profiles in the original style."""
    observed, simulated = _ensure_same_profile_length(
        log_gs_length, sim_log_gs_lengths
    )
    cn11 = np.asarray(cn11_log_gs_lengths, dtype=float)
    if cn11.shape != observed.shape:
        raise ValueError('CN11 and observed profiles have different lengths')

    t = 0
    x = [_ + 1 for _ in range(t, dimension)]
    plt.figure(figsize=(15, 10), dpi=100)
    plt.scatter(
        x,
        observed[t:dimension],
        marker='^',
        c='none',
        edgecolors='k',
    )
    plt.scatter(
        [_ + 1 for _ in range(t, simulated.size)],
        simulated[t:],
        marker='x',
        c='r',
    )
    plt.scatter(
        [_ + 1 for _ in range(t, cn11.size)],
        cn11[t:],
        marker='o',
        c='none',
        edgecolors='b',
    )
    plt.xlabel('Index')
    plt.ylabel('log GS norm')
    plt.legend(
        ['Experimental', 'PnJBKZ simulator (LLL)', 'BKZ2.0 simulator']
    )

    output = Path(directory) / (
        'Comparing figure, n=%d, alpha = %s, dimension = %d, '
        'Blocksize = %d, Jump = %d, Current Tours= %d.png'
        % (n, alpha_, dimension, blocksize, jump, tour_index + 1)
    )
    plt.savefig(output)
    plt.close()
    return output


def show_gs_slope_figure3(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error1,
    blocksize,
    jump,
    tour_index,
    test_number,
):
    """Plot the observed and PnJBKZ-simulated profiles in the original style."""
    observed, simulated = _ensure_same_profile_length(
        log_gs_length, sim_log_gs_lengths
    )
    t = 0
    plt.figure(figsize=(15, 10), dpi=100)
    plt.scatter(
        [_ + 1 for _ in range(t, dimension)],
        observed[t:dimension],
        marker='^',
        c='none',
        edgecolors='k',
    )
    plt.scatter(
        [_ + 1 for _ in range(t, simulated.size)],
        simulated[t:],
        marker='x',
        c='r',
    )
    plt.xlabel('Index', fontsize=24)
    plt.ylabel('log GS norm', fontsize=24)
    plt.legend(['Experimental', 'PnJBKZ simulator (LLL)'], fontsize=24)
    plt.title(
        'Dimension = %d, Blocksize = %d, Jump = %d, '
        'Current Tours= %d,SimError = %f'
        % (
            dimension,
            blocksize,
            jump,
            tour_index + 1,
            square_error1,
        ),
        fontsize=20,
    )
    plt.tick_params(labelsize=22)

    output = Path(directory) / (
        'Comparing figure, n=%d, alpha = %s, dimension = %d, '
        'Blocksize = %d, Jump = %d, Current Tours= %d.png'
        % (n, alpha_, dimension, blocksize, jump, tour_index + 1)
    )
    plt.savefig(output)
    plt.close()
    return output


def show_gs_slope_figure4(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error1,
    blocksize,
    jump,
    tour_index,
    test_number,
):
    """Plot one-tour GS-norm ratios with all data points in the axis range."""
    ratio = real_to_simulated_norm_ratio(
        log_gs_length, sim_log_gs_lengths
    )
    x = [_ + 1 for _ in range(ratio.size)]
    avg_ratio = average_ratio(ratio)
    mare = mean_absolute_relative_error(ratio)

    print(
        'Average Ratio %f, Ratio Abs Average %f, Current Tours = %d'
        % (avg_ratio, mare, tour_index + 1)
    )

    plt.figure(figsize=(8, 10), dpi=100)
    plt.scatter(
        x,
        ratio,
        s=12,
        marker='*',
        c='k',
        edgecolors='k',
    )
    plt.ylim(*ratio_axis_limits(ratio))
    plt.plot([1, dimension], [1, 1], c='k')
    plt.plot([1, dimension], [1.05, 1.05], c='k', linestyle='--')
    plt.plot([1, dimension], [0.95, 0.95], c='k', linestyle='--')
    plt.legend([RATIO_LATEX], fontsize=16)
    plt.xlabel('Index', fontsize=16)
    plt.ylabel('Ratio', fontsize=16)
    plt.title(
        'Dimension = %d, Blocksize = %d, Jump = %d, Current Tours= %d'
        % (dimension, blocksize, jump, tour_index + 1),
        fontsize=14,
    )

    output = Path(directory) / (
        'Ratio figure, n=%d, alpha = %s, dimension = %d, '
        'Blocksize = %d, Jump = %d, Current Tours= %d.png'
        % (n, alpha_, dimension, blocksize, jump, tour_index + 1)
    )
    plt.savefig(output)
    plt.close()
    return output


def ratio_txt(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error1,
    blocksize,
    jump,
    tour_index,
    test_number,
):
    """Save Real(||b_i^*||) / Sim(||b_i^*||) for one tour."""
    ratio = real_to_simulated_norm_ratio(
        log_gs_length, sim_log_gs_lengths
    )
    if ratio.size != dimension:
        raise ValueError(
            "Expected %d ratio entries, found %d" % (dimension, ratio.size)
        )

    current_tour = tour_index + 1
    output = ratio_file_path(
        directory,
        n,
        alpha_,
        dimension,
        blocksize,
        jump,
        current_tour,
    )
    np.savetxt(output, ratio.reshape(1, -1), fmt='%.16e')
    return output


def load_ratio_file(
    directory,
    n,
    alpha_,
    dimension,
    blocksize,
    jump,
    current_tour,
):
    path = ratio_file_path(
        directory,
        n,
        alpha_,
        dimension,
        blocksize,
        jump,
        current_tour,
    )
    if not path.is_file():
        raise FileNotFoundError(
            "Ratio data for tour %d were not found: %s"
            % (current_tour, path)
        )
    ratio = np.atleast_1d(np.loadtxt(path, dtype=float))
    if ratio.size != dimension:
        raise ValueError(
            "Ratio file %s contains %d values; expected %d"
            % (path, ratio.size, dimension)
        )
    return ratio


def select_representative_tours(total_tours):
    """Select the first, middle, and last tours without duplicates."""
    if total_tours < 1:
        raise ValueError("total_tours must be at least 1")
    return sorted({1, (total_tours + 1) // 2, total_tours})


def plot_ratio_across_tours(
    directory,
    n,
    dimension,
    alpha_,
    blocksizes,
    jump,
    total_tours,
    selected_tours=None,
):
    """Plot first/middle/last tour ratios while preserving the original style.

    The numerical fixes are retained:
      * the selected tour files are distinct;
      * ratios are Real(||b_i^*||) divided by Sim(||b_i^*||);
      * both average ratio and mean absolute relative error are reported.

    The original colors, marker shapes, shaded band, line styles, and figure
    size are preserved.  The axis range includes every ratio, and the legend
    is placed automatically to reduce overlap with data points.
    """
    if len(blocksizes) != total_tours:
        raise ValueError(
            'blocksizes must contain one entry per tour: expected %d, found %d'
            % (total_tours, len(blocksizes))
        )

    tours = (
        select_representative_tours(total_tours)
        if selected_tours is None
        else sorted(set(int(tour) for tour in selected_tours))
    )
    if not tours or tours[0] < 1 or tours[-1] > total_tours:
        raise ValueError('selected_tours must lie in [1, total_tours]')
    if len(tours) > 3:
        raise ValueError(
            'The original visual style supports at most three selected tours'
        )

    ratios = []
    statistics = []
    for tour in tours:
        blocksize = blocksizes[tour - 1]
        ratio = load_ratio_file(
            directory,
            n,
            alpha_,
            dimension,
            blocksize,
            jump,
            tour,
        )
        ratios.append(ratio)

        avg_ratio = average_ratio(ratio)
        mare = mean_absolute_relative_error(ratio)
        statistics.append((tour, blocksize, avg_ratio, mare))
        print(
            'Average Ratio %f, Ratio Abs Average %f, Current Tours = %d'
            % (avg_ratio, mare, tour)
        )

    x = [_ for _ in range(dimension)]
    plt.figure(figsize=(16, 9), dpi=100)
    plt.xlabel('Index', fontsize=24)
    plt.ylabel('Ratio', fontsize=24)
    plt.fill_between(
        [0, dimension],
        0.951,
        1.049,
        facecolor='gainsboro',
        alpha=0.7,
    )

    original_styles = [
        dict(s=24, marker='*', c='k', edgecolors='k'),
        dict(s=16, marker='^', c='none', edgecolors='b'),
        dict(s=16, marker='o', c='r', edgecolors='r'),
    ]

    for ratio, tour, style in zip(ratios, tours, original_styles):
        plt.scatter(
            x,
            ratio,
            label=(
                RATIO_LATEX[:-1]
                + r',\;\mathrm{Tours}=%d$' % tour
            ),
            **style
        )

    plt.ylim(*ratio_axis_limits(np.concatenate(ratios)))
    plotted_dimension = len(ratios[0])
    plt.plot([1, plotted_dimension], [1, 1], c='k')
    plt.plot(
        [1, plotted_dimension],
        [1.1, 1.1],
        c='k',
        linestyle=':',
    )
    plt.plot(
        [1, plotted_dimension],
        [0.9, 0.9],
        c='k',
        linestyle=':',
    )
    plt.plot(
        [1, plotted_dimension],
        [1.05, 1.05],
        c='k',
        linestyle='--',
    )
    plt.plot(
        [1, plotted_dimension],
        [0.95, 0.95],
        c='k',
        linestyle='--',
    )
    plt.tick_params(labelsize=20)
    plt.legend(fontsize=20, loc='best')

    unique_blocksizes = sorted(set(blocksizes))
    if len(unique_blocksizes) == 1:
        title_blocksize = unique_blocksizes[0]
        filename_blocksize = 'Blocksize = %d' % unique_blocksizes[0]
        title = 'Dimension = %d, Blocksize = %d, Jump = %d, Average Error Ratio = %.2f%%' % (
            plotted_dimension,
            title_blocksize,
            jump,
            100. * mare
        )
    else:
        filename_blocksize = 'Variable blocksizes'
        title = 'Dimension = %d, Variable Blocksizes, Jump = %d' % (
            plotted_dimension,
            jump,
        )

    plt.title(title, fontsize=20)

    output = Path(directory) / (
        'Ratio with different tours, Dimension = %d, %s, Jump = %d.png'
        % (plotted_dimension, filename_blocksize, jump)
    )
    plt.savefig(output)
    plt.close()

    stats_output = output.with_name(output.stem + ', statistics.txt')
    with stats_output.open('w', encoding='utf-8') as handle:
        handle.write('tour\tblocksize\taverage_ratio\tMARE_percent\n')
        for tour, blocksize, avg_ratio, mare in statistics:
            handle.write(
                '%d\t%d\t%.12f\t%.8f\n'
                % (tour, blocksize, avg_ratio, 100.0 * mare)
            )

    return output


def Ratio_figure(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error1,
    blocksize,
    jump,
    total_tours,
    test_number,
):
    """Backward-compatible wrapper for a constant block size."""
    return plot_ratio_across_tours(
        directory=directory,
        n=n,
        dimension=dimension,
        alpha_=alpha_,
        blocksizes=[blocksize] * total_tours,
        jump=jump,
        total_tours=total_tours,
    )


def Ratio_figure_high_jump(
    directory,
    log_gs_length,
    sim_log_gs_lengths,
    n,
    dimension,
    alpha_,
    square_error1,
    blocksize,
    jump,
    total_tours,
    test_number,
):
    """Compatibility wrapper; the generic plotter also handles high jumps."""
    return Ratio_figure(
        directory,
        log_gs_length,
        sim_log_gs_lengths,
        n,
        dimension,
        alpha_,
        square_error1,
        blocksize,
        jump,
        total_tours,
        test_number,
    )


def compute_square_error(list1, list2, flag=1):
    """Squared error in the log-GS domain."""
    first = np.asarray(list1, dtype=float)
    second = np.asarray(list2, dtype=float)
    overlap = min(first.size, second.size)
    start = 0
    if flag == 1 and second.size > 1 and second[1] - second[0] > 1:
        start = 1

    square_error = float(np.sum((first[start:overlap] - second[start:overlap]) ** 2))
    if first.size > overlap:
        square_error += float(np.sum(first[overlap:] ** 2))
    elif second.size > overlap:
        square_error += float(np.sum(second[overlap:] ** 2))
    return square_error


#=========================================input==============================================

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Experiment configuration
    # ------------------------------------------------------------------
    
    
    
    
    
    
    
    # n, d, alpha_, Pumpdown = 70, 263, "010", "True"
    
    
    # n, d, alpha_, Pumpdown = 50, 194, "015", "True"
    
    
    
    
    #60-010
    #n, d, alpha_, Pumpdown = 60, 222, "010", "True"
    #slope = -0.04963
    # Total_test_number = 10
    # Blocksize, jump, Tours, test_number, extra_dim4free = 80, 1, 9, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 80, 20, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 80, 40, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 80, 80, 13, 1, 0
    #Total_test_number = 20
    # Blocksize, jump, Tours, test_number, extra_dim4free = 85, 1, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 85, 12, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 100, 1, 5, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 100, 12, 13, 1, 0
    
    
    #75-005
    # n, d, alpha_, Pumpdown = 75, 252, "005", "True"
    # slope = -0.04339
    # Total_test_number = 20
    # Blocksize, jump, Tours, test_number, extra_dim4free = 75, 1, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 80, 3, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 80, 6, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 85, 3, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 100, 6, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 85, 9, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 95, 12, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 95, 9, 13, 1, 0
    
    
    #70-005
    n, d, alpha_, Pumpdown = 70, 235, "005", "True"
    slope = -0.04927
    Total_test_number = 20
    # Blocksize, jump, Tours, test_number, extra_dim4free = 70, 3, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 60, 6, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 65, 6, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 70, 6, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 60, 9, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 65, 9, 13, 1, 0
    # Blocksize, jump, Tours, test_number, extra_dim4free = 70, 9, 13, 1, 0
    Blocksize, jump, Tours, test_number, extra_dim4free = 90, 12, 13, 1, 0
    
    #50-015
    #slope = -0.05127
    
    
    
    
    
    
    
    Blocksizes = [Blocksize for _ in range(Tours)]
    Sim_Blocksizes = [Blocksize for _ in range(Tours)]

    # Example of a variable block-size schedule:
    # Blocksizes = [80, 82, 84, 86, 88, 90]
    # Sim_Blocksizes = Blocksizes.copy()
    # Tours = len(Blocksizes)

    # ------------------------------------------------------------------
    # d4f helper functions
    # ------------------------------------------------------------------
    def Theory1_d4f(beta):
        return math.floor(beta * log(4 / 3) / log(beta / (2 * pi * exp(1))))

    def Theory2_d4f(beta):
        return math.floor(beta * log(4 / 3) / log(beta / (2 * pi)))

    def accs_2023_d4f(current_slope):
        return math.floor(log(sqrt(4 / 3)) / (-current_slope / 4))

    def d4f_gap(beta):
        if beta < 40:
            gap = 0
        elif beta < 76:
            gap = (beta - 40) / 2 - Theory1_d4f(beta)
        else:
            gap = (11.5 + 0.075 * beta) - Theory1_d4f(beta)
        return round(gap)

    def d4f_gap2(beta):
        if beta < 40:
            gap = 0
        elif beta < 76:
            gap = (beta - 40) / 2 - Theory2_d4f(beta)
        else:
            gap = (11.5 + 0.075 * beta) - Theory2_d4f(beta)
        return round(gap)

    def d4f_gap_accs(beta, current_slope):
        if beta < 40:
            gap = 0
        elif beta < 75:
            gap = (beta - 40) / 2 - min(
                Theory1_d4f(beta), accs_2023_d4f(current_slope)
            )
        else:
            gap = (11.5 + 0.075 * beta) - min(
                Theory1_d4f(beta), accs_2023_d4f(current_slope) - 1
            )
        return round(gap)

    # Keep the currently selected simulation adjustment from the original code.
    d4f_value = "No d4f value and no extra d4f"
    for index, beta in enumerate(Sim_Blocksizes):
        if jump / float(Blocksize) < 0.75:
            Sim_Blocksizes[index] = beta - 3
        else:
            Sim_Blocksizes[index] = beta
    print("d4f setting: %s" % d4f_value)

    # ------------------------------------------------------------------
    # Average the observed profiles over all independent experiments.
    # Each input file is expected to contain Tours+1 profiles: the initial
    # profile followed by one output profile for each tour.
    # ------------------------------------------------------------------
    
    experiment_profiles = []
    directory = None

    for experiment_number in range(1, Total_test_number + 1):
        profiles, _, directory = read_file(
            n,
            alpha_,
            jump,
            Pumpdown,
            Blocksize,
            Tours,
            experiment_number,
            d,
        )
        profile_array = np.asarray(profiles, dtype=float)
        expected_shape = (Tours + 1, d)
        if profile_array.shape != expected_shape:
            raise ValueError(
                "Experiment %d has profile shape %s; expected %s"
                % (experiment_number, profile_array.shape, expected_shape)
            )
        experiment_profiles.append(profile_array)

    # G6K stores 2*log(||b_i^*||); divide by two after averaging.
    log_GS_lengths = np.mean(
        np.stack(experiment_profiles, axis=0), axis=0
    ) / 2.0

    # ------------------------------------------------------------------
    # Simulate one tour at a time, save per-tour ratios, and create figures.
    # ------------------------------------------------------------------
    mid_log_gs = None
    cn11_mid_log_gs = None

    for tour_index, (blocksize, sim_blocksize) in enumerate(
        zip(Blocksizes, Sim_Blocksizes)
    ):
        input_profile = (
            log_GS_lengths[0]
            if tour_index == 0
            else mid_log_gs
        )
        cn11_input_profile = (
            log_GS_lengths[0]
            if tour_index == 0
            else cn11_mid_log_gs
        )

        if blocksize > 44:
            sim_log_gs_lengths = np.asarray(
                pnjBKZ_simulator(
                    input_profile.tolist(), sim_blocksize, 1, d, jump
                ),
                dtype=float,
            )
        else:
            sim_log_gs_lengths = np.asarray(
                CN11_simulator(input_profile.tolist(), sim_blocksize, 1, d),
                dtype=float,
            )

        cn11_sim = np.asarray(
            CN11_simulator(cn11_input_profile.tolist(), blocksize, 1, d),
            dtype=float,
        )

        mid_log_gs = sim_log_gs_lengths.copy()
        cn11_mid_log_gs = cn11_sim.copy()
        observed_profile = log_GS_lengths[tour_index + 1]

        square_error1 = compute_square_error(
            sim_log_gs_lengths, observed_profile, 1
        )
        square_error2 = compute_square_error(
            cn11_sim, observed_profile, 1
        )
        print(
            "\nTour %d: our squared log-profile error = %.6f; "
            "CN11 error = %.6f"
            % (tour_index + 1, square_error1, square_error2)
        )

        show_gs_slope_figure3(
            directory,
            observed_profile,
            sim_log_gs_lengths,
            n,
            d,
            alpha_,
            square_error1,
            blocksize,
            jump,
            tour_index,
            test_number,
        )
        show_gs_slope_figure4(
            directory,
            observed_profile,
            sim_log_gs_lengths,
            n,
            d,
            alpha_,
            square_error1,
            blocksize,
            jump,
            tour_index,
            test_number,
        )
        ratio_txt(
            directory,
            observed_profile,
            sim_log_gs_lengths,
            n,
            d,
            alpha_,
            square_error1,
            blocksize,
            jump,
            tour_index,
            test_number,
        )

    # Compare tours 1, middle, and last.  For Tours=5, this is 1, 3, and 5.
    combined_ratio_figure = plot_ratio_across_tours(
        directory=directory,
        n=n,
        dimension=d,
        alpha_=alpha_,
        blocksizes=Blocksizes,
        jump=jump,
        total_tours=Tours,
    )
    print("Combined ratio figure saved to: %s" % combined_ratio_figure)
