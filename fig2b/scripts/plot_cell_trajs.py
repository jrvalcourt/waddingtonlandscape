import sys
from PIL import Image
import os
import numpy as np
import math
import pickle
import scipy.stats
import math
from sklearn.metrics import mutual_info_score
from sklearn.feature_selection import mutual_info_regression
from sklearn import mixture
import scipy.stats as stats
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import colors as mplcolors
import seaborn as sb

# define some global constants for this time lapse
START_FRAME = 71
END_FRAME = 402
SIGNAL_FRAME = 302
DIFF_START_FRAME = 86
MIN_PER_FRAME = 15
MIN_PER_HOUR = 60
LOG_DIVIDING_LINE = -1.1
DIVIDING_LINE = 2 ** LOG_DIVIDING_LINE
END_T = (END_FRAME * MIN_PER_FRAME - 
        SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR
START_T = (START_FRAME * MIN_PER_FRAME - 
        SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR
DIFF_START_T = (DIFF_START_FRAME * MIN_PER_FRAME - 
        SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR
SIGNAL_T = 0

###############################################################################
#
# function definitions
#
###############################################################################

def yelblu_cmap(n, line=-1.1):
    if n > line:
        return (209/256.0, 171/256.0, 43/256.0)
    else:
        return (40/256.0, 126/256.0, 194/256.0)

def log2(n):
    return np.log(n) / np.log(2)

def fsigmoid(x, a, b):
   return 1.0 / (1.0 + np.exp(-a*(x-b)))

def calc_mi(cs_at_signal_uncommitted, cs_at_signal_committed, log_div_line,
        pseudocount = np.finfo(float).eps):
    n_unc = len(cs_at_signal_uncommitted)
    n_com = len(cs_at_signal_committed)
    
    counts = {}
    counts['unc'] = {}
    counts['com'] = {}
    counts['unc']['unc'] = sum([1 for x in cs_at_signal_uncommitted \
                                if x >= log_div_line])
    counts['unc']['com'] = sum([1 for x in cs_at_signal_uncommitted \
                                if x <  log_div_line])
    counts['com']['com'] = sum([1 for x in cs_at_signal_committed \
                                if x <  log_div_line])
    counts['com']['unc'] = sum([1 for x in cs_at_signal_committed \
                                if x >= log_div_line])

    total = 0
    for a in counts:
        for b in counts[a]:
            total += counts[a][b]
    total += 4 * pseudocount
    
    prob_pre  = {}
    prob_pre['unc'] = (counts['unc']['unc'] + counts['unc']['com'] \
            + 2*pseudocount) / total
    prob_pre['com'] = (counts['com']['unc'] + counts['com']['com'] \
            + 2*pseudocount) / total
    
    prob_post = {}
    prob_post['unc'] = (counts['unc']['unc'] + counts['com']['unc'] \
            + 2*pseudocount) / total
    prob_post['com'] = (counts['unc']['com'] + counts['com']['com'] \
            + 2*pseudocount) / total
    
    prob = {}
    prob['unc'] = {}
    prob['com'] = {}
    for a in counts:
        for b in counts[a]:
            prob[a][b] = (counts[a][b] + pseudocount) / total
    
    
    mi = 0
    for a in prob:
        for b in prob:
            mi += prob[a][b] * np.log(prob[a][b] / \
                    (prob_pre[a] * prob_post[b])) / np.log(2)
    return mi

def plot_osr(data):
    fig1, ax1 = plt.subplots()
    plt.xlabel('time after BMP4 + Activin A addition (h)')
    plt.ylabel('log2 [OCT4:RFP / SOX2:YFP]')
    cs_at_end = []
    
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        ratio = [a[4] / a[5] for a in data[track]]
        if not END_T in t:
            continue
        c_at_end = ratio[t.index(END_T)]
        cs_at_end.append(c_at_end)
        color = yelblu_cmap(c_at_end, line=DIVIDING_LINE)
        plt.plot(t, np.log(ratio)/np.log(2), c=color)
    plt.plot([0] * 2, [-100000,100000], 'k--')
    plt.plot([DIFF_START_T] * 2, [-100000,100000], 'k--')
    ylim = [-5, 2.5]
    plt.ylim(ylim)
    plt.xlim([START_T, END_T])
    plt.savefig("plots/log2_ratio_over_time.png", bbox_inches='tight', dpi=1000)
    plt.close()

def plot_oct4(data):
    fig1, ax1 = plt.subplots()
    plt.xlabel('time after BMP4 + Activin A addition (h)')
    plt.ylabel('OCT4:RFP')
    cs_at_end = []
    
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        ratio = [a[4] / a[5] for a in data[track]]
        oct4  = [a[4] for a in data[track]]
        if not END_T in t:
            continue
        c_at_end = ratio[t.index(END_T)]
        cs_at_end.append(c_at_end)
        color = yelblu_cmap(c_at_end, line=DIVIDING_LINE)
        plt.plot(t, oct4, c=color)
    plt.plot([0] * 2, [-100000,100000], 'k--')
    plt.plot([DIFF_START_T] * 2, [-100000,100000], 'k--')
    ylim = [0, 1200]
    plt.ylim(ylim)
    plt.xlim([START_T, END_T])
    plt.savefig("plots/oct4_over_time.png", bbox_inches='tight', dpi=1000)
    plt.close()

def plot_sox2(data):
    fig1, ax1 = plt.subplots()
    plt.xlabel('time after BMP4 + Activin A addition (h)')
    plt.ylabel('SOX2:YFP')
    cs_at_end = []
    
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        ratio = [a[4] / a[5] for a in data[track]]
        sox2  = [a[5] for a in data[track]]
        if not END_T in t:
            continue
        c_at_end = ratio[t.index(END_T)]
        cs_at_end.append(c_at_end)
        color = yelblu_cmap(c_at_end, line=DIVIDING_LINE)
        plt.plot(t, sox2, c=color)
    plt.plot([0] * 2, [-100000,100000], 'k--')
    plt.plot([DIFF_START_T] * 2, [-100000,100000], 'k--')
    ylim = [0, 1200]
    plt.ylim(ylim)
    plt.xlim([START_T, END_T])
    plt.savefig("plots/sox2_over_time.png", bbox_inches='tight', dpi=1000)
    plt.close()

def plot_histograms(data):
    plt.figure(figsize=(8,4))
    cs_at_end = []
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        ratio = [a[4] / a[5] for a in data[track]]
        if not END_T in t:
            continue
        c_at_end = ratio[t.index(END_T)]
        cs_at_end.append(c_at_end)
    log_cs_at_end = np.log(cs_at_end)/np.log(2)
    yellow_cs = log_cs_at_end[log_cs_at_end >  LOG_DIVIDING_LINE]
    blue_cs   = log_cs_at_end[log_cs_at_end <= LOG_DIVIDING_LINE]
    ylim = [-5, 2.5]
    plt.hist(blue_cs,   bins=np.linspace(ylim[0], ylim[1], num=20), 
            orientation='vertical', density=True, 
            color=(40/256.0, 126/256.0, 194/256.0))
    plt.hist(yellow_cs, bins=np.linspace(ylim[0], ylim[1], num=20), 
            orientation='vertical', density=True, 
            color=(209/256.0, 171/256.0, 43/256.0))
    xs = np.linspace(ylim[0], ylim[1], 200)
    clf = mixture.GaussianMixture(n_components=2)
    clf.fit((np.log(cs_at_end) / np.log(2)).reshape(-1,1))
    means, covs = zip(*sorted(zip(clf.means_, clf.covariances_), 
        key=lambda x: x[0]))
    plt.xlim(ylim)
    plt.yticks([0,0.5], ['0','0.5'])
    plt.savefig("plots/log2_ratio_hist.png", dpi=300)
    plt.close()
    
    addition_bins = np.linspace(ylim[0], ylim[1], num=50)
    cs_at_signal_committed   = []
    cs_at_signal_uncommitted = []
    oct4_at_signal_committed   = []
    sox2_at_signal_committed   = []
    oct4_at_signal_uncommitted = []
    sox2_at_signal_uncommitted = []
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        logratio = [np.log(a[4] / a[5])/np.log(2) for a in data[track]]
        oct4 = [a[4] for a in data[track]]
        sox2 = [a[5] for a in data[track]]
        cs_before_signal = [c for time, c in zip(t, logratio) if time <= 0]
        oct4_before_signal = [c for time, c in zip(t, oct4) if time <= 0]
        sox2_before_signal = [c for time, c in zip(t, sox2) if time <= 0]
        if not END_T in t:
            continue
        if len(cs_before_signal) < 3:
            continue
        end_idx = t.index(END_T)
        end_c = logratio[end_idx]
        c_at_signal = np.mean(cs_before_signal[-1:])
        oct4_at_signal = np.mean(oct4_before_signal[-1:])
        sox2_at_signal = np.mean(sox2_before_signal[-1:])
        if end_c > LOG_DIVIDING_LINE:
            cs_at_signal_uncommitted.append(c_at_signal)
            oct4_at_signal_uncommitted.append(oct4_at_signal)
            sox2_at_signal_uncommitted.append(sox2_at_signal)
        else:
            cs_at_signal_committed.append(c_at_signal)
            oct4_at_signal_committed.append(oct4_at_signal)
            sox2_at_signal_committed.append(sox2_at_signal)
   
    n_uncommitted, _, _ = plt.hist(cs_at_signal_uncommitted, 
            bins=addition_bins, density=True)
    n_committed,   _, _ = plt.hist(cs_at_signal_committed,   
            bins=addition_bins, density=True)
    plt.close()
    plt.close()
    
    bar_width = (addition_bins[1] - addition_bins[0]) * 0.95
    plt.figure()
    plt.bar(addition_bins[:-1], height=n_uncommitted, width=bar_width, 
            align='edge', color=(209/255, 171/255, 43/255))
    plt.bar(addition_bins[:-1], height=n_committed,   width=bar_width, 
            align='edge', color=(40/255,  126/255, 194/255))
    plt.xlim([max(addition_bins), min(addition_bins)])
    plt.savefig("plots/stacked_bar_color_at_signal_addition.png", 
            bbox_inches='tight', dpi=300)
    plt.close()
    
def get_osr_pairs_for_mi(data):
    cs_at_signal_committed   = []
    cs_at_signal_uncommitted = []
    oct4_at_signal_committed   = []
    sox2_at_signal_committed   = []
    oct4_at_signal_uncommitted = []
    sox2_at_signal_uncommitted = []
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        logratio = [np.log(a[4] / a[5])/np.log(2) for a in data[track]]
        oct4 = [a[4] for a in data[track]]
        sox2 = [a[5] for a in data[track]]

        cs_before_signal = [c for time, c in zip(t, logratio) if time <= 0]
        oct4_before_signal = [c for time, c in zip(t, oct4) if time <= 0]
        sox2_before_signal = [c for time, c in zip(t, sox2) if time <= 0]

        if not END_T in t:
            continue
        if len(cs_before_signal) < 3:
            continue

        end_idx = t.index(END_T)
        end_c = logratio[end_idx]

        c_at_signal = np.mean(cs_before_signal[-1:])
        oct4_at_signal = np.mean(oct4_before_signal[-1:])
        sox2_at_signal = np.mean(sox2_before_signal[-1:])

        if end_c > LOG_DIVIDING_LINE:
            cs_at_signal_uncommitted.append(c_at_signal)
            oct4_at_signal_uncommitted.append(oct4_at_signal)
            sox2_at_signal_uncommitted.append(sox2_at_signal)
        else:
            cs_at_signal_committed.append(c_at_signal)
            oct4_at_signal_committed.append(oct4_at_signal)
            sox2_at_signal_committed.append(sox2_at_signal)
    return cs_at_signal_uncommitted, cs_at_signal_committed

def get_oct4_pairs_for_mi(data):
    oct4_at_signal_committed   = []
    oct4_at_signal_uncommitted = []
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        logratio = [np.log(a[4] / a[5])/np.log(2) for a in data[track]]
        oct4 = [a[4] for a in data[track]]

        cs_before_signal = [c for time, c in zip(t, logratio) if time <= 0]
        oct4_before_signal = [c for time, c in zip(t, oct4) if time <= 0]

        if not END_T in t:
            continue
        if len(cs_before_signal) < 3:
            continue

        end_idx = t.index(END_T)
        end_c = logratio[end_idx]

        c_at_signal = np.mean(cs_before_signal[-1:])
        oct4_at_signal = np.mean(oct4_before_signal[-1:])

        if end_c > LOG_DIVIDING_LINE:
            oct4_at_signal_uncommitted.append(oct4_at_signal)
        else:
            oct4_at_signal_committed.append(oct4_at_signal)
    return oct4_at_signal_uncommitted, oct4_at_signal_committed

def get_sox2_pairs_for_mi(data):
    sox2_at_signal_committed   = []
    sox2_at_signal_uncommitted = []
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        logratio = [np.log(a[4] / a[5])/np.log(2) for a in data[track]]
        sox2 = [a[5] for a in data[track]]

        cs_before_signal = [c for time, c in zip(t, logratio) if time <= 0]
        sox2_before_signal = [c for time, c in zip(t, sox2) if time <= 0]

        if not END_T in t:
            continue
        if len(cs_before_signal) < 3:
            continue

        end_idx = t.index(END_T)
        end_c = logratio[end_idx]

        sox2_at_signal = np.mean(sox2_before_signal[-1:])

        if end_c > LOG_DIVIDING_LINE:
            sox2_at_signal_uncommitted.append(sox2_at_signal)
        else:
            sox2_at_signal_committed.append(sox2_at_signal)
    return sox2_at_signal_uncommitted, sox2_at_signal_committed

def plot_cell_cycle(data, cell_divs):
    fig1, ax1 = plt.subplots()
    plt.xlabel('fraction through cell cycle')
    plt.ylabel('responds to signal')
    cs_at_end = []    
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) / MIN_PER_HOUR \
                for a in data[track]]
        ratio = [a[4] / a[5] for a in data[track]]
        if not END_T in t:
            continue
        if not track in cell_divs:
            continue
        divs = [((x - SIGNAL_FRAME) * MIN_PER_FRAME) / MIN_PER_HOUR \
                for x in cell_divs[track]]
        div_before = None
        div_after = None
        for d in sorted(divs):
            if d < SIGNAL_T:
                div_before = d
            if (d >= SIGNAL_T) and (div_after is None):
                div_after = d
                break
        if (div_before is None) or (div_after is None):
            continue
        frac_through_at_signal = (SIGNAL_T - div_before) / \
                (div_after - div_before)
        c_at_end = ratio[t.index(END_T)]
        cs_at_end.append(c_at_end)
        color = yelblu_cmap(c_at_end, line=DIVIDING_LINE)
        plt.plot(frac_through_at_signal, c_at_end > DIVIDING_LINE, 
                'o', c=color)
    ylim = [-0.05, 1.05]
    plt.ylim(ylim)
    plt.xlim([-0.05,1.05])
    plt.savefig("plots/cell_cycle.png", bbox_inches='tight', dpi=1000)
    plt.close()

if __name__ == '__main__':

    # load the data
    data      = pickle.load(open(sys.argv[1], "rb"))
    cell_divs = pickle.load(open(sys.argv[2], "rb"))
    
    # plot OSR, OCT4, and SOX2 over time
    plot_osr(data)
    plot_oct4(data)
    plot_sox2(data)
    
    # plot OSR histograms at signal addtion and at the end of the time lapse
    plot_histograms(data)

    # cell divisions
    plot_cell_cycle(data, cell_divs)

    ###########################################################################
    #
    # calculate the mutual info for OSR or OCT4 or SOX2 with final fate
    #
    ###########################################################################

    osr_at_signal_uncommitted, osr_at_signal_committed = \
            get_osr_pairs_for_mi(data)
    oct4_at_signal_uncommitted, oct4_at_signal_committed = \
            get_oct4_pairs_for_mi(data)
    sox2_at_signal_uncommitted, sox2_at_signal_committed = \
            get_sox2_pairs_for_mi(data)
    
    print(f"MI with final fate")
    mi = calc_mi(osr_at_signal_uncommitted, osr_at_signal_committed,
            LOG_DIVIDING_LINE)
    print(f"OSR: {mi:0.3f} bits")
    mi = calc_mi(oct4_at_signal_uncommitted, oct4_at_signal_committed, 228)
    print(f"OCT4: {mi:0.3f} bits")
    mi = calc_mi(sox2_at_signal_uncommitted, sox2_at_signal_committed, 575)
    print(f"SOX2: {mi:0.3f} bits")

    ###########################################################################
    # 
    #
    #
    ###########################################################################
    
    color_pairs = []
    cells_at_signal = []
    es_colors = []
    es_pairs = []
    for track in data:
        t     = [(a[1] - SIGNAL_FRAME * MIN_PER_FRAME) \
                / MIN_PER_HOUR for a in data[track]]
        logred   = [np.log(a[4]) / np.log(2) for a in data[track]]
        loggreen = [np.log(a[5]) / np.log(2) for a in data[track]]
        logratio = [np.log(a[4] / a[5]) / np.log(2) for a in data[track]]
        cs_before_signal = [c for time, c in zip(t, logratio) if time <= 0]
        logred_before_signal = [c for time, c in zip(t, logred) if time <= 0]
        loggreen_before_signal = [c for time, c in zip(t, loggreen) \
                if time <= 0]
        cs_during_es = [c for time, c in zip(t, logratio) \
                if time <= DIFF_START_T]
        logred_during_es = [c for time, c in zip(t, logred) \
                if time <= DIFF_START_T]
        loggreen_during_es = [c for time, c in zip(t, loggreen) \
                if time <= DIFF_START_T]
        if not END_T in t:
            continue
        if len(cs_before_signal) < 1:
            continue
        if len(cs_during_es) > 0:
            es_colors.append(np.mean(cs_during_es))
        if len(cs_during_es) > 0:
            es_pairs.append((np.mean(logred_during_es), 
                np.mean(loggreen_during_es)))
        end_idx = t.index(END_T)
        end_c = logratio[end_idx]
        c_at_signal = np.mean(cs_before_signal[-1:])
        logred_at_signal = np.mean(logred_before_signal[-1:])
        loggreen_at_signal = np.mean(loggreen_before_signal[-1:])
        if end_c > LOG_DIVIDING_LINE:
            final_color = 1
        else:
            final_color = 0
        color_pairs.append((c_at_signal, final_color))
        cells_at_signal.append((logred_at_signal, loggreen_at_signal, 
            final_color))
    
    ###########################################################################
    # 
    # fit and plot a sigmoid curve to p(mesendo|OSR)
    #
    ###########################################################################

    plt.figure()
    xdata = [x[0] for x in color_pairs]
    ydata = [x[1] for x in color_pairs]
    mean_color_atsignal = np.mean(xdata)
    mean_es_color = np.mean(es_colors)
    std_es_color = np.std(es_colors)
    xs = np.linspace(-5, 2, num=200)
    popt, pcov = curve_fit(fsigmoid, xdata, ydata, 
                           method='trf', 
                           bounds=([0.0, -5.0], [30.0, 1.0])
                          )
    
    # now draw a bunch of sigmoid parameters that match 
    # the mean and covariance matrix we estimated
    # see https://stats.stackexchange.com/questions/ \
    #        120179/generating-data-with-a-given-sample-covariance-matrix
    n_curves = 1000
    X = np.matmul(np.random.randn(n_curves, 2), np.linalg.cholesky(pcov));
    popt_tile = np.tile(popt, (X.shape[0], 1))
    X = X + popt_tile
    # trim any sigmoids outside of 1 sigma
    perr = np.sqrt(np.diag(pcov))
    X = X[np.logical_and(X[:,0] < popt[0] + perr[0], 
        X[:,0] > popt[0] - perr[0]), :]
    X = X[np.logical_and(X[:,1] < popt[1] + perr[1], 
        X[:,1] > popt[1] - perr[1]), :]
    
    es_rect = Rectangle(xy=(mean_es_color - std_es_color, -0.5), 
            width=2 * std_es_color,
            height=2,
            edgecolor=None, 
            fc=(16/255, 167/255, 74/255, 1.0))
    plt.gca().add_artist(es_rect)
    max_sigmoid = None
    min_sigmoid = None
    for ii in range(X.shape[0]):
        plt.plot(xs, [fsigmoid(x, X[ii,0], X[ii,1]) for x in xs], 
                color=(0.3,0.3,0.3,1/255))
    plt.plot(xs, [fsigmoid(x, popt[0], popt[1]) for x in xs], 'k')
    plt.plot(xdata, ydata, 'o', markerfacecolor="None", markeredgecolor='k')
    plt.ylim([-0.05,1.05])
    plt.xlim([-3.5,0.5])
    plt.savefig("plots/sigmoid.png", bbox_inches='tight', dpi=300)
    plt.close()
