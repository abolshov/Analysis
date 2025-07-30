import matplotlib.pyplot as plt
import os
import numpy as np
from MiscUtils import PredPeak, PredWidth, ground_truth_map, pretty_vars, objects


def PlotMetric(history, model, metric, plotting_dir=None):
    plt.plot(history.history[metric], label=f'train_{metric}')
    plt.plot(history.history[f'val_{metric}'], label=f'val_{metric}')
    plt.title(f'{model} {metric}')
    plt.ylabel(metric)
    plt.xlabel('Epoch')
    plt.legend(loc='upper right')
    plt.grid(True)
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f"{metric}_{model}.pdf"), bbox_inches='tight')
    else:
        plt.savefig(f"{metric}_{model}.pdf", bbox_inches='tight')
    plt.clf()


def PlotCompare2D(target, output, var, obj, bins=None, title='', xlabel='', ylabel='', quantile=None, plotting_dir=None):
    plt.grid(False)

    plot_name_tokens = ['cmp2d', obj, var]
    if quantile:
        plot_name_tokens.append(f'q{quantile}')
    
    if bins is None:
        min_target = np.min(target)
        min_output = np.min(output)
        bin_left = np.min([min_target, min_output]) - 1.0

        max_target = np.max(target)
        max_output = np.max(output)
        bin_right = np.max([max_target, max_output]) + 1.0

        bins = np.linspace(bin_left, bin_right, 100)

    plt.hist2d(target, output, bins=bins)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plot_name = '_'.join(plot_name_tokens)
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f'{plot_name}.pdf'), bbox_inches='tight')
    else:
        plt.savefig(f'{plot_name}.pdf', bbox_inches='tight')
    plt.clf()


def PlotHist(bins=np.linspace(0, 250, 100), **kwargs):
    plt.hist(kwargs['data'], bins=bins)
    plt.title(kwargs['title'])
    plt.ylabel(kwargs['ylabel'])
    plt.xlabel(kwargs['xlabel'])
    plt.grid(True)

    text_y = 0.8
    if 'peak' in kwargs.keys() and kwargs['peak']:
        plt.figtext(0.75, text_y, f"peak: {PredPeak(kwargs['data']):.2f}")
        text_y -= 0.05
    if 'width' in kwargs.keys() and kwargs['width']:
        plt.figtext(0.75, text_y, f"width: {PredWidth(kwargs['data']):.2f}")
        text_y -= 0.05
    if 'count' in kwargs.keys() and kwargs['count']:
        plt.figtext(0.75, text_y, f"count: {len(kwargs['data'])}")
        text_y -= 0.05
    if 'pos_frac' in kwargs.keys():
        plt.figtext(0.75, text_y, f"pos: {kwargs['pos_frac']:.2f}")
        text_y -= 0.05

    if '.' in kwargs['file_name']:
        plt.savefig(os.path.join(kwargs['plotting_dir'], kwargs['file_name']), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(kwargs['plotting_dir'], f"{kwargs['file_name']}.pdf"), bbox_inches='tight')
    plt.clf()