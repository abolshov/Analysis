import matplotlib.pyplot as plt
import os
import numpy as np
from MiscUtils import PredPeak, PredWidth, target_names


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

def PlotCompare2D(target, output, quantity, plotting_dir=None):
    plt.grid(False)
    min_bin = 0 if quantity[-1] == 'E' else -1200
    bins = np.linspace(min_bin, 1200, 100)
    plt.hist2d(target, output, bins=bins)
    var = quantity.split('_')[-1]
    plt.title(f'{var} comparison')
    plt.ylabel(f'predicted {target_names[quantity]}')
    plt.xlabel(f'true {target_names[quantity]}')
    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f"cmp2d_{quantity}.pdf"), bbox_inches='tight')
    else:
        plt.savefig(f"cmp2d_{quantity}.pdf", bbox_inches='tight')
    plt.clf()


def PlotHist(bins=np.linspace(0, 250, 100), pos_frac=None, **kwargs):
    plt.hist(kwargs['data'], bins=bins)
    plt.title(kwargs['title'])
    plt.ylabel(kwargs['ylabel'])
    plt.xlabel(kwargs['ylabel'])
    plt.figtext(0.75, 0.8, f"peak: {PredPeak(kwargs['data']):.2f}")
    plt.figtext(0.75, 0.75, f"width: {PredWidth(kwargs['data']):.2f}")
    if pos_frac:
        plt.figtext(0.75, 0.7, f"pos: {pos_frac:.2f}")
    plt.grid(True)
    if '.' in kwargs['file_name']:
        plt.savefig(os.path.join(kwargs['plotting_dir'], kwargs['file_name']), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(kwargs['plotting_dir'], f"{kwargs['file_name']}.pdf"), bbox_inches='tight')
    plt.clf()