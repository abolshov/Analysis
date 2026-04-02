import matplotlib.pyplot as plt
import os
import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_curve, precision_recall_curve, ConfusionMatrixDisplay
from MiscUtils import PredPeak, PredWidth, ground_truth_map, pretty_vars, objects

from numpy.typing import NDArray
from typing import Literal, List, Tuple, Dict, Any

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
    plt.close()


def PlotCovarMtrx(mtrx, method, labels, plotdir, display_values=True):
    plt.matshow(mtrx)
    plt.xticks(ticks=np.arange(mtrx.shape[0]), labels=labels, fontsize=5, rotation=45)
    plt.yticks(ticks=np.arange(mtrx.shape[0]), labels=labels, fontsize=5, rotation=45)

    if display_values:
        for (i, j), val in np.ndenumerate(mtrx):
            plt.text(j, i, f'{val:.1e}', ha='center', va='center', color='white', fontsize=5)

    cb = plt.colorbar()
    plt.title('Covariance Matrix')
    ax = plt.gca()
    ax.xaxis.set_ticks_position('bottom')
    plt.savefig(os.path.join(plotdir, f'{method}_cov_mtrx.pdf'), bbox_inches='tight')
    plt.clf()
    plt.close()

def PlotHistStack(*,
                  data: List[NDArray], 
                  hist_params: Dict[str, Dict[str, Any]], 
                  file_name: str | os.PathLike | pathlib.Path, 
                  val_range: Tuple[float, float],
                  weights: List[NDArray] = None, 
                  density: bool = False, 
                  stacked: bool = False, 
                  num_bins: int = 50, 
                  title: str = None, 
                  xlabel: str = None, 
                  ylabel: str = None, 
                  plotting_dir: str | os.PathLike | pathlib.Path = None) -> None:
    
    assert len(data) == len(hist_params), f'Mismatch between number of input arrays ({len(data)}) and sets of parameters ({len(hist_params)})'

    begin, end = val_range
    bins = np.linspace(begin, end, num_bins)

    for i, (label, params) in enumerate(hist_params.items()):
        plt.hist(data[i], 
                 label=label, 
                 weights=weights[i] if weights else None,
                 bins=bins, 
                 density=density,
                 range=val_range,
                 color=params['color'],
                 histtype='step',
                 stacked=stacked,
                 linewidth=params['linewidth'])

    if title:
        plt.title(title)
    if ylabel:
        plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    plt.grid(True)
    plt.legend()

    if plotting_dir:
        plt.savefig(os.path.join(plotting_dir, f'{file_name}.pdf'), bbox_inches='tight')
    else:
        plt.savefig(f'{file_name}.pdf', bbox_inches='tight')
    plt.clf()
    plt.close()

def PlotConfMatrix(*,
                   labels: NDArray, 
                   predictions: NDArray,
                   threshold: float,
                   plotdir: str | os.PathLike | pathlib.Path,
                   normalize: Literal['true', 'pred', 'all'] ='true',
                   tick_labels: List[str] = None,
                   file_name: str = None):
    cm = confusion_matrix(labels, predictions > threshold, normalize=normalize)
    if tick_labels is None:
        tick_labels = np.unique(labels)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=tick_labels)
    disp.plot(cmap=plt.cm.viridis, values_format='.4f')
    plt.title('Confusion matrix @{:.2f}'.format(threshold))
    if file_name is None:
        plt.savefig(os.path.join(plotdir, 'conf_mtrx.pdf'), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(plotdir, f'{file_name}.pdf'), bbox_inches='tight')
    plt.clf()
    plt.close()

def PlotROC(*,
            labels, 
            predictions,
            plotdir,
            xlim=[-0.5, 20],
            ylim=[80, 100.5],
            **kwargs):
    fp, tp, _ = roc_curve(labels, predictions)

    plt.plot(100*fp, 100*tp, linewidth=2, **kwargs)
    plt.title('Reciever Operation Characteristic (ROC)')
    plt.xlabel('False positives [%]')
    plt.ylabel('True positives [%]')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(os.path.join(plotdir, 'roc.pdf'), bbox_inches='tight')
    plt.clf()
    plt.close()

# implement ROC comparator:
# takes lists of labels, predictions and plots on the same canvas
def PlotCompareROC():
    pass

def PlotPRC(*,
            labels, 
            predictions,
            plotdir,
            **kwargs):
    "Plots Precision Recall Curve => PRC"
    precision, recall, _ = precision_recall_curve(labels, predictions)

    plt.plot(precision, recall, linewidth=2, **kwargs)
    plt.title('Precision-Recall Curve (PRC)')
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(os.path.join(plotdir, 'prc.pdf'), bbox_inches='tight')
    plt.clf()
    plt.close()

# implement PRC comparator:
# takes lists of labels, predictions and plots on the same canvas
def PlotComparePRC():
    pass

def PlotPrecisionAtK(*,
                     y_true : NDArray, 
                     y_scores: NDArray,
                     plot_dir: str | os.PathLike, 
                     weights : NDArray | None = None, 
                     bins: int = 10):
    """
    Plots weighted counts using horizontal lines for bin width, 
    markers for the mean, and vertical lines for error bars.
    """
    if weights is None:
        weights = np.ones_like(y_true)
        
    df = pd.DataFrame({
        'label': y_true,
        'score': y_scores,
        'weight': weights,
        'weight_sq': np.square(weights)
    })
    
    # 1. Calculate Quantiles and Bin Ranges
    # We need the actual score boundaries for the horizontal lines
    df['bin'], bin_edges = pd.qcut(df['score'], q=bins, labels=False, retbins=True, duplicates='drop')
    
    # 2. Calculate Weighted Counts and Errors
    counts = df.groupby(['bin', 'label'])['weight'].sum().unstack(fill_value=0)
    errors = np.sqrt(df.groupby(['bin', 'label'])['weight_sq'].sum().unstack(fill_value=0))
    
    plt.figure(figsize=(12, 7))
    colors = {0: '#e74c3c', 1: '#2ecc71'}
    labels = {0: 'Background', 1: 'Signal'}
    
    # 3. Plot each class separately
    for label in [0, 1]:
        if label not in counts.columns:
            continue
            
        # X-coordinates: the midpoint of the bin for the marker
        bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Y-coordinates: the weighted counts
        y_vals = counts[label].values
        y_errs = errors[label].values
        
        # Plot vertical error bars and markers
        plt.errorbar(bin_midpoints, y_vals, yerr=y_errs, fmt='o', 
                     color=colors[label], label=labels[label],
                     capsize=0, elinewidth=2, markersize=8, markeredgecolor='white')
        
        # Plot horizontal lines representing the bin width
        for i in range(len(bin_edges) - 1):
            plt.hlines(y=y_vals[i], xmin=bin_edges[i], xmax=bin_edges[i+1], 
                       color=colors[label], linewidth=3, alpha=0.4)

    plt.title('DNN Score', fontsize=14)
    plt.xlabel('DNN Score', fontsize=12)
    if weights:
        plt.ylabel('Weighted Count', fontsize=12)
    else:
        plt.ylabel('Count', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()
    
    # Adjust x-axis to fit the score range [0, 1]
    plt.xlim(max(0, bin_edges.min() - 0.05), min(1, bin_edges.max() + 0.05))
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'precision_at_k.pdf'), bbox_inches='tight')
    plt.clf()
    plt.close()
