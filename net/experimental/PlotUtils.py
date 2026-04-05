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
    if 'peak' in kwargs.keys():
        if isinstance(kwargs['peak'], bool) and kwargs['peak']:
            plt.figtext(0.75, text_y, f"width: {PredPeak(kwargs['data']):.2f}")
        elif isinstance(kwargs['peak'], float):
            plt.figtext(0.75, text_y, f"peak: {kwargs['peak']:.2f}")
        text_y -= 0.05
    if 'width' in kwargs.keys():
        if isinstance(kwargs['width'], bool) and kwargs['width']:
            plt.figtext(0.75, text_y, f"width: {PredWidth(kwargs['data']):.2f}")
        elif isinstance(kwargs['width'], float):
            plt.figtext(0.75, text_y, f"width: {kwargs['width']:.2f}")
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

def PlotQuantileBinnedDistr(*,
                            y_true : NDArray, 
                            y_pred: NDArray,
                            plotdir: str | os.PathLike, 
                            weights : NDArray | None = None, 
                            bins: int = 10,
                            verbose: bool = False,
                            title: str | None = None,
                            xlabel: str | None = None,
                            ylabel: str | None = "Arbitrary Unit",
                            xrange: Tuple[float, float] = (-0.05, 1.05),
                            xticks: List[float] | None | NDArray[float] = None,
                            plot_name: str = "quant_binned_score"):
    if weights is None:
        weights = np.ones_like(y_true)

    signal_mask = y_true == 1  
    sig_score = y_pred[signal_mask]
    bkg_score = y_pred[~signal_mask]

    weights_sig = weights[signal_mask]
    weights_bkg = weights[~signal_mask]
    
    signal_quantiles = np.zeros(bins + 1)
    percentiles = np.linspace(0, 1, bins + 1)
    signal_quantiles[1:-1] = np.quantile(
        sig_score.ravel(), 
        percentiles[1:-1], 
        weights=weights_sig.ravel(), 
        method='inverted_cdf'
    )
    signal_quantiles[-1] = 1.0

    sig_counts = []
    bkg_counts = []
    sig_sumw2 = []
    bkg_sumw2 = []
    for i in range(bins):
        low = signal_quantiles[i]
        high = signal_quantiles[i + 1]
        sig_bin_mask = (sig_score < high) & (sig_score >= low)
        bkg_bin_mask = (bkg_score < high) & (bkg_score >= low)
        sig_count = np.sum(weights_sig[sig_bin_mask])
        bkg_count = np.sum(weights_bkg[bkg_bin_mask])
        sig_counts.append(sig_count)
        bkg_counts.append(bkg_count)
        sig_sumw2.append(np.sum(weights_sig[sig_bin_mask]**2))
        bkg_sumw2.append(np.sum(weights_bkg[bkg_bin_mask]**2))
        if verbose:
            print(f"Bin {i} [{low:.2f}, {high:.2f}]: N(sig)={sig_count:.2f}, N(bkg)={bkg_count:.2f}")

    xs = []
    for i in range(bins):
        start = percentiles[i]
        end = percentiles[i + 1]
        xc = (end + start)/2
        xs.append(xc)

    xs = np.array(xs)
    sig_sumw2 = np.array(sig_sumw2)
    bkg_sumw2 = np.array(bkg_sumw2)
    sig_err = np.sqrt(sig_sumw2)
    bkg_err = np.sqrt(bkg_sumw2)
    plt.errorbar(
        x=xs, 
        y=sig_counts, 
        yerr=sig_err, 
        marker='s', 
        label="signal", 
        linestyle='none',
        markersize=3,
        color="blue"
    )
    plt.errorbar(
        x=xs, 
        y=bkg_counts, 
        yerr=bkg_err, 
        marker='o', 
        label="background", 
        linestyle='none',
        markersize=3,
        color="red"
    )
    if title:
        plt.title(title, fontsize=14)
    if xlabel: 
        plt.xlabel(xlabel=xlabel, fontsize=12)
    plt.ylabel(ylabel=ylabel, fontsize=12)
    plt.yscale('log')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()

    if xticks is not None:
        plt.xticks(xticks)
    
    # Adjust x-axis to fit the score range [0, 1]
    xmin, xmax = xrange
    plt.xlim(xmin, xmax)
    plt.tight_layout()
    out_path = os.path.join(plotdir, f'{plot_name}.pdf')
    plt.savefig(out_path, bbox_inches='tight')
    plt.clf()
    plt.close()
