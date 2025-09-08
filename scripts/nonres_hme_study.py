import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

def EqualBinningPlot(x, y, var_name,
                     output_name,
                     n_bins=20, 
                     statistic='mean',
                     color="tab:blue", 
                     title=None,
                     xlabel=None,
                     ylabel=None):
    """
    Plot mean ± IQR (q16–q84) of rel_delta vs mass using equal-population bins.

    Parameters:
        x : array-like
            The array of mass values (x-axis).
        y : array-like
            The array of y values (y-axis).
        n_bins : int
            Number of bins (default 20).
        color : str
            Color for the plot and shaded IQR band.
        title : str
            Optional title for the plot.
    """
    x = np.asarray(x)
    t = np.asarray(t)

    # Equal-population (quantile) bin edges
    bins = np.quantile(mass, np.linspace(0, 1, n_bins + 1))
    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    # Bin assignment and chosen variable
    st, _, binnumber = binned_statistic(x, y, statistic=statistic, bins=bins)

    # Compute q16 and q84 for each bin
    q16 = np.full(n_bins, np.nan)
    q84 = np.full(n_bins, np.nan)
    for i in range(1, n_bins + 1):
        mask = binnumber == i
        if np.any(mask):
            q16[i-1], q84[i-1] = np.percentile(y[mask], [16, 84])

    # Plot mean and IQR band
    plt.plot(bin_centers, st, marker="o", color=color, label=f"{statistic.capitalize()} {var_name}", linestyle='--')
    plt.fill_between(bin_centers, q16, q84, color=color, alpha=0.3, label="IQR (q16–q84)")

    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    plt.legend()
    plt.savefig(f'{output_name}.pdf', bbox_inches='tight')
    plt.clf()
    plt.close()
