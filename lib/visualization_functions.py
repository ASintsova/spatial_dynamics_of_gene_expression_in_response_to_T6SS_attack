
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
import os
import scipy as ss
import scipy.stats
#%matplotlib inline
import itertools
import seaborn as sns
from sklearn.decomposition import PCA

def invnorm(x):
    return scipy.stats.norm.ppf((x.rank() - 0.5) / x.count())


def findTwoComponents(df, meta):
    df = df.T
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(df)
    pDf = pd.DataFrame(data=principalComponents
                       , columns=['PC1', 'PC2'])
    pDf.set_index(df.index, inplace=True)
    pc1_var = round(pca.explained_variance_ratio_[0] * 100, 2)
    pc2_var = round(pca.explained_variance_ratio_[1] * 100, 2)
    pDf2 = pDf.merge(meta, left_index=True, right_index=True)
    return pDf2, pc1_var, pc2_var


def plotPCA(pDf, pc1_var, pc2_var, colorby, c="", nameby="", title="", filename='',
            el=False):  # , xlimits, ylimits, labels): # colorby is column in pDf
    sns.set_style("ticks")
    group = pDf[colorby].unique()

    assert len(group) < 5
    if c:
        colrs = c
    else:
        colrs = colors[:len(group) + 1]

    fig = plt.figure(figsize=(8, 8))
    for g, c in zip(group, colrs):
        df = pDf[pDf[colorby] == g]
        x, y = df[["PC1"]].values, df[["PC2"]].values
        pts = np.asarray([[float(a), float(b)] for a, b in zip(x, y)])
        ax = plt.scatter(x, y, c=c, s=150, label=g)
        if el:
            plot_point_cov(pts, nstd=2, alpha=0.1, color=c)
        # ax = plt.scatter(df.PC1, df.PC2, c= c, s = 150,label = g)
        plt.legend(fontsize=18, frameon=True)
        if nameby:
            labels = df[nameby]
            for label, pc1, pc2 in zip(labels, df.PC1, df.PC2):
                plt.annotate(label, xy=(pc1, pc2), xytext=(-5, 7), textcoords="offset points",
                             fontsize=14)
        plt.xlabel('Principal Component 1, {} %'.format(pc1_var), fontsize=15)
        plt.ylabel('Principal Component 2, {} %'.format(pc2_var), fontsize=15)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        if title:
            plt.title('2 component PCA', fontsize=20)
    if filename:
        fig.savefig(filename, dpi=300)
    return fig



def plot_point_cov(points, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)

def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


def remove_hypotheticals(df, col_name="Function"):
    keep = []
    for i, x in zip(df.index, df[col_name]):
        if not type(x) == str:
            continue
        if not "hypothetical" in x and not "phage" in x:
            keep.append(i)
    return df.loc[keep]

def get_subset_genes(df, key, col_return,column_name="function"):
    keep = []
    for i, x in zip(df.index, df[column_name]):
        if not type(x) == str:
            continue
        if key in x:
            keep.append(i)
    return df.loc[keep][col_return]


def draw_heatmap_of_subset(df_index, meta, title, rpkms, samples,
                           my_cmap, fs=(4,4), draw=True):
    """
    if draw is True returns figure, if False returns dataframe

    """
    dline_meta = meta[meta["group.ID"].isin(["Case7", "Case8", "Case11", "Case12"])]
    subset_rpkms = dline_meta.merge(rpkms.loc[df_index][dline_meta.index].T, left_index=True, right_index=True)
    subset_means = {}
    for gene in df_index:
        subset_means[gene] = {}
        for case in subset_rpkms["group.ID"].unique():
            m = round(subset_rpkms[subset_rpkms["group.ID"] == case][gene].mean(), 2)
            subset_means[gene][case] = m
    t = pd.DataFrame(subset_means).T
    t = t[["Case12", "Case7", "Case11", "Case8"]]
    t.rename(index=str, columns={c: samples[c] for c in t.columns}, inplace=True)
    if draw:
        fig = plt.figure(figsize=fs)
        s = sns.heatmap(np.log2(t + 1), cmap=my_cmap,
                        linewidths=0.5, linecolor='black',
                        cbar_kws={'label': 'Log2 RPKMs'});

        s.set_title(title)
        return fig
    return t


def join_subset_means(label_subset_dict, meta, rpkms, samples, my_cmap):
    means_list = []
    for label, subset in label_subset_dict.items():
        means_list.append(find_subset_mean(label, subset, meta, rpkms, samples, my_cmap))

    return pd.concat(means_list, axis=1, keys=[s.name for s in means_list]).T


def find_subset_mean(label, subset, meta, rpkms, samples, my_cmap):
    cts = draw_heatmap_of_subset(subset.index, meta, label,
                                 rpkms, samples, my_cmap, draw=False).mean()
    cts.name = label
    return cts