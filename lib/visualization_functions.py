import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from settings import *
from sklearn.decomposition import PCA


def invnorm(x):
    return scipy.stats.norm.ppf((x.rank() - 0.5) / x.count())


def findTwoComponents(df, meta):
    df = df.T
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(df)
    pca_df = pd.DataFrame(data=principal_components,
                          columns=['PC1', 'PC2'])
    pca_df.set_index(df.index, inplace=True)
    pc1_var = round(pca.explained_variance_ratio_[0] * 100, 2)
    pc2_var = round(pca.explained_variance_ratio_[1] * 100, 2)
    pca_df2 = pca_df.merge(meta, left_index=True, right_index=True)
    return pca_df2, pc1_var, pc2_var


def plotPCA(pca_df, pc1_var, pc2_var, colorby, clrs, nameby="", title="", filename='',
            el=False):  # , xlimits, ylimits, labels): # colorby is column in pDf
    sns.set_style("ticks")
    group = pca_df[colorby].unique()
    fig = plt.figure(figsize=(8, 8))
    for g, c in zip(group, clrs):
        df = pca_df[pca_df[colorby] == g]
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
    author:
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
    author:
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


def draw_heatmap_of_subset(cnts, info, genes,  names, col_id='case',
                           cases= ["Case12", "Case7", "Case11", "Case8"], subset_name="test",
                           draw=False, fs=(4,4), my_cmap=''):
    """
    given a list of genes and cases calculate mean for each gene + for whole subset
    """
    assert col_id in info.columns
    info = info[info[col_id].isin(cases)]
    data = info.merge(cnts.T, left_index=True, right_index=True)
    subset_means = []
    gene_means = []
    groups = data.groupby(col_id)
    for case in cases:
        df = groups.get_group(case)
        gene_means.append(pd.Series(df[genes].mean(), name=names[case]))
        subset_means.append((names[case], df[genes].mean().mean()))
    subset_df = pd.DataFrame.from_records(subset_means, columns=["Sample", subset_name], index="Sample")
    subset_df.index.name=""
    gene_df = pd.concat(gene_means, axis=1)
    if draw:
        fig = plt.figure(figsize=fs)
        s = sns.heatmap(np.log2(gene_df + 1), cmap=my_cmap, linewidths=0.5, linecolor='black',
                        cbar_kws={'label': 'Log2 TPMs'})
        return fig, ''
    return gene_df, subset_df


def join_subset_means(subset_dict, cnts, info, cases, names, draw=False, fs=(4,4), my_cmap=''):
    subsets=[]
    for label, subset in subset_dict.items():
        _, subset_df = draw_heatmap_of_subset(cnts, info, subset, cases=cases,
                                              names=names, subset_name=label)
        subsets.append(subset_df)
    subset_df = pd.concat(subsets, axis=1).T
    if draw:
        fig = plt.figure(figsize=fs)
        sns.heatmap(np.log2(subset_df + 1), cmap=my_cmap, linewidths=0.5, linecolor='black',
                        cbar_kws={'label': 'Log2 TPMs'})
        return fig
    return subset_df

#
# def draw_heatmap_of_subset(genes, meta, title, cnts, samples,
#                            my_cmap, fs=(4, 4), col = 'case',
#                            cases=["Case12", "Case7", "Case11", "Case8"], draw=True):
#     """
#     if draw is True returns figure, if False returns dataframe
#
#     """
#     case_meta = meta[meta[col].isin(cases)]
#     samples = case_meta.index
#     case_cnts = case_meta.join(cnts.loc[genes].T, how='inner')
#     subset_means = {}
#     for gene in genes:
#         subset_means[gene] = {}
#         for case in case_cnts[col].unique():
#             m = round(case_cnts[case_cnts[col] == case][gene].mean(), 2)
#             subset_means[gene][case] = m
#     t = pd.DataFrame(subset_means).T
#     t = t[cases]
#     t.rename(index=str, columns={c: samples[c] for c in t.columns}, inplace=True)
#     if draw:
#         fig = plt.figure(figsize=fs)
#         s = sns.heatmap(np.log2(t + 1), cmap=my_cmap, linewidths=0.5, linecolor='black',
#                         cbar_kws={'label': 'Log2 TPMs'})
#         s.set_title(title)
#         return fig
#     return t

# Should not need these functions below


