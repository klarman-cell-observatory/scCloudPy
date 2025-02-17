import time
import numpy as np
from scipy.sparse import issparse
from anndata import AnnData
import logging

logger = logging.getLogger("sccloud")

from sccloud.tools import estimate_feature_statistics, select_features


def set_group_attribute(data: AnnData, attribute_string: str) -> None:
    """Set group attributes used in batch correction.

    Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, *sccloud* will only perform batch correction for channels within each group.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    attribute_string: ``str``
        Attributes used to construct groups:

        * ``None``
            Assume all channels are from one group.

        * ``attr`` 
            Define groups by sample attribute ``attr``, which is a keyword in ``data.obs``.

        * ``att1+att2+...+attrn`` 
            Define groups by the Cartesian product of these *n* attributes, which are keywords in ``data.obs``.

        * ``attr=value_11,...value_1n_1;value_21,...value_2n_2;...;value_m1,...,value_mn_m``
            In this form, there will be *(m+1)* groups. A cell belongs to group *i* (*i > 1*) if and only if its sample attribute ``attr``, which is a keyword in ``data.obs``, has a value among ``value_i1``, ... ``value_in_i``. A cell belongs to group 0 if it does not belong to any other groups.

    Returns
    -------
    None

        Update ``data.obs``:
        
        * ``data.obs["Group"]``: Group ID for each cell.

    Examples
    --------

    >>> scc.set_group_attribute(adata, attr_string = "Individual")

    >>> scc.set_group_attribute(adata, attr_string = "Individual+assignment")

    >>> scc.set_group_attribute(adata, attr_string = "Channel=1,3,5;2,4,6,8")    
    """
    
    if attribute_string.find("=") >= 0:
        attr, value_str = attribute_string.split("=")
        assert attr in data.obs.columns
        values = value_str.split(";")
        data.obs["Group"] = "0"
        for group_id, value in enumerate(values):
            vals = value.split(",")
            idx = np.isin(data.obs[attr], vals)
            data.obs.loc[idx, "Group"] = str(group_id + 1)
    elif attribute_string.find("+") >= 0:
        attrs = attribute_string.split("+")
        assert np.isin(attrs, data.obs.columns).sum() == len(attrs)
        data.obs["Group"] = data.obs[attrs].apply(lambda x: "+".join(x), axis=1)
    else:
        assert attribute_string in data.obs.columns
        data.obs["Group"] = data.obs[attribute_string]


def estimate_adjustment_matrices(data: AnnData) -> bool:
    """ Estimate adjustment matrices
    """

    if ("gmeans" not in data.varm) or ("gstds" not in data.varm):
        estimate_feature_statistics(data, True)

    if data.uns["Channels"].size == 1:
        logger.warning(
            "Warning: data only contains 1 channel. Batch correction disabled!"
        )
        return False

    nchannel = data.uns["Channels"].size

    plus = np.zeros((data.shape[1], nchannel))
    muls = np.zeros((data.shape[1], nchannel))

    ncells = data.uns["ncells"]
    means = data.varm["means"]
    partial_sum = data.varm["partial_sum"]
    gmeans = data.varm["gmeans"]
    gstds = data.varm["gstds"]
    c2gid = data.uns["c2gid"]
    for i in range(data.uns["Channels"].size):
        if ncells[i] > 1:
            muls[:, i] = (partial_sum[:, i] / (ncells[i] - 1.0)) ** 0.5
        outliers = muls[:, i] < 1e-6
        normals = np.logical_not(outliers)
        muls[outliers, i] = 1.0
        muls[normals, i] = gstds[normals, c2gid[i]] / muls[normals, i]
        plus[:, i] = gmeans[:, c2gid[i]] - muls[:, i] * means[:, i]

    data.varm["plus"] = plus
    data.varm["muls"] = muls

    return True


def correct_batch_effects(data: AnnData, keyword: str, features: str = None) -> None:
    """ Apply calculated plus and muls to correct batch effects for a dense matrix
    """
    X = data.uns[keyword]
    m = X.shape[1]
    if features is not None:
        selected = data.var[features].values
        plus = data.varm["plus"][selected, :]
        muls = data.varm["muls"][selected, :]
    else:
        selected = np.ones(data.shape[1], dtype=bool)
        plus = data.varm["plus"]
        muls = data.varm["muls"]

    for i, channel in enumerate(data.uns["Channels"]):
        idx = np.isin(data.obs["Channel"], channel)
        if idx.sum() == 0:
            continue
        X[idx] = X[idx] * np.reshape(muls[:, i], newshape=(1, m)) + np.reshape(
            plus[:, i], newshape=(1, m)
        )
    # X[X < 0.0] = 0.0


def correct_batch(data: AnnData, features: str = None) -> None:
    """Batch correction on data.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    features: `str`, optional, default: ``None``
        Features to be included in batch correction computation. If ``None``, simply consider all features.

    Returns
    -------
    ``None``

    Update ``data.X`` by the corrected count matrix.

    Examples
    --------
    >>> scc.correct_batch(adata, features = "highly_variable_features")
    """

    tot_seconds = 0.0

    # estimate adjustment parameters
    start = time.time()
    can_correct = estimate_adjustment_matrices(data)
    end = time.time()
    tot_seconds += end - start
    logger.info("Adjustment parameters are estimated.")

    # select dense matrix
    keyword = select_features(data, features)
    logger.info("Features are selected.")

    if can_correct:
        start = time.time()
        correct_batch_effects(data, keyword, features)
        end = time.time()
        tot_seconds += end - start
        logger.info(
            "Batch correction is finished. Time spent = {:.2f}s.".format(tot_seconds)
        )
