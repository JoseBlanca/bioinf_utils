from sklearn.preprocessing import power_transform, quantile_transform
from pandas import Series, DataFrame


def _normalize_series_boxcox(series, method):
    normed_values = power_transform(series.values.reshape(-1, 1), method=method)
    return Series(normed_values.flat, index=series.index)


def normalize_series(series, method="yeo-johnson"):
    if method in ("box-cox", "yeo-johnson"):
        return _normalize_series_boxcox(series, method=method)


def normalize_dframe(dframe, method):
    return DataFrame(
        power_transform(dframe, method=method),
        index=dframe.index,
        columns=dframe.columns,
    )


def quantile_transform_series(series, n_quantiles):
    return quantile_transform(series.values.reshape(-1, 1), n_quantiles=n_quantiles)


def quantile_transform_dframe(dframe, n_quantiles):
    return DataFrame(
        quantile_transform(dframe, n_quantiles=n_quantiles),
        index=dframe.index,
        columns=dframe.columns,
    )
