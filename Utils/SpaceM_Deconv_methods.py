"""
Implementations of deconvolution methods

See :doc:`/user-guide/background/deconvolution-methods`
"""

from collections.abc import Callable
from typing import Final
from functools import singledispatch
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix

import numpy as np
from typing import Any, TYPE_CHECKING, Literal, Protocol, TypeVar, Union
import xarray as xr
import dask
import dask.array as da



# AblationMarkValuesListType = np.ndarray[Shape[NumAblationMarks], float]
# CellValuesListType = np.ndarray[Shape[NumCells], float]
# OverlapMatrixType = XDataArrayType[Shape2d[NumAblationMarks, NumCells], float]



A = TypeVar("A", bound=Union[np.ndarray, xr.DataArray, da.Array, pd.DataFrame])

@singledispatch
def fillna(array: A, value: Any) -> A:
    """
    Fill NaN in different array implementations

    This is a workaround for that Numpy and Dask have ``nan_to_num`` while others have ``fillna``.
    :py:func:`numpy.nan_to_num` accepts other arrays but always returns a dense Numpy array.
    :py:meth:`xarray.DataArray.fillna` accepts other arrays (except DataFrame) and returns the
    same type of array.

    Args:
        array: An array
        value: The value to fill in for NaN

    Returns:
        An array of the same type as the input
    """
    raise NotImplementedError


@fillna.register(np.ndarray)
@fillna.register(xr.DataArray)
def _(array, value):
    # Unlike np.nan_to_num, xarray.DataArray preserves the input array type, except for pd.DataFrame
    result = xr.DataArray.fillna(array, value)
    assert isinstance(result, type(array))
    return result


@fillna.register
def _(array: pd.DataFrame, value) -> pd.DataFrame:
    # Only pandas.DataFrame.fillna returns again a DataFrame
    return array.fillna(value=value)

@singledispatch
def ensure_dense(
    array: Union[
        np.ndarray, xr.DataArray, da.Array, csc_matrix, csr_matrix, pd.DataFrame
    ],
):
    """
    Make sparse arrays dense

    Conversion to dense can be needed for functions that don't support sparse arrays or a specific
    sparse implementation. This function supports sparse array implementations of
    :py:mod:`scipy.sparse` and :py:mod:`sparse`.

    Args:
        array: An array

    Returns:
        A dense array. If the input array is already dense, it is returned unchanged. If the input
        array wraps a sparse array, a new instance of the same type wrapping a dense array is
        returned. If the type cannot be preserved, a Numpy array is returned.
    """
    raise NotImplementedError


@ensure_dense.register
def _(array: np.ndarray) -> np.ndarray:
    return array


@ensure_dense.register
def _(array: xr.DataArray) -> xr.DataArray:
    # Workaround for that the todense method is not directly on the DataArray object,
    # but nested on the sparse data object.
    if isinstance(array.data, (csc_matrix, csr_matrix, da.Array)):
        dense_array = array.copy()
        dense_array.data = ensure_dense(array.data)
        return dense_array
    return array


@ensure_dense.register
def _(array: da.Array) -> da.Array:
    return da.from_delayed(
        dask.delayed(lambda a: (ensure_dense(a),), nout=1)(array),
        shape=array.shape,
        dtype=array.dtype,
    )


@ensure_dense.register(csc_matrix)
@ensure_dense.register(csr_matrix)
def _(array: Union[csc_matrix, csr_matrix]) -> np.ndarray:
    return array.todense()  # noqa (sparse.SparseArray -> np.ndarray)


@ensure_dense.register
def _(array: pd.DataFrame) -> pd.DataFrame:
    return array.sparse.to_dense()




def inverse_sampling_proportion_weighted_by_sampling_specificity(
    overlap_matrix,
    am_areas,
    cell_areas
):
    r"""
    This method is close to the one published in the paper

        [The intensity assigned to a cell for a given metabolite] is calculated as a weighted mean
        of the metabolite intensities from the ablation marks sampling that cell. The intensities
        are divided by the sampling proportion to account for differences in amount of sampled
        cellular material between ablation marks. To increase the contribution of ablation marks
        which sample the cell of interest more than other ablation marks, their intensities are
        weighted by the sampling specificity.

        .. image:: /user-guide/background/images/deconvolution-spacem-paper.png
           :align: center

    - **Filtering** (separate step): mostly off-cell ablation marks
    - **Correction**: mostly on-cell ablation marks.
      Assuming all (on-cell) ions come from the cell, the ablation mark should have a higher
      ion count if it was fully overlapping the cell.
    - **Weighting**: by sampling specificity

    .. math::
       sampling\;proportion(a_n) = \frac{\sum_{c_i \in cells}{area(a_n \cap c_i)}}{area(a_n)}

    .. math::
       sampling\;specificity(a_n, c_k) = \frac{area(a_n \cap c_k)}{\sum_{c_i \in cells}{area(a_n \cap c_i)}}

    .. math::
       intensity(c_k) = \frac{
           \sum_{a_n \in a.marks}{
               intensity(a_n) \cdot \dfrac{1}{sampling\;proportion(a_n)} \cdot sampling\;specificity(a_n, c_k)
           }
       }{
           \sum_{a_n \in a.marks}{sampling\;specificity(a_n, c_k)}
       }
    """
    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    am_sampling_proportion = am_intracellular_sampling_areas / am_areas
    sampling_specificity = (overlap_matrix.T / am_intracellular_sampling_areas).T
    coefficients_matrix = ((sampling_specificity * 1.0).T / am_sampling_proportion).T / np.sum(
        fillna(sampling_specificity, value=0.0), axis=0
        )
    return fillna(coefficients_matrix, value=0.0)

def inverse_sampling_proportion_weighted_by_sampling_specificity_corrected(
    overlap_matrix,
    am_areas,
    cell_areas
):

    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    am_sampling_proportion = am_intracellular_sampling_areas / am_areas
    sampling_specificity = (overlap_matrix.T / am_intracellular_sampling_areas).T
    coefficients_matrix = ((sampling_specificity * 1.0).T * am_sampling_proportion).T / np.sum(
        fillna(sampling_specificity, value=0.0), axis=0
        )
    return fillna(coefficients_matrix, value=0.0)

def biggest_overlap_raw(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    No deconvolution, select for each cell the ablation mark with the largest overlap

    This is the most basic method without normalization. For a cell hit by one or more ablation
    marks, only the one ablation mark with maximum overlap area is selected and all its signal is
    assigned directly to the cell, without proportional scaling.

    This is a baseline avoiding the complexity of interpolation, under the ideal assumption that
    each cell is sampled by exactly one ablation mark with full overlap.

    .. math::
       intensity(c_k) = intensity \left(
           \underset{a_n \in a.marks}{\operatorname{arg\,max}}\, \bigl(
               area(a_n \cap c_k)
           \bigr)
       \right)

    .. image:: /user-guide/background/images/deconvolution_biggest_overlap_raw.svg
       :align: center
    """
    # Note: xarray max does not support argument initial=0.0
    max_overlap_per_cell = np.max(overlap_matrix, axis=0, keepdims=True)
    # This is also true if a cell is not overlapped at all (max = 0), therefore the second condition
    cell_and_am_overlap_is_max = overlap_matrix == np.broadcast_to(
        ensure_dense(max_overlap_per_cell), overlap_matrix.shape
    )
    cell_and_am_have_overlap = overlap_matrix > 0
    coefficients_matrix = cell_and_am_have_overlap & cell_and_am_overlap_is_max
    return fillna(coefficients_matrix, value=0.0)


def mean_sampling_proportion(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    Weighted by mean sampling proportion (original)

    .. math::
       sampling\;proportion(a_n) = \frac{\sum_{c_i \in cells}{area(a_n \cap c_i)}}{area(a_n)}

    .. math::
       \delta_{area(a_n \cap c_k) > 0} = \begin{cases}
       1, & \text{if } area(a_n \cap c_k) > 0\\
       0, & \text{otherwise}
       \end{cases}

    .. math::
       intensity(c_k) = \sum_{a_n \in a.marks}{
           intensity(a_n) \cdot \delta_{area(a_n \cap c_k) > 0} \cdot \dfrac{ 1 }{ sampling\;proportion(a_n) }
       }

    .. deprecated:: This method is included from the original prototype for completeness,
       but has not been successfully validated
    """


    # Ablation mark sampling area is the total area of an ablation mark that overlaps cells
    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    am_sampling_proportion = am_intracellular_sampling_areas / am_areas
    cell_and_am_that_have_overlap = overlap_matrix > 0
    # Temporarily setting sparse fill_value zero to NaN ensures that the next operation works,
    # where otherwise uniqueness of fill value would not be preserved.
    # (fill_value = 0; [fill_value, fill_value] / [1, 0] = [0, NaN])
    cell_and_am_that_have_overlap = xr.where(cell_and_am_that_have_overlap, 1.0, np.nan)
    # For every cell consider all its overlapping ablation marks equally only weighted by their
    # sampling ratio (not specific for each cell).
    coefficients_matrix = (cell_and_am_that_have_overlap.T / am_sampling_proportion).T
    return fillna(coefficients_matrix, value=0.0)


def mean_sampling_proportion_corrected(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    Weighted by mean sampling proportion (inverted weight)

    .. math::
       intensity(c_k) = \sum_{a_n \in a.marks}{
           intensity(a_n) \cdot \delta_{area(a_n \cap c_k) > 0} \cdot sampling\;proportion(a_n)
       }

    .. deprecated:: This method is included from the original prototype for completeness,
       but has not been successfully validated
    """

    # Ablation mark sampling area is the total area of an ablation mark that overlaps cells
    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    am_sampling_proportion = am_intracellular_sampling_areas / am_areas
    cell_and_am_that_have_overlap = overlap_matrix > 0
    # Temporarily setting sparse fill_value zero to NaN ensures that the next operation works,
    # where otherwise uniqueness of fill value would not be preserved.
    cell_and_am_that_have_overlap = xr.where(
        cell_and_am_that_have_overlap, cell_and_am_that_have_overlap, np.nan
    )
    # For every cell consider all its overlapping ablation marks equally only weighted by their
    # sampling ratio (not specific for each cell).
    coefficients_matrix = (cell_and_am_that_have_overlap.T * am_sampling_proportion).T
    return fillna(coefficients_matrix, value=0.0)

def mean_sampling_area(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    Weighted by mean samping area

    .. math::
       ablation\;mark\;sampling\;area(a_n) = \sum_{c_i \in cells}{area(a_n \cap c_i)}

    .. math::
       intensity(c_k) = \sum_{a_n \in a.marks}{
           intensity(a_n) \cdot \delta_{area(a_n \cap c_k) > 0} \dfrac{ 1 }{ \cdot ablation\;mark\;sampling\;area(a_n) }
       }

    .. deprecated:: This method is included from the original prototype for completeness,
       but has not been successfully validated
    """
    # Ablation mark sampling area is the total area of an ablation mark that overlaps cells
    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    cell_and_am_that_have_overlap = overlap_matrix > 0
    # For every cell consider all its overlapping ablation marks equally only weighted by their
    # sampling area (not specific for each cell).
    coefficients_matrix = (cell_and_am_that_have_overlap.T / am_intracellular_sampling_areas).T
    return fillna(coefficients_matrix, value=0.0)


def mass_density_weighted_by_overlap_areas(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    Ion density weighted by sampling area of the ablation mark and sampled area of the cell

    .. math::
       cell\;sampling\;area(c_k) = \sum_{a_i \in cells}{area(a_i \cap c_k)}

    .. math::
        intensity(c_k) = \sum_{a_n \in a.marks}{
            intensity(a_n) \cdot \dfrac{ area(a_n \cap c_k) }{ area(a_n) } \cdot \dfrac{ area(a_n \cap c_k) }{ cell\;sampling\;area(c_k) }
        }

    .. deprecated:: This method is included from the original prototype for completeness,
       but has not been successfully validated
    """

    cell_sampled_areas = np.sum(fillna(overlap_matrix, value=0.0), axis=0)
    # Temporarily setting sparse fill_value zero to NaN ensures that the next operation works,
    # where otherwise uniqueness of fill value would not be preserved.
    overlap_matrix = xr.where(overlap_matrix != 0, overlap_matrix, np.nan)
    coefficients_matrix = (overlap_matrix.T / am_areas).T * (overlap_matrix / cell_sampled_areas)
    return fillna(coefficients_matrix, value=0.0)


def weighted_mean_sampling_area_mark_cell_overlap_int(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    Sampling specificity weighted by sampled area of the cell

    .. math::
        intensity(c_k) = \sum_{a_n \in a.marks}{
            intensity(a_n) \cdot \dfrac{ sampling\;specificity(a_n, c_k) }{ cell\;sampling\;area(c_k) }
        }

    .. deprecated:: This method is included from the original prototype for completeness,
       but has not been successfully validated
    """
    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    cell_sampled_areas = np.sum(overlap_matrix, axis=0)
    sampling_specificity = (overlap_matrix.T / am_intracellular_sampling_areas).T
    coefficients_matrix = sampling_specificity / cell_sampled_areas
    return fillna(coefficients_matrix, value=0.0)


def weighted_by_overlap_area(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    r"""
    Extrapolate ion intensities to the total cell area, assuming uniform density

    The ion densities of the ablation marks are averaged, weighted by the proportion of each overlap
    relative to the total cell sampling area (total of ablation mark overlaps with that cell).
    The average ion density is scaled proportionally to the whole cell area.

    .. math::
       intensity(c_k) = \frac{
           \sum_{a_n \in a.marks}{
               intensity(a_n) \cdot \frac{area(a_n \cap c_k)}{area(a_n)}
           }
       }{
           \sum_{a_n \in a.marks}{ \frac{ area(a_n \cap c_k) }{ area(c_k) } }
       }

    .. image:: /user-guide/background/images/deconvolution_weighted_by_sampling_proportion.svg
       :align: center
    """

    # The ratio of all sampled cell area relative to the cell.
    cell_sampled_areas = np.sum(fillna(overlap_matrix, value=0.0), axis=0)
    # Temporarily setting sparse fill_value zero to NaN ensures that the next operation works,
    # where otherwise uniqueness of fill value would not be preserved.
    overlap_matrix = xr.where(overlap_matrix != 0, overlap_matrix, np.nan)
    # The ratio of each overlap relative to the ablation mark.
    # This is the fraction of the ion count in the overlap area.
    am_overlap_ratio = (overlap_matrix.T / am_areas).T
    cell_sampled_proportion = cell_sampled_areas / np.asarray(cell_areas)
    # The factor that each ablation mark contributes to each cell.
    coefficients_matrix = am_overlap_ratio / cell_sampled_proportion
    return fillna(coefficients_matrix, value=0.0)


def inverse_sampling_proportion_weighted_by_overlap_area(
    overlap_matrix,
    am_areas,
    cell_areas,
):
    """
    This method is inspired by to the one published in the paper

    - **Filtering** (separate step): mostly off-cell ablation marks
    - **Correction**: mostly on-cell ablation marks.
      Assuming all (on-cell) ions come from the cell, the ablation mark should have a higher
      ion count if it was fully overlapping the cell.
    - **Weighting**: by overlap.
      This method does not weight by sampling specificity but is strictly proportional by area.
    """

    # The ratio of all sampled cell area relative to the cell.
    cell_sampled_areas = np.sum(fillna(overlap_matrix, value=0.0), axis=0)
    # The ratio of all sampling ablation mark area relative to the whole ablation mark
    am_intracellular_sampling_areas = np.sum(overlap_matrix, axis=1)
    am_sampling_proportion = am_intracellular_sampling_areas / am_areas
    # Temporarily setting sparse fill_value zero to NaN ensures that the next operation works,
    # where otherwise uniqueness of fill value would not be preserved.
    overlap_matrix = xr.where(overlap_matrix != 0, overlap_matrix, np.nan)
    # The ratio of each overlap relative to the ablation mark.
    # This is the fraction of the ion count in the overlap area.
    am_overlap_ratio = (overlap_matrix.T / am_areas).T
    cell_sampled_proportion = cell_sampled_areas / np.asarray(cell_areas)
    # The factor that each ablation mark contributes to each cell.
    coefficients_matrix = (
        (am_overlap_ratio / cell_sampled_proportion).T * (1.0 / am_sampling_proportion)
    ).T
    return fillna(coefficients_matrix, value=0.0)


DECONVOLUTION_METHODS = {
    "inverse_sampling_proportion_weighted_by_sampling_specificity": inverse_sampling_proportion_weighted_by_sampling_specificity,
    "inverse_sampling_proportion_weighted_by_sampling_specificity_corrected" : inverse_sampling_proportion_weighted_by_sampling_specificity_corrected,
    "biggest_overlap_raw": biggest_overlap_raw,
    "mean_sampling_proportion": mean_sampling_proportion,
    "mean_sampling_proportion_corrected": mean_sampling_proportion_corrected,
    "mean_sampling_area": mean_sampling_area,
    "mass_density_weighted_by_overlap_areas": mass_density_weighted_by_overlap_areas,
    "weighted_mean_sampling_area_mark_cell_overlap_int": weighted_mean_sampling_area_mark_cell_overlap_int,
    "weighted_by_overlap_area": weighted_by_overlap_area,
    "inverse_sampling_proportion_weighted_by_overlap_area": inverse_sampling_proportion_weighted_by_overlap_area,
}