import awkward as ak
from ..Utils import mask_indices

def test_mask_indices():
    data = ak.Array([[10, 20, 30], [40, 50]])
    result = ak.Array([[True, False, True], [True, False]])
    indices = [[0, 2]]
    assert ak.all(mask_indices(indices, data) == result)