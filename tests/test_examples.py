from pathlib import Path
import sys
import numpy as np

this_file_loc = Path(__file__).parent.resolve()
examples_loc = this_file_loc.parent / 'examples'
sys.path.insert(0, str(examples_loc))

def test_examples_run():
    """
    simply imports all examples: since they are scripts, they execute on import.
    The variables of each example are available under the import name space
    e.g. everything created in _2_marker_matching is under _2_marker_matching.relevant_variable
    """

    from examples import _2_marker_matching
    from examples import _3_field_calculation
    from examples import _4_fit_harmonics
    from examples import _5_reporting
    from examples import _6_all_together_now

    assert np.allclose(_5_reporting.report._MatchedMarkerVolume.abs_dis, _6_all_together_now.report._MatchedMarkerVolume.abs_dis)
    assert np.mean(_5_reporting.report._MatchedMarkerVolume.abs_dis) < 10  # was 8.73, could change a bit with changes to other algorithms
    assert np.max(_5_reporting.report._MatchedMarkerVolume.abs_dis) < 100  # was 71

