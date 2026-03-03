def test_import():
    from ..FileOrganizer import FileOrganizer
    assert True

from ..FileOrganizer import FileOrganizer

def test_init():
    pattern_map = {'SL': r'^(?=.*(?:2B2JLNu|2B2Q1L1Nu)).*M(?:X)?_(\d+)',
                   'DL': r'^(?=.*2B2L2Nu).*M(?:X)?_(\d+)'}
    for _, pattern in pattern_map.items():
        _ = FileOrganizer(start_dir='../train_data', pattern=pattern)
    assert True