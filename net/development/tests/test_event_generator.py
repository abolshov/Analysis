def test_import():
    from ..EventGenerator import EventGenerator
    assert True

from ..EventGenerator import EventGenerator

def test_init():
    generator = EventGenerator(masspoint=800,
                               file_list=[],
                               branches_to_load=[],
                               features=[],
                               targets=[],
                               step_size="50 MB")
    assert True