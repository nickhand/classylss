from classylss.cosmology import Cosmology

def test_cosmology_init():
    c = Cosmology()
    c.get_transfer(z=0)
    dir(c)
