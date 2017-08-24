from classylss.correlation import CorrelationFunction
from classylss.cosmology import Cosmology
from classylss.power import LinearPower, HalofitPower, ZeldovichPower
import numpy

def test_linear():

    # linear power
    Plin = LinearPower(Cosmology(), z=0)

    # desired separation (in Mpc/h)
    r = numpy.logspace(0, numpy.log10(150), 500)

    # linear correlation
    CF = CorrelationFunction(Plin)

    xi = CF(r)

def test_halofit():

    # nonlinear power
    Pnl = HalofitPower(Cosmology(), z=0)

    # desired separation (in Mpc/h)
    r = numpy.logspace(0, numpy.log10(150), 500)

    # nonlinear correlation
    CF = CorrelationFunction(Pnl)

    xi = CF(r)

def test_zeldovich():

    # zeldovich power
    Pzel = ZeldovichPower(Cosmology(), z=0)

    # desired separation (in Mpc/h)
    r = numpy.logspace(0, numpy.log10(150), 500)

    # zeldovich correlation
    CF = CorrelationFunction(Pzel)

    xi = CF(r)
