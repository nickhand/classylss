import numpy

class HalofitPower(object):

    def __init__(self, cosmo, z):

        # make sure nonlinear is enabled in CLASS
        if not cosmo.nonlinear:
            cosmo = cosmo.clone(**{'non linear':'halofit'})
        self.cosmo = cosmo
        self.z = z

    def __call__(self, k):

        k = numpy.asarray(k)
        kmin = self.cosmo.P_k_min
        inrange = k > 1.00001*kmin

        Pk = numpy.zeros_like(k)
        k_in = k[inrange]; k_out = k[~inrange]

        Pk[inrange] = self.cosmo.get_pk(k=k_in, z=self.z)
        Pk[~inrange] = self.cosmo.get_pklin(k=k_out, z=self.z)
        return Pk
