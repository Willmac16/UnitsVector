import numpy as np
from Two_Dimensional_Fit.twoDimFit import rref, solve

# Widest Class: should be broader than MKS/SI units
class UnitsVector:
    def __init__(self, vec, val, val_scale):
        self.value = val
        self.vector = vec
        self.value_scale = val_scale

    def __add__(self, other: UnitsVector):
        if np.equal(self.vector, other.vector).all():
            abs_value = self.value * self.value_scale + other.value * other.value_scale
            relative_val = abs_value / self.value_scale

            return UnitsVector(self.vector, relative_val, self.value_scale)

        else:
            raise Exception("UnitsVector addition error: units do not match")

    def __sub__(self, other: UnitsVector):
        if np.equal(self.vector, other.vector).all():
            abs_value = self.value * self.value_scale - other.value * other.value_scale
            relative_val = abs_value / self.value_scale

            return UnitsVector(self.vector, relative_val, self.value_scale)

        else:
            raise Exception("UnitsVector subtraction error: units do not match")

    def __mul__(self, other: UnitsVector):
        out_vec = self.vector + other.vector
        out_val = self.value * other.value * self.value_scale * other.value_scale

        return UnitsVector(out_vec, out_val/self.value_scale, self.value_scale)

    def __truediv__(self, other: UnitsVector):
        out_vec = self.vector - other.vector
        out_val = (self.value * self.value_scale) / (other.value * other.value_scale)

        return UnitsVector(out_vec, out_val / self.value_scale, self.value_scale)

    # I plan on printing in mks unless explicitly casted to another unit system
    def __str__(self):
        out_str = str(self.value * self.value_scale)
        out_str += " (s^{0[0]} m^{0[1]} kg^{0[2]} A^{0[3]} K^{0[4]} mol^{0[5]} cd^{0[6]})"
        return out_str.format(self.vector)



class MKS:
    def __init__(self, value, t, l, m, j, k, mol, cd):
        self.value = value
        self.value_scale = 1

        self.vector = np.array([t, l, m, j, k, mol, cd])

        # self.seconds = t
        # self.meters = l
        # self.kilograms = m
        # self.amps = j
        # self.kelvins = k
        # self.moles = mol
        # self.candela = cd

# Imperial stuff
seconds_to_seconds = 1
feet_to_meters = 12 * 0.0254
pound_mass_to_kilograms = 0.45359237
amps_to_amps = 1
rankine_to_kelvins = 5 / 9
pound_moles_to_moles = 453.59237
candela_to_candela = 1

imperial_to_MKS_matrix = np.identity(7)
imperial_scale_vector = np.array(seconds_to_seconds, feet_to_meters, pound_mass_to_kilograms, amps_to_amps, rankine_to_kelvins, pound_moles_to_moles, candela_to_candela)


class Imperial:
    def __init__(self, value, t, l, lbm, j, R, mol, cd):
        self.value = value

        base_vector = np.array([t, l, lbm, j, R, lb-mol, cd])


        # self.seconds = t
        # self.feet = l
        # self.pound-mass = lbm
        # self.amps = j
        # self.rankine = R
        # self.moles = mol
        # self.candela = cd

        self.vector = base_vector
        self.value_scale = imperial_scale_vector @ base_vector

# Imperial to MKS and Back
def toImperial(value: UnitsVector):
    out_value = value * value.value_scale / (imperial_scale_vector @ value.vector)
    out_vector = value.vector

    return Imperial(out_value, out_vector[0], out_vector[1], out_vector[2], out_vector[3], out_vector[4], out_vector[5], out_vector[6])

def toMKS(value: UnitsVector):
    out_value = value * value.value_scale
    out_vector = value.vector

    return MKS(out_value, out_vector[0], out_vector[1], out_vector[2], out_vector[3], out_vector[4], out_vector[5], out_vector[6])


