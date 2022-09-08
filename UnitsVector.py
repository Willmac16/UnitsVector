### Will MacCormack 2022
### Inspired by https://youtu.be/bI-FS7aZJpY and a desire not to write bad code


# TODO: Switch from dumb units conversion to "extract basis"
# Store alternatives in terms of a units vector and a scale vector of MKS
import numpy as np
import math

mks_units = ["s", "m", "kg", "A", "K", "mol", "cd"]


# Widest Class: should be broader than MKS/SI units
class UnitsVector:
    def __init__(self, vec, val, val_scale):
        self.value = val
        self.vector = vec
        self.value_scale = val_scale

    def __units__(self):
        out = ""
        for i in range(len(self.vector)):
            if self.vector[i] != 0:
                out += mks_units[i] + "^" + str(self.vector[i]) + " "

        return out

    def __neg__(self):
        return UnitsVector(self.vector, -self.value, self.value_scale)

    def __add__(self, other):
        if isinstance(other, UnitsVector):
            if np.equal(self.vector, other.vector).all():
                abs_value = self.value * self.value_scale + other.value * other.value_scale
                relative_val = abs_value / self.value_scale

                return UnitsVector(self.vector, relative_val, self.value_scale)
            else:
                raise Exception("UnitsVector addition error: units do not match ({} | {})".format(self.__units__(), other.__units__()))
        elif isinstance(other, (int, float)):
            if np.linalg.norm(self.vector) == 0:
                out_val = self.value * self.value_scale + other
                return UnitsVector(self.vector, out_val / self.value_scale, self.value_scale)
            else:
                raise Exception("UnitsVector addition error: cannot add something with Units to a number ({})".format(self.__units__()))
        else:
            raise Exception("UnitsVector addition error: other is not a UnitsVector")

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, UnitsVector):
            if np.equal(self.vector, other.vector).all():
                abs_value = self.value * self.value_scale - other.value * other.value_scale
                relative_val = abs_value / self.value_scale

                return UnitsVector(self.vector, relative_val, self.value_scale)

            else:
                raise Exception("UnitsVector subtraction error: units do not match ({} | {})".format(self.__units__(), other.__units__()))
        elif isinstance(other, (int, float)):
            if np.linalg.norm(self.vector) == 0:
                out_val = self.value * self.value_scale - other
                return UnitsVector(self.vector, out_val / self.value_scale, self.value_scale)
            else:
                raise Exception("UnitsVector subtraction error: cannot subtract a number from something with Units")
        else:
            raise Exception("UnitsVector subtraction error: other is not a UnitsVector")

    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            if np.linalg.norm(self.vector) == 0:
                out_val = other - self.value * self.value_scale
                return UnitsVector(self.vector, out_val / self.value_scale, self.value_scale)
            else:
                raise Exception("UnitsVector subtraction error: cannot subtract something with Units from a number")
        else:
            raise Exception("UnitsVector subtraction error: other is not a UnitsVector")

    def __mul__(self, other):
        if isinstance(other, UnitsVector):
            out_vec = self.vector + other.vector
            out_val = self.value * other.value
            out_scale = self.value_scale * other.value_scale

            return UnitsVector(out_vec, out_val, out_scale)
        elif isinstance(other, (int, float)):
            out_val = self.value * other
            return UnitsVector(self.vector, out_val, self.value_scale)
        else:
            raise Exception("UnitsVector multiplication error: other is not a UnitsVector or a number")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, UnitsVector):
            out_vec = self.vector - other.vector
            out_val = self.value / other.value
            out_scale = self.value_scale / other.value_scale

            return UnitsVector(out_vec, out_val, out_scale)
        elif isinstance(other, (int, float)):
            out_val = self.value / other
            return UnitsVector(self.vector, out_val, self.value_scale)
        else:
            raise Exception("UnitsVector division error: other is not a UnitsVector or a number")

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            out_val = other / self.value
            return UnitsVector(self.vector * -1, out_val, 1 / self.value_scale)
        else:
            raise Exception("UnitsVector division error: other is not a UnitsVector or a number")

    # I don't feel like dealing with the ternary form with modulo
    def __pow__(self, power):
        if isinstance(power, (int, float, np.int32, np.int64, np.float32, np.float64, np.float16)):
            out_val = self.value ** power
            return UnitsVector(self.vector * power, out_val, self.value_scale ** power)
        elif isinstance(power, UnitsVector):
            if np.linalg.norm(power.vector) == 0:
                out_val = self.value ** (power.value * power.value_scale)
                return UnitsVector(self.vector * power.value, out_val, self.value_scale ** power.value)
            else:
                raise Exception("UnitsVector power error: power must be unitless")
        else:
            print(type(power))
            raise Exception("UnitsVector power error: power is not a number")

    def __rpow__(self, base):
        if isinstance(base, (int, float)):
            if np.linalg.nomr(self.vector) == 0:
                out_val = base ** (self.value * self.value_scale)
                return UnitsVector(np.zeros(7), out_val, 1)
            else:
                raise Exception("UnitsVector power error: power must be unitless")
        else:
            raise Exception("UnitsVector power error: base is not a number")

    def __value__(self):
        return self.value * self.value_scale

    def __lt__(self, other):
        if isinstance(other, UnitsVector):
            other = other.__value__()
        elif not isinstance(other, (int, float)):
            raise Exception("UnitsVector comparison error: other is not a UnitsVector or a number")

        return self.__value__() < other

    def __le__(self, other):
        if isinstance(other, UnitsVector):
            other = other.__value__()
        elif not isinstance(other, (int, float)):
            raise Exception("UnitsVector comparison error: other is not a UnitsVector or a number")

        return self.__value__() <= other

    def __gt__(self, other):
        if isinstance(other, UnitsVector):
            other = other.__value__()
        elif not isinstance(other, (int, float)):
            raise Exception("UnitsVector comparison error: other is not a UnitsVector or a number")

        return self.__value__() > other

    def __ge__(self, other):
        if isinstance(other, UnitsVector):
            other = other.__value__()
        elif not isinstance(other, (int, float)):
            raise Exception("UnitsVector comparison error: other is not a UnitsVector or a number")

        return self.__value__() >= other

    def __abs__(self):
        if self < 0:
            return self * -1
        else:
            return self

    def sqrt(self):
        return UnitsVector(self.vector / 2, math.sqrt(self.value), math.sqrt(self.value_scale))

    # I plan on printing in mks unless explicitly casted to another unit system
    def __repr__(self):
        out_str = str(self.value * self.value_scale)
        out_str += " ("
        spaces = False
        for i in range(len(self.vector)):
            if self.vector[i] != 0:
                if spaces:
                    out_str += " "
                out_str += mks_units[i] + "^" + str(self.vector[i])
                spaces = True
        out_str += ")"
        return out_str

    def __str__(self):
        out_str = "{:.5f}".format(self.value * self.value_scale)
        out_str += " ("
        spaces = False
        for i in range(len(self.vector)):
            if self.vector[i] != 0:
                if spaces:
                    out_str += " "
                out_str += mks_units[i] + "^" + str(self.vector[i])
                spaces = True
        out_str += ")"
        return out_str

    def __float__(self):
        if np.linalg.norm(self.vector) == 0:
            return self.value * self.value_scale
        else:
            raise Exception("UnitsVector float error: cannot convert to float if not unitless")

    def __format__(self, spec):
        return f'{self.__value__():{spec}}'  + " ( " + self.__units__() + ")"


class MKS(UnitsVector):
    def __init__(self, value, s, m, kg, A, k, mol, cd):
        self.value = value
        self.value_scale = 1

        self.vector = np.array([s, m, kg, A, k, mol, cd])

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
imperial_scale_vector = np.array((seconds_to_seconds, feet_to_meters, pound_mass_to_kilograms, amps_to_amps, rankine_to_kelvins, pound_moles_to_moles, candela_to_candela))

imp_units = ["s", "ft", "lbm", "A", "R", "mol", "cd"]

class Imperial(UnitsVector):
    def __init__(self, value, s, ft, lbm, A, R, lb_mol, cd):
        self.value = value

        base_vector = np.array([s, ft, lbm, A, R, lb_mol, cd])


        # self.seconds = s
        # self.feet = ft
        # self.pound-mass = lbm
        # self.amps = A
        # self.rankine = R
        # self.moles = lb_mol
        # self.candela = cd

        self.vector = base_vector
        self.value_scale = imperial_scale_vector @ base_vector

    def __str__(self):
        out_str = str(self.value * self.value_scale)
        out_str += " ("
        spaces = False
        for i in range(len(self.vector)):
            if self.vector[i] != 0:
                if spaces:
                    out_str += " "
                out_str += imp_units[i] + "^" + str(self.vector[i])
                spaces = True
        out_str += ")"
        return out_str

# Imperial to MKS and Back
def toImperial(value: UnitsVector):
    out_value = value * value.value_scale / (imperial_scale_vector @ value.vector)
    out_vector = value.vector

    return Imperial(out_value, out_vector[0], out_vector[1], out_vector[2], out_vector[3], out_vector[4], out_vector[5], out_vector[6])

def toMKS(value: UnitsVector):
    out_value = value * value.value_scale
    out_vector = value.vector

    return MKS(out_value, out_vector[0], out_vector[1], out_vector[2], out_vector[3], out_vector[4], out_vector[5], out_vector[6])

# Specific Units Subclasses for Ease of Use
# Gonna be MKS extensions
# If I want an imperial printout I'll cast later
class Seconds(MKS):
    def __init__(self, t):
        super().__init__(t, 1, 0, 0, 0, 0, 0, 0)

class Meters(MKS):
    def __init__(self, x):
        super().__init__(x, 0, 1, 0, 0, 0, 0, 0)

class Centimeters(Meters):
    def __init__(self, x):
        super().__init__(x/100.0)

class Millimeters(Meters):
    def __init__(self, x):
        super().__init__(x/1000.0)

class Kilometers(Meters):
    def __init__(self, x):
        super().__init__(x*1000.0)

class Inches(Meters):
    def __init__(self, x):
        super().__init__(x * 0.0254)

class Feet(Inches):
    def __init_(self, x):
        super().__init__(x * 12.0)

class Kilograms(MKS):
    def __init__(self, m):
        super().__init__(m, 0, 0, 1, 0, 0, 0, 0)

class Grams(Kilograms):
    def __init__(self, m):
        super().__init__(m/1000.0)

class Pounds(Kilograms):
    def __init__(self, m):
        super().__init__(m * pound_mass_to_kilograms)

class Omega(MKS):
    def __init__(self, w):
        super().__init__(w, -1, 0, 0, 0, 0, 0, 0)

class RPM(Omega):
    def __init__(self, w):
        super().__init__(w * 2 * math.pi / 60)

class Newtons(MKS):
    def __init__(self, f):
        super().__init__(f, -2, 1, 1, 0, 0, 0, 0)

class PoundsForce(Newtons):
    def __init__(self, f):
        super().__init__(f / 4.4482216152605)

class MetersPerSecond(MKS):
    def __init__(self, v):
        super().__init__(v, -1, 1, 0, 0, 0, 0, 0)

def kilometersPerHour(v):
    if isinstance(v, UnitsVector):
        vel_vec = np.array([-1, 1, 0, 0, 0, 0, 0])
        if (v.vector == vel_vec).all():
            return v.value * v.value_scale * 3.6
        else:
            raise Exception("UnitsVector error: cannot convert disimilar units")

class KilometersPerHour(MetersPerSecond):
    def __init__(self, v):
        super().__init__(v / 3.6)

class MilesPerHour(KilometersPerHour):
    def __init__(self, v):
        super().__init__(v * 1.60934)

class Radians(MKS):
    def __init__(self, t):
        super().__init__(t, 0, 0, 0, 0, 0, 0, 0)

class Degrees(Radians):
    def __init__(self, t):
        super().__init__(t * math.pi / 180)

class MetersPerSecondSquared(MKS):
    def __init__(self, a):
        super().__init__(a, -2, 1, 0, 0, 0, 0, 0)

class Watts(MKS):
    def __init__(self, p):
        super().__init__(p, -3, 1, 1, 0, 0, 0, 0)


def test():
    print("Units testing")

    ## Solve for how far a mass goes in 10s under constant force

    mass_pounds = 240.0
    force_pounds = 100.0
    time = 10.0

    # Basic Float Sanity Check
    mass_kilograms = mass_pounds * pound_mass_to_kilograms
    force_newtons = force_pounds / 4.4482216152605

    print("F:", 0.5 * force_newtons / mass_kilograms * time ** 2)

    # Pure Metric Sanity Check
    m = Kilograms(mass_kilograms)
    f = Newtons(force_newtons)
    t = Seconds(time)

    print("MKS:", 0.5 * f / m * t ** 2)

    # Imperial Sanity Check
    mI = Pounds(mass_pounds)
    fI = PoundsForce(force_pounds)
    tI = Seconds(time)

    print("I:", 0.5 * fI / mI * tI ** 2)

    # Unitlessness Check
    # As Prof Rashidi Always Says trig and exponential functions take unitless arguments
    # invalid
    try:
        print(math.sin(Seconds(2*math.pi)))
    except:
        print("Arg wasn't unitless")

    print(math.sin(2 * math.pi * Seconds(20) / Seconds(60)))









if __name__ == "__main__":
    test()
