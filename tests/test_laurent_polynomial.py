from laurentpolynomials.laurent_polynomial import LaurentPolynomial


def test__init__():
    a = LaurentPolynomial('A')

    # Creating an object with no arguments should create the indeterminate
    assert repr(a) == 'A'

    # Creating objects with the same exponents should simplify
    a1 = LaurentPolynomial('A', (1, 1, 1), (1, 1, 1))
    assert a1 == 3*a


def test__str__():
    pass


def test__repr__():
    pass


def test__eq__():
    a = LaurentPolynomial('A')

    assert a == a
    assert a == a + 0
    assert a == a - 0


def test__ne__():
    a = LaurentPolynomial('A')

    assert a != 0
    assert a != 1
    assert a != 1 + a


def test__neg__():
    a = LaurentPolynomial('A')

    assert -a == -a
    assert -a == a * (-1)
    assert -a == a - 2 * a


def test__add__():
    a = LaurentPolynomial('A')

    assert a + a == 2 * a

    # Associative rule
    assert (a + 1) + 2 == a + (1 + 2)
    assert (a + 1) + 2 == (a + 2) + 1
    assert (a + 1) + 2 == a + 3
    assert (a + 1) + 2 == 3 + a

    # Identity rule
    assert a + 0 == 0 + a == a

    # Additive inverse rule
    assert a + (-a) == 0
    assert a + (-a) == a - a
    assert a + (-a) == (-a) + a


def test__radd__():
    pass


def test__iadd__():
    pass


def test__sub__():
    a = LaurentPolynomial('A')

    assert a - a == 0
    assert a - 1 == a - 1
    assert a - 1 == a + (-1)
    assert a - 1 == a + (-1)


def test__rsub__():
    a = LaurentPolynomial('A')

    assert 1 - a == 1 - a
    assert 1 - a == 1 + (-a)
    assert 1 - (-a) == 1 + a
    assert 1 - 2*a**2 + 2*a == -2*a**2 + 2*a + 1
    assert 1 - (2*a**2 + 2*a) == -2*a**2 - 2*a + 1


def test__isub__():
    pass


def test__mul__():
    a = LaurentPolynomial('A')

    assert a * a == a ** 2
    assert a * 2 == 2 * a
    assert a * 2 == a + a

    # Identity rule
    assert a * 1 == 1 * a == a


def test__mul__compound():
    a = LaurentPolynomial('A')

    assert a * (2*a) == 2*a**2
    assert 2*a * (2*a) == 4*a**2
    assert (2*a + 1) * (2*a + 1) == 4*a**2 + 4*a + 1



def test__rmul__():

    # TODO Fix for zero

    a = LaurentPolynomial('A')

    # assert a * 0 == LaurentPolynomial('A',(1,),(0,))

    # assert a * 1 == 1

    assert 2 * a == 2 * a
    assert 2 * a == a * 2
    assert 2 * a == a + a


def test__rmul__compound():
    a = LaurentPolynomial('A')
    assert 2 * 2 * a == 4 * a
    assert 2 * (2*a) == 4 * a
    assert 2 * (2*a + 1) == 4 * a + 2
    assert 2 * (2*a + 1) == 4 * (a + 1/2)
    # assert 0 * a == 0


def test__imul__():
    pass


def test__truediv__():
    a = LaurentPolynomial('A')

    assert a / a == 1
    assert a / 2 == a * (1 / 2)
    # assert a/ (2*a) == 1 / 2
    # assert a / (2*a) == 1 / (2*a)

def test__truediv__compound():
    pass


def test__rtruediv__():
    a = LaurentPolynomial('A')

    assert 2 / a == 2 * (1 / a)
#
# def test__itruediv__():
#     pass

def test__pow__():
    # TODO Test raising exceptions
    a = LaurentPolynomial('A')

    assert a ** 2 == a * a
    assert a ** 3 == a * a * a
    assert a ** 0 == 1
    assert a ** -1 == 1 / a
    assert (a+1) ** 2 == (a+1) * (a+1)

def test_simplify_expression():
    a = LaurentPolynomial('A')

    assert a + a == 2 * a
    assert a - a == 0
    assert a * a == a ** 2
    assert a / a == 1
    assert a ** 2 == a * a
    assert a ** 3 == a * a * a
    assert a ** 0 == 1
    assert a ** -1 == 1 / a
    assert (a+1) * (a+1) == (a+1) ** 2
    assert (a+1) * (a+1) == a**2 + 2*a + 1
    assert (a+1) * a == a**2 + a
    assert (a+1) * 2 == 2*a + 2

def test_sort():
    pass


