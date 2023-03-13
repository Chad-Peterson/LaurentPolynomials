"""

"""

from typing import Union
from warnings import warn
import copy
import sympy
from sympy import symbols, Poly, Expr


class InputValidation:
    """
    TODO Can exponents be any real number? Only integers?
    """
    def __init__(self,
                 indeterminate:  str,
                 coefficients:   tuple[Union[int, float]],
                 exponents:      tuple[Union[int, float]],
                 order:          str):

        self.indeterminate = self._validate_indeterminate(indeterminate)
        self.coefficients  = self._validate_coefficients(coefficients)
        self.exponents     = self._validate_exponents(exponents)
        self.order         = self._validate_order(order)

    @staticmethod
    def _validate_indeterminate(indeterminate: str) -> str:

        if type(indeterminate) != str:
            raise TypeError("The indeterminate must be a string")

        if indeterminate == 'x':
            raise ValueError("Indeterminate 'x' is reserved. Please choose another indeterminate")

        if len(indeterminate) != 1:
            raise ValueError("The indeterminate must be a single character")

        return indeterminate

    @staticmethod
    def _validate_coefficients(coefficients: tuple[Union[int, float]]) -> tuple[Union[int, float]]:

        if type(coefficients) in [int, float]:
            coefficients = (coefficients,)

        if type(coefficients) != tuple:

            if type(coefficients) == list:
                # warn("The coefficients should be a tuple, not a list. Converting to a tuple")
                coefficients = tuple(coefficients)
            else:
                raise TypeError("The coefficients must be a tuple")

        for coefficient in coefficients:
            if type(coefficient) not in [int, float]:
                raise TypeError("The coefficients must be integers or floats")

        return coefficients

    @staticmethod
    def _validate_exponents(exponents: tuple[Union[int, float]]) -> tuple[Union[int, float]]:

        # TODO Can Laurent polynomials have negative exponents?
        # TODO Can Laurent polynomials have zero exponents?
        # TODO Can Laurent polynomials have non-integer exponents?

        if type(exponents) in [int, float]:
            # warn("The exponents should be a tuple, not a single number. Converting to a list")
            exponents = (exponents,)

        if type(exponents) != tuple:

            if type(exponents) == list:
                # warn("The exponents should be a tuple, not a list. Converting to a tuple")
                exponents = tuple(exponents)
            else:
                raise TypeError("The exponents must be a tuple")

        for exponent in exponents:
            if type(exponent) not in [int, float]:
                raise TypeError("The exponents must be integers or floats")

        return exponents

    @staticmethod
    def _validate_order(order: str) -> str:

        if type(order) != str:
            raise TypeError("The order must be a string")

        if order not in ['increasing', 'decreasing']:
            raise ValueError("The order must be either 'increasing' or 'decreasing'")

        return order

    def _format_polynomial(self, user_input):
        """_check_input_format Converts input to LaurentPolynomial for arithmetic operations

        :param user_input: _description_
        :type user_input: _type_
        :raises TypeError: _description_
        :return: _description_
        :rtype: _type_
        """

        if isinstance(user_input, LaurentPolynomial):
            return user_input

        elif isinstance(user_input, int) or isinstance(user_input, float):
            # warn("Converting input to a LaurentPolynomial with a single term")
            coefficients = (user_input,)
            exponents    = (0,)
            return LaurentPolynomial(self.indeterminate, coefficients=coefficients, exponents=exponents)

        else:
            raise TypeError("Polynomial must be of type LaurentPolynomial, int, or float")


class LaurentPolynomial(InputValidation):

    def __init__(self,
                 indeterminate: str,
                 coefficients:  tuple[Union[int, float]] = (1,),
                 exponents:     tuple[Union[int, float]] = (1,),
                 order:         str = 'increasing',
                 normalize:     bool = False):
        """
        LaurentPolynomial Class for representing Laurent polynomials

        TODO Make Laurent Polynomials subscriptable by monomial
        TODO Can we use real numbers in the expression?
        TODO Implement __floordiv__ method?
        TODO Implement reverse division methods
        TODO Implement _normalize method
        TODO finish switching lists to tuples in the various methods
        TODO __rpow__ method: Can we raise a number to an indeterminate?
        TODO Fix (1-a) -> (a-1) something is wrong with ordering.

        :param indeterminate: The character used to represent the polynomial
        :type indeterminate: str
        :param coefficients: The coefficients of the polynomial
        :type coefficients: tuple
        :param exponents: The orders of the polynomial
        :type exponents: tuple
        param order: The order in which the polynomial is represented (i.e., "increasing" or "decreasing")
        :type order: str
        :param normalize: Whether or not to normalize the polynomial
        :type normalize: bool
        """
        super(LaurentPolynomial, self).__init__(indeterminate, coefficients, exponents, order)

        self._sort()

        self._simplify_expression()

        self._normalize()


    def __repr__(self):

        # TODO Implement the self.order kwargs

        # If the polynomial only contains a single indeterminate with a coefficient of zero, then the polynomial is zero
        if self.coefficients == [0] and self.exponents == [0]:
            return "0"

        else:

            term_index = 0  # The index of the current term
            polynomial = ""

            for coefficient, exponent in zip(self.coefficients, self.exponents):

                # If the coefficient is zero, then don't show the term


                # First, determine whether a leading + or - sign is needed


                # If this is the first term, then don't show the addition sign (e.g., +1+x**2 --> 1+x**2)
                # Note: It is okay to display a leading minus sign
                if term_index == 0 and coefficient >= 0:
                    pass
                # elif coefficient > 0:
                #     polynomial += "+"
                elif coefficient == 0:
                    pass
                elif coefficient < 0:
                    polynomial += "-"
                else:
                    polynomial += "+"


                # Next, determine how to show the coefficient


                # If the coefficient is zero, then the whole term is zero --> don't show it
                if coefficient == 0:
                    pass

                # Don't show a coefficient of one if the exponent is nonzero (e.g., 1x**2 --> x**2)
                elif coefficient == 1 and exponent != 0:
                    pass

                elif coefficient == -1 and exponent != 0:
                    pass

                # Otherwise show the coefficient (strip minus signs since this is already accounted for)
                else:
                    polynomial += str(abs(coefficient))

                # Next, determine how to show the term

                # If the exponent is zero, then don't show the term (e.g., 3x**0 --> 3)
                if exponent == 0:
                    pass

                else:
                    polynomial += self.indeterminate


                # Finally, determine how to show the exponent


                # Don't show an exponent of zero or one (e.g., A**0 --> A, A**1 --> A)
                if exponent == 0 or exponent == 1:
                    pass
                else:
                    polynomial += "**" + str(exponent)

                term_index += 1

            return polynomial

    def __eq__(self, other):

        other = self._format_polynomial(other)

        if self.coefficients == other.coefficients and self.exponents == other.exponents:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __neg__(self):
        return LaurentPolynomial(self.indeterminate, coefficients=tuple([-1 * coefficient for coefficient in self.coefficients]), exponents=self.exponents)

    def __gt__(self, other):

        other = self._format_polynomial(other)

        if self._is_constant() and other._is_constant():
            return self.coefficients[0] > other.coefficients[0]

        return self.degree > other.degree

    def __add__(self, addend):

        # Ensure the input being added (the addend) is a LaurentPolynomial or throw an error
        addend = self._format_polynomial(addend)
        
        sum_coefficients = self.coefficients.copy()
        sum_exponents = self.exponents.copy()

        for addend_coefficient, addend_exponent in zip(addend.coefficients, addend.exponents):

            if addend_exponent in sum_exponents:

                # Get the index of the term in the polynomial
                term_index = sum_exponents.index(addend_exponent)

                # Increment the corresponding coefficient
                sum_coefficients[term_index] += addend_coefficient

            else:
                sum_coefficients.append(addend_coefficient)
                sum_exponents.append(addend_exponent)

        return LaurentPolynomial(self.indeterminate, coefficients=sum_coefficients, exponents=sum_exponents)

    def __radd__(self, addend):
        return self.__add__(addend)

    def __sub__(self, subtrahend):
        """
        Subtraction is just addition with a negative sign.
        E.g., x-1 == x + (-1)
        """
        return self.__add__(-1*subtrahend)

    def __rsub__(self, subtrahend):
        """
        When subtracting in reverse, the subtrahend is the negative of the minuend.
        E.g., 1 - x == (-1)*x + 1
        """
        self = (-1) * self
        return self.__add__(subtrahend)

    def __mul__(self, multiple):

        # Ensure the input being multiple (the multiple) is a LaurentPolynomial or throw an error
        multiple = self._format_polynomial(multiple)

        product_coefficients = []
        product_exponents = []

        for multiple_coefficient, multiple_exponent in zip(multiple.coefficients, multiple.exponents):


            for coefficient, exponent in zip(self.coefficients, self.exponents):

                product_coefficients.append(multiple_coefficient * coefficient)
                product_exponents.append(multiple_exponent + exponent)

        return LaurentPolynomial(self.indeterminate, coefficients=product_coefficients, exponents=product_exponents)

    def __rmul__(self, multiple):
        return self.__mul__(multiple)


    def __pow__(self, power):
        """


        :param power:
        :return:
        """

        power = self._format_polynomial(power)

        if not power._is_monomial():
            raise TypeError("Power must be a monomial")

        if not power._is_integer():
            raise TypeError("Power must be an integer")

        if power._is_positive():
            sp = self._as_sympy_poly()
            power = power._as_constant()
            sp = Poly(sp.pow(power))
            return self._as_laurent_poly(sp)

        elif power._is_zero():
            return LaurentPolynomial(self.indeterminate, coefficients=[1], exponents=[0])

        elif power._is_negative():
            sp = self._as_sympy_poly()
            power = power._as_constant()
            sp = Poly(sp.pow(-1 * power))
            return self._as_laurent_poly(sp, sign='negative')
        else:
            raise TypeError("IDK")


    def __rpow__(self, base):
        raise NotImplementedError


    def __truediv__(self, divisor):

        divisor_copy = divisor
        # TODO Come up with more reliable logic

        if isinstance(divisor, LaurentPolynomial) and divisor._is_monomial():
            return self.__mul__(divisor.__pow__(-1))

        if isinstance(divisor, LaurentPolynomial):

            divisor = divisor._as_sympy_poly()
            sp = self._as_sympy_poly()

            # If the resulting expression is a constant, then return a constant (int or float)
            if isinstance(sp/divisor, sympy.core.numbers.Integer):
                return int(sp/divisor)

            elif isinstance(sp/divisor, sympy.core.numbers.Float):
                return float(sp/divisor)

            elif isinstance(sp/divisor, sympy.core.mul.Mul):
                # TODO Handle negative terms

                num, denom, = (sp/divisor).as_numer_denom()

                # If numerator is a constant, then the denominator is
                if num.is_number:
                    self = self.__pow__(-1)
                    return self.__mul__(1/divisor_copy)

                sp = Poly(sp/divisor)
                return self._as_laurent_poly(sp)

            elif isinstance(sp/divisor, sympy.polys.polytools.Poly):
                return self._as_laurent_poly(sp/divisor)

            else:
                sp = float(sp/divisor)
                return self._as_laurent_poly(sp)

        else:

            # TODO Fix __truediv__ method

            # TODO is this necessary?
            # Ensure the input being divided (the divisor) is a LaurentPolynomial or throw an error
            # divisor = self._format_polynomial(divisor)

            if repr(divisor) == "0":
                raise ZeroDivisionError("Cannot divide by zero")

            sp = self._as_sympy_poly()
            sp = Poly(sp / divisor)

            # TODO just do 1/divisor and then multiply by self?
            return self._as_laurent_poly(sp)

    def __rtruediv__(self, divisor):
        """
        Reverse multiplication is equivalent to multiplying the constant by an inverted indeterminate.
        """
        self = self.__pow__(-1)
        return self.__mul__(divisor)


    def _is_monomial(self):
        """
        :return: True if the polynomial is a monomial, False otherwise
        """
        if len(self.exponents) == 1:
            return True
        else:
            return False

    def _is_constant(self):
        """
        :return: True if the polynomial is a constant, False otherwise
        """

        if not self._is_monomial():
            raise ValueError("Polynomial is not a monomial")

        if self.exponents == [0]:
            return True
        else:
            return False

    def _is_zero(self):
        """
        :return: True if the polynomial is zero, False otherwise
        """
        if self.coefficients == [0]:
            return True
        else:
            return False

    def _is_positive(self):
        """
        :return: True if the polynomial is positive, False otherwise
        """

        if len(self.coefficients) > 1:
            raise ValueError("Polynomial is not a constant")

        if self.coefficients[0] > 0:
            return True
        else:
            return False

    def _is_negative(self):
        """
        :return: True if the polynomial is negative, False otherwise
        """

        if len(self.coefficients) > 1:
            raise ValueError("Polynomial is not a constant")

        if self.coefficients[0] < 0:
            return True
        else:
            return False

    def _is_integer(self):
        """
        :return: True if the polynomial is an integer, False otherwise
        """

        if len(self.coefficients) > 1:
            raise ValueError("Polynomial is not a constant and therefore cannot be an integer")

        if isinstance(self.coefficients[0], int):
            return True

        else:
            return False

    def _is_rational(self):
        """
        :return: True if the polynomial is rational, False otherwise
        """
        if self._is_constant():
            if self._is_integer():
                return True
            else:
                return False
        else:
            return False

    def _as_constant(self):
        """
        :return: The polynomial as an integer
        """
        if self._is_constant():
            return int(self.coefficients[0])
        else:
            raise ValueError("Polynomial is not a constant")

    @staticmethod
    def _as_laurent_poly(poly, sign='positive'):

        if isinstance(poly, sympy.polys.polytools.Poly):

            poly_dict = poly.as_dict()

            if sign == 'positive':
                coefficients = tuple(float(coeff) for coeff in list(poly_dict.values()))
                exponents = tuple(int(exp[0]) for exp in list(poly_dict.keys()))
            elif sign == 'negative':
                coefficients = tuple(1/float(coeff) for coeff in list(poly_dict.values()))
                exponents = tuple(int(-exp[0]) for exp in list(poly_dict.keys()))
            else:
                raise ValueError("Sign must be either 'positive' or 'negative'")


            return LaurentPolynomial('A', coefficients, exponents)

        else:
            # TODO Other than int?
            return float(poly)



    def _as_sympy_poly(self):
        x = symbols('x')
        lp_dict = {(exp,): coeff for exp, coeff in zip(self.exponents, self.coefficients)}
        return Poly(lp_dict, x)

    def _sort(self):
        """
        Sort the coefficients and exponents by exponent
        :return:
        """
        self.exponents, self.coefficients = zip(*sorted(zip(self.exponents, self.coefficients)))

    def _simplify_expression(self):
        
        new_coefficients = []
        new_exponents = []

        for coefficient, exponent in zip(self.coefficients, self.exponents):

            if coefficient == 0 and exponent == 0:
                pass

            elif coefficient != 0 and exponent not in new_exponents:
                new_coefficients.append(coefficient)
                new_exponents.append(exponent)

            elif coefficient != 0 and exponent in new_exponents:
                new_coefficients[new_exponents.index(exponent)] += coefficient

            else:
                pass

        # If the lists are completely empty, add a zero term
        if len(new_coefficients) == 0 and len(new_exponents) == 0:
            new_coefficients.append(0)
            new_exponents.append(0)

        self.coefficients, self.exponents = new_coefficients, new_exponents

    def _normalize(self):
        """
        TODO Implement the _normalize method
        """
        pass


# A = LaurentPolynomial('A')
#
# aa = (A+1)*(A+1)
#
# print(aa)

# print('A', A)
#
# a = A+2
#
# print('a', a)
#
# a - A
#
# print('a', a)
#
# b = 2+A
#
#
# print('b', b)
#
# A**2

# b = 3*A

# aa = LaurentPolynomial([(1,1), (1,2), (1,3)])

# bb = LaurentPolynomial([(1,1), (1,2), (1,3)])

# print('polynomial', aa)

# print(type(aa))

# cc = aa + bb

# print('add self  ',cc)
# print(type(cc))

# print(cc.coeff)

# dd = aa * bb

# print('mult self ',dd)


print('Done')