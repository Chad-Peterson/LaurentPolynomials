"""

"""

from typing import Union
from warnings import warn
import copy
from sympy import symbols, Poly


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

        """
        (1+x) * (1+x) = (1+x)**2 = 1 + 2x + x**2

        TODO check if need to expand...
        TODO fix object creation statement
        TODO Fix 0*self should be 0 (and simplify like terms before returning new object)
        
        """
        
        # Ensure the input being multiple (the multiple) is a LaurentPolynomial or throw an error
        multiple = self._format_polynomial(multiple)

        # If either polynomial is zero, then the product is zero
        if repr(multiple) == "0":
            return LaurentPolynomial(self.term, coefficients=[0], exponents=[0])

        product_coefficients = []
        product_exponents = []

        for multiple_coefficient, multiple_exponent in zip(multiple.coefficients, multiple.exponents):


            for coefficient, exponent in zip(self.coefficients, self.exponents):

                product_coefficients.append(multiple_coefficient * coefficient)
                product_exponents.append(multiple_exponent + exponent)

        return LaurentPolynomial(self.indeterminate, coefficients=product_coefficients, exponents=product_exponents)


    def __imul__(self, multiple):
        return self.__mul__(multiple)

    def __rmul__(self, multiple):
        return self.__mul__(multiple)


    def __pow__(self, power):
        """

        TODO Verify negative exponentiation works properly
        TODO Is deep copy necessary?
        TODO Remove special case for power = 0
        (1+x)**2 = 1 + 2x + x**2
        :param power:
        :return:
        """

        sp = self._as_sympy_poly()

        if power >= 0:
            sp = Poly(sp.pow(power))
            return self._as_laurent_poly(sp)
        else:
            sp = Poly(sp.pow(-1 * power))
            # sp = 1 / sp
            return self._as_laurent_poly(sp, sign='negative')


    def __rpow__(self, base):
        raise NotImplementedError


    def __truediv__(self, divisor):
        # TODO Fix __truediv__ method

        # TODO is this necessary?
        # Ensure the input being divided (the divisor) is a LaurentPolynomial or throw an error
        # divisor = self._format_polynomial(divisor)

        if repr(divisor) == "0":
            raise ZeroDivisionError("Cannot divide by zero")

        sp = self._as_sympy_poly()
        sp = sp / divisor

        # TODO just do 1/divisor and then multiply by self?
        return self._as_laurent_poly(sp)

    def __rtruediv__(self, divisor):


        return self.__pow__(-divisor)


    @staticmethod
    def _as_laurent_poly(poly, sign='positive'):

        poly_dict = poly.as_dict()

        if sign == 'positive':
            coefficients = tuple(int(coeff) for coeff in list(poly_dict.values()))
            exponents = tuple(int(exp[0]) for exp in list(poly_dict.keys()))
        elif sign == 'negative':
            coefficients = tuple(int(coeff) for coeff in list(poly_dict.values()))
            exponents = tuple(int(-exp[0]) for exp in list(poly_dict.keys()))
        else:
            raise ValueError("Sign must be either 'positive' or 'negative'")



        return LaurentPolynomial('A', coefficients, exponents)

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