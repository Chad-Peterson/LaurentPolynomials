from laurentpolynomials.laurent_polynomial import LaurentPolynomial
from sympy import symbols, Poly, Expr

a = LaurentPolynomial('A')




# b = (a+1)**0
#
# print(b)

x = symbols('x')

# eq1 = Poly((x+1)*(x-1)*(x+1)**(-1))



#
# lp1 = as_laurent_poly(eq1)
#
# print(lp1)
#
# print(as_sympy_poly(lp1.coefficients, lp1.exponents))

# eq1 = (1 - a)
# eq2 = (1 + -a)
# eq3 = (1 - (-a))
#
#
# eq1._sort()
# eq2._sort()
# eq3._sort()
#
# print('eq1:', eq1)
# print('eq2:', eq2)
# print('eq3:', eq3)



# a / (2 * a)

x1, x2 = symbols('x1 x2')
syms = (x1, x2)    #  specify the desired order of symbols: do not rely on default order being what you want
expr = 3*x1**(-2) - x1*x2 + 5*x2


# todo fix as lp for this...
Poly(expr).as_dict()