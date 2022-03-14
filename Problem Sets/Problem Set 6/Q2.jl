using LinearAlgebra

# a) Calculate k̄, c̄
# Parameters
alpha = 0.3
delta = 0.1
beta = 0.95
rho = 0.9
s_bar = 1 # !
# log v_t+1
# disturbance_to_productivity_variance = 0.007^2
k_bar = ((1 / beta - 1 + delta) / alpha)^(1 / (alpha - 1))
c_bar = s_bar * (k_bar^alpha) - delta * k_bar


# b) calculate a1,a2,b1,b2,b3 coefficients

# Wrap up all the coefficients in a struct
mutable struct LinearizedApproximationCoeffictients
    a1::Float64
    a2::Float64
    b1::Float64
    b2::Float64
    b3::Float64
end

linearizedApproximationCoefficients = LinearizedApproximationCoeffictients(0, 0, 0, 0, 0)
# Define formulas for coefficients
# a1 = βα(α - 1)s̄k̄^α - 1
linearizedApproximationCoefficients.a1 = beta * alpha * ((alpha - 1) * s_bar * (k_bar^(alpha - 1)))
# a2 = βαs̄k̄^α-1
linearizedApproximationCoefficients.a2 = beta * alpha * (s_bar * (k_bar^(alpha - 1)))
# b1 = 1 - δ + αs̄k̄^α-1
linearizedApproximationCoefficients.b1 = 1 - delta + (alpha * (s_bar * (k_bar^(alpha - 1))))
# b2 = s̄k̄^α-1
linearizedApproximationCoefficients.b2 = s_bar * (k_bar^(alpha - 1))
# b3 = - c̄ / k̄
linearizedApproximationCoefficients.b3 = -(c_bar / k_bar)

# c) compute the inverse of the matrix with coefficients, e.g.
#  [-1 0 0; b3 b1 b2; 0 0 p]^-1

b3, b1, b2 = (linearizedApproximationCoefficients.b3,
linearizedApproximationCoefficients.b1,
linearizedApproximationCoefficients.b2)

b_coefficient_matrix = [
    -1 0 0
    b3 b1 b2
    0 0 rho
]

b_coefficient_matrix_inverse = inv(coefficient_matrix)

# d) get the matrix A
a1, a2 = (
    linearizedApproximationCoefficients.a1,
    linearizedApproximationCoefficients.a2
)
a_coefficient_matrix = [
    -1 a1 a2;
    0 1 0;
    0 0 1
]

A = b_coefficient_matrix_inverse * a_coefficient_matrix

# e) decompose the matrix A into
# A = QΛQ-¹
Lambda_Temp, Q_temp = eigen(A)
Lambda_diagonal = diagm(Lambda_Temp)

Q_inverse = inv(Q_temp)

# Use this to compute motions, coefficients and such. This is the toolkit of DSGE.
# Unfortunatly, I can't solve the problem to it's end (Problem set 6, Question 2)