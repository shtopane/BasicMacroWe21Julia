# This implements the solution of a
# multivariate rational expectations model,
# using the example of a basic RBC model,
# i.e., the model economy used in the lectures
# to explain the general solution procedure.
# The computations performed in the following code therefore
# correspond to pages 41 - 51 of the book
# by Roger E. A. Farmer,
# "Macroeconomics of Self-fulfilling Prophecies", 2nd edition, MIT Press, 1999.
# The comments refer to page numbers and equation numbers from the book.

using LinearAlgebra
using Plots;
using Random
plotly();
include("./helpers/DSGEHelperFunctions.jl")
using DSGEHelperFunctions

# Parameters used in the model
beta = 0.95 # discount factor
alpha = 0.3 # capital share
delta = 0.1 # depreciation rate
rho = 0.9 # autocorrelation of productivity

# Getting non-stochastic steady-state values
# See function implementation for details. Function form is used from Problem Set 6, Question 2
k_bar = DSGEHelperFunctions.getSteadyStateCapital(alpha, beta, delta)
s_bar = DSGEHelperFunctions.getSteadyStateProductivity(1)
c_bar = DSGEHelperFunctions.getSteadyStateConsumption(k_bar, alpha, delta, s_bar)

# Describing relationships underlying linear coefficients
# For functional form see DSGEHelperFunctions.jl
a1, a2, b1, b2, b3 = DSGEHelperFunctions.getLinearCoefficients(
    alpha, beta, delta, k_bar, c_bar
)

# Writing dynamic equilibrium conditions into matrix form
# the matrix m371 describes the matrix multiplying period-t variables Part 9, Slide 15
# -1 0 0
# b3 b1 b2
# 0 0 rho
m371 = zeros(3, 3)
m371[1, 1] = -1
m371[2, :] = [b3, b1, b2]
m371[3, 3] = rho

# m372 stands for the second matrix, multiplying t+1 variables: Part 9, Slide 15
m372 = zeros(3, 3)
m372[1, :] = [-1 a1 a2]
m372[2, 2] = 1
m372[3, 3] = 1

# m38 is the inverse of m371
m38 = inv(m371)
# Matrix A: Part 9, Slide 16
# Note that this encodes all the relevant dynamics for the solution
# This turns the solution in the form of (L1)
# ĉt      ĉt+1      v̂t+1
# k̂t  = A k̂t+1 +  B wᶜt+1
# ŝt      ŝt+1      wᵏt+1
#                   wˢt+1
# where w^c, w^k, w^s are expectational errors for the variables c,k,s
A = m38 * m372

# Calculate eigenvectors and eigenvalues of A
# the idea is to get an equation of the form A = QΛQ^-1
Lambda_temp, Q_temp = eigen(A)
Lambda_temp_sorted = sort(Lambda_temp)
Lambda = diagm(Lambda_temp_sorted) # Making diagonal matrix of the eigenValues
# TODO: Test
# Idea - the eigenvector associated with the stable eigenvalue(<1) should be at first position
Q_sorted = sort(Q_temp, dims = ndims(Q_temp))

# The inverse of Q as used in the decomposition a = QΛQ^-1,
# Part 10, Slide 4
Q_inv = inv(Q_sorted)

# The key restriction: getting a free variable(consumption)
# as an equilibrium function of predetermined variables(capital, productivity-state)
# Part 10, Slide 9
# coefficient of rational expectation ( capital ? )
# TODO: TEST
Q_new = [-0.4456 0.4649 0.5289
    0.8952 0.8854 0.8303
    0 0 0.1759]
Q_new_inv = inv(Q_new)
# TODO: END TEST
c_coeff_k = ((-Q_new_inv[1, 2]) / Q_new_inv[1, 1])
# coefficient of rational expectation ( consumption ? )
c_coeff_s = ((-Q_new_inv[1, 3]) / Q_new_inv[1, 1])

# log of disturbances to productivity variance
v_variance = 0.007^2
T = 100
Random.seed!(123)

shock_realizations_wrapper = zeros(T, 1)
shock_realizations_wrapper[2:end] = sqrt(v_variance) * randn(length(shock_realizations_wrapper) - 1, 1)

# sequence of autocorrelated productivity states
# This is doing a matrix not a vector!
productivity_states_wrapper = NaN * zeros(T, 1)
productivity_states_wrapper[1] = 0

for t = 2:T
    productivity_states_wrapper[t] = rho * productivity_states_wrapper[t-1] + shock_realizations_wrapper[t]
end

plot(productivity_states_wrapper)

# Generating time series of variables, which are equilibria for the given shock reaizations
c_t = NaN * zeros(T, 1)
k_t = NaN * zeros(T + 1, 1)
k_t[1] = 0

for j = 1:T
    # Imposing the key equilibrium restriction from rational expectations:
    # The free variable (consumption) is expressed as a function of
    # the predetermined variables (capital and productivity state)
    # Part 10, Slide 9.
    # This is the recursive structure described in Part 10, Slide 12.

    # This is the code from MATLAB
    # c_t[j] = (-Q_new_inv[1, 2] / Q_new_inv[1, 1]) * k_t[j] + (-Q_new_inv[1, 3] / Q_new_inv[1, 1]) * productivity_states_wrapper[j]
    # Alternative:
    c_t[j] = c_coeff_k * k_t[j] + c_coeff_s * productivity_states_wrapper[j]
    k_t[j+1] = b1 * k_t[j] + b2 * productivity_states_wrapper[j] + b3 * c_t[j]
end

# Plotting!
plot(productivity_states_wrapper, color=:red, label="Productivity state")
plot!(c_t,  color=:blue, label="Consumption")
plot!(k_t, color=:green, label="Capital")
title!("DGSE - Productivity, consumption and capital")
