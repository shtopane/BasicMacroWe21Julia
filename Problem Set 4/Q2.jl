using LinearAlgebra
using Plots; plotly()

include("./MacroHelperFunctions.jl")
using .MacroHelperFunctions

# Parameters
alpha = 0.3
beta = 0.95
sigma = 2
delta_old = 0.1
delta_new = 1
theta = (1/beta) - 1
Ā = 1
# Utility fun is CRRA
# U(ct) = c^1-sigma - 1 / 1 - sigma
# U`(c) = 1/c^σ

## Compute INITIAL Steady state values
# From MATLAB ( (1 - beta + delta_old*beta) / (alpha * beta * Ā))^(1/(alpha - 1))
# function getSteadyStateConsumption(steady_state_capital, alpha, delta)
#     (steady_state_capital^alpha) - (delta * steady_state_capital)
# end

# function getSteadyStateCapital(alpha, delta, theta)
#     (alpha / (delta + theta))^(1 / (1 - alpha)) # F'(k) = alpha*k^(alpha-1) = 1/beta + delta -1 
# end

k_star_old = getSteadyStateCapital(alpha, delta_old, theta)# (alpha / (delta_old + theta))^(1 / (1 - alpha)) # F'(k) = alpha*k^(alpha-1) = 1/beta + delta -1 
c_star_old = getSteadyStateConsumption(k_star_old, alpha, delta_old) # (k_star_old^alpha) - (delta_old * k_star_old)

## Compute NEW Steady state values after change in parameter

k_star_new = getSteadyStateCapital(alpha, delta_new, theta)
c_star_new = getSteadyStateConsumption(k_star_new, alpha, delta_new)

# Calculate deviation to use in the policy function
capital_deviation = k_star_old - k_star_new

# Getting entries in the A matrix
# U`(c*)
consumption_first_der = 1/(c_star_new^sigma)
# U``(c*)
consumption_second_der = -sigma*(c_star_new^(-sigma - 1))
# F``(k*)
capital_second_der = alpha * (alpha - 1) * k_star_new^(alpha - 2)

matrix_A = zeros(2,2)
matrix_A[1,1] = 1 + beta * (consumption_first_der * capital_second_der) / consumption_second_der
matrix_A[2,1] = -1
matrix_A[1,2] = -((consumption_first_der*capital_second_der)/ consumption_second_der)
matrix_A[2,2] = 1 + theta

eigenValues, eigenVectors = eigen(matrix_A)
stableEigenIndex = 0
stableEigenValue = 0

for i =1:length(eigenValues)
    if(eigenValues[i] < 1)
        global stableEigenIndex = i
        global stableEigenValue = eigenValues[i]
    end
end

T = 50
eigenVectorDivision = (eigenVectors[1, stableEigenIndex] / eigenVectors[2, stableEigenIndex])
simulation_wrapper = fill(0.0, (2,T))
simulation_wrapper[2,1] = capital_deviation
simulation_wrapper[1,1] = eigenVectorDivision * simulation_wrapper[2,1]

for i = 2:T
    simulation_wrapper[2,i] = dot(matrix_A[2, :],simulation_wrapper[:, i-1])
    simulation_wrapper[1, i] = eigenVectorDivision * simulation_wrapper[2,i]
end

con_cap_dev = plot(1:T, simulation_wrapper[1, :], label="Consumption deviation")
plot!(simulation_wrapper[2, :], label="Capital deviation")