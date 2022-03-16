# linearized solution of a Ramsey model with an unaticipated and 
# permanent shock to depreciation
# CRRA utility function
# u(c) = (c^(1 - \sigma) - 1)/(1 - \sigma),
# Cobb-Douglas production function
using LinearAlgebra
using Plots;
plotly();

include("../Problem Sets/Problem Set 4/MacroHelperFunctions.jl")
using .MacroHelperFunctions

# Parameters
alpha = 0.3
beta = 0.95
sigma = 2.0
delta_old = 0.1
delta_1 = 0.2
delta_2 = 0.3 # From 5 period onwards
theta = (1 / beta) - 1
Ā = 1

T = 50

# INITIAL Steady state
k_star_old = getSteadyStateCapital(alpha, delta_old, theta)
c_star_old = getSteadyStateConsumption(k_star_old, alpha, delta_old)

# NEW Steady state
k_star_1 = getSteadyStateCapital(alpha, delta_1, theta)
c_star_1 = getSteadyStateConsumption(k_star_1, alpha, delta_1)

# Steady state 2 - for delta 0.3
k_star_2 = getSteadyStateCapital(alpha, delta_2, theta)
c_star_2 = getSteadyStateConsumption(k_star_2, alpha, delta_2)


# Compite deviation of capital - needed for policy function
capital_deviation_1 = k_star_old - k_star_1
capital_deviation_2 = k_star_1 - k_star_2

# Get matrix entries for NEW regime
# U`(c*)
consumption_first_der = 1 / (c_star_1^sigma)
# U``(c*)
consumption_second_der = -sigma * (c_star_1^(-sigma - 1))
# F``(k*)
capital_second_der = alpha * (alpha - 1) * k_star_1^(alpha - 2)
matrix_A_1 = zeros(2, 2)
matrix_A_1[1, 1] = 1 + beta * (consumption_first_der * capital_second_der) / consumption_second_der
matrix_A_1[2, 1] = -1
matrix_A_1[1, 2] = -((consumption_first_der * capital_second_der) / consumption_second_der)
matrix_A_1[2, 2] = 1 + theta

eigenValues, eigenVectors = eigen(matrix_A_1)
stableEigenIndex = 0
stableEigenValue = 0

for i = 1:length(eigenValues)
    if (eigenValues[i] < 1)
        global stableEigenIndex = i
        global stableEigenValue = eigenValues[i]
    end
end

eigenVectorDivision = (eigenVectors[1, stableEigenIndex] / eigenVectors[2, stableEigenIndex])
simulation_wrapper = fill(0.0, (2, T))
simulation_wrapper[2, 1] = capital_deviation_1
simulation_wrapper[1, 1] = eigenVectorDivision * simulation_wrapper[2, 1]

for t = 2:5
    simulation_wrapper[2, t] = dot(matrix_A_1[2, :], simulation_wrapper[:, t-1])
    simulation_wrapper[1, t] = eigenVectorDivision * simulation_wrapper[2, t]
end

simulation_wrapper[2,6] = capital_deviation_2
simulation_wrapper[1, 6] = eigenVectorDivision * simulation_wrapper[2, 6]

for t = 7:T
    simulation_wrapper[2, t] = dot(matrix_A_1[2, :], simulation_wrapper[:, t-1])
    simulation_wrapper[1, t] = eigenVectorDivision * simulation_wrapper[2, t]
end

# Compute levels
consumption_levels_1 = fill(NaN, T)
capital_levels_1 = fill(NaN, T)

consumption_levels_1[1:6] = simulation_wrapper[1, 1:6] .+ c_star_1
capital_levels_1[1:6] = simulation_wrapper[2, 1:6] .+ k_star_1
# Compute levels for delta = 0.3
consumption_levels_1[6:T] = simulation_wrapper[1, 6:T] .+ c_star_2
capital_levels_1[6:T] = simulation_wrapper[1, 6:T] .+ k_star_2


# PLOTS
con_cap_dev = plot(1:T, simulation_wrapper[1, :], label = "Consumption deviation")
title!("Capital Consumption Deviation")
con_cap_dev_1 = plot!(simulation_wrapper[2, :], label = "Capital deviation")

con_cap_levels = plot(1:T, [consumption_levels_1, capital_levels_1], labels = ["Consumption level" "Capital level"])

# x = [1,2,3,4,5]
# y = [10,20,30,40,50]
# p1 = plot(x,y)
# p2 = scatter(x,y)

# plot(p1,p2,legend=false)

transition_path_capital = plot(capital_levels_1 ./ k_star_old, lw = 3, title = "Transition path capital k", xlabel = "Time period t", ylabel = "Capital k")

transition_path_consumption = plot(consumption_levels_1 ./ c_star_old, lw = 3, title = "Transition path consumption c", xlabel = "Time period t", ylabel = "Consumtion c")

plot(transition_path_capital, transition_path_consumption, legend = false)

# Figure 7 - output
output_levels = capital_levels_1.^alpha
plot(1:T, output_levels)