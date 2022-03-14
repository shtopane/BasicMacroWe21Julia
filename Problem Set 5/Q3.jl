# linearized solution of a Ramsey model with an unaticipated and 
# permanent shock to depreciation
# CRRA utility function
# u(c) = (c^(1 - \sigma) - 1)/(1 - \sigma),
# Cobb-Douglas production function
using LinearAlgebra
using Plots;
plotly();

include("../Problem Set 4/MacroHelperFunctions.jl")
using .MacroHelperFunctions

# Parameters
alpha = 0.3
beta = 0.95
sigma = 2.0
delta_old = 0.1
delta_new = 0.2
Ā = 1

T = 50

# INITIAL Steady state
k_star_old = getSteadyStateCapital(alpha, delta_old, theta)
c_star_old = getSteadyStateConsumption(k_star_old, alpha, delta_old)

rate_of_return_assets_old = 1 / beta - 1 # steady_state + Euler = beta * U(c)*(1 + r) => 1/beta = 1 + r
rental_price_capital_old = rate_of_return_assets_old + delta_old # households face depreciation imposed on their assets?
# wage_rate_old = (k_star_old^alpha) - ((alpha * (k_star_old^(alpha - 1)) * k_star_old)) # w = f(k) - f'(k)*k
wage_rate_old = getWageRate(alpha, k_star_old) # w = f(k) - f'(k)*k

# NEW Steady state
k_star_new = getSteadyStateCapital(alpha, delta_new, theta)
c_star_new = getSteadyStateConsumption(k_star_new, alpha, delta_new)
rate_of_return_assets_new = 1 / beta - 1 # SAME as delta_old
rental_price_capital_new = rate_of_return_assets_new + delta_new
wage_rate_new = getWageRate(alpha, k_star_new)

# Compite deviation of capital - needed for policy function
capital_deviation = k_star_old - k_star_new

# Get matrix entries for NEW regime
# U`(c*)
consumption_first_der = 1 / (c_star_new^sigma)
# U``(c*)
consumption_second_der = -sigma * (c_star_new^(-sigma - 1))
# F``(k*)
capital_second_der = alpha * (alpha - 1) * k_star_new^(alpha - 2)
matrix_A = zeros(2, 2)
matrix_A[1, 1] = 1 + beta * (consumption_first_der * capital_second_der) / consumption_second_der
matrix_A[2, 1] = -1
matrix_A[1, 2] = -((consumption_first_der * capital_second_der) / consumption_second_der)
matrix_A[2, 2] = 1 + theta

eigenValues, eigenVectors = eigen(matrix_A)
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
simulation_wrapper[2, 1] = capital_deviation
simulation_wrapper[1, 1] = eigenVectorDivision * simulation_wrapper[2, 1]

for t = 2:T
    simulation_wrapper[2, t] = dot(matrix_A[2, :], simulation_wrapper[:, t-1])
    simulation_wrapper[1, t] = eigenVectorDivision * simulation_wrapper[2, t]
end

# Compute levels
consumption_levels = simulation_wrapper[1, :] .+ c_star_new
capital_levels = simulation_wrapper[2, :] .+ k_star_new

# Compute prices?
# Levels.
rental_price_capital_levels = alpha .* (capital_levels .^ (alpha - 1)) # Rental price of capital, R = F`(k)
real_interest_rate_levels = rental_price_capital_levels .- delta_new # Real interest rate?
wage_rate_levels = [getWageRate(alpha, capital) for capital in capital_levels] # wage rate

# Compute deviations from new steady state
rental_price_capital_deviation = rental_price_capital_levels .- rental_price_capital_new
real_interest_rate_deviation = real_interest_rate_levels .- rate_of_return_assets_new
wage_rate_deviation = wage_rate_levels .- wage_rate_new

# Labor income share
labor_share = wage_rate_levels ./ (capital_levels .^ alpha)


# PLOTS
con_cap_dev = plot(1:T, simulation_wrapper[1, :], label = "Consumption deviation")
title!("Capital Consumption Deviation")
con_cap_dev_1 = plot!(simulation_wrapper[2, :], label = "Capital deviation")

con_cap_levels = plot(1:T, [consumption_levels, capital_levels], labels = ["Consumption level" "Capital level"])

# x = [1,2,3,4,5]
# y = [10,20,30,40,50]
# p1 = plot(x,y)
# p2 = scatter(x,y)

# plot(p1,p2,legend=false)

transition_path_capital = plot(capital_levels ./ k_star_old, lw = 3, title = "Transition path capital k", xlabel = "Time period t", ylabel = "Capital k")

transition_path_consumption = plot(consumption_levels ./ c_star_old, lw = 3, title = "Transition path consumption c", xlabel = "Time period t", ylabel = "Consumtion c")

plot(transition_path_capital, transition_path_consumption, legend = false)

# Figure 4
real_interest_rate_price_of_capital_level_plot = plot(1:T, [real_interest_rate_levels, rental_price_capital_levels], lw = 3, xlabel = "Time", ylabel = "Levels", label = ["Real interest rate, level" "Rental price of capital, level"], color = [:red :coral])
wage_rate_level_plot = plot(1:T, wage_rate_levels, lw = 3, xlabel = "Time", ylabel = "Levels", label = "Wage rate, level")
plot(
    real_interest_rate_price_of_capital_level_plot,
    wage_rate_level_plot
)

# Figure 5
real_interest_rate_price_of_capital_deviation_plot = plot(
    1:T, [real_interest_rate_deviation, rental_price_capital_deviation],
    label = ["Real interest rate, deviation" "Rental price of capital, deviation"],
    color = [:black :cyan]
)

wage_rate_deviation_plot = plot(
    1:T, wage_rate_deviation,
    xlabel = "time",
    ylabel = "Deviations from new w*",
    label = "Wage rate"
)

plot(
    real_interest_rate_price_of_capital_deviation_plot,
    wage_rate_deviation_plot
)

# Figure 6 - Not need for figure - labor share, which is constant but okay
plot(1:T, [round(x, digits=3) for x in labor_share], ylims=(0.0,1.0), label="Labor share", color=:pink, lw=3)

# Figure 7 - output
output_levels = capital_levels.^alpha
plot(1:T, output_levels)