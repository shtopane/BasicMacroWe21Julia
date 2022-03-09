using LinearAlgebra
using Plots; plotly()
# e)
# Calculate steady state values: c*, k*
# Parameters
alpha = 0.5
beta = 0.9
theta = (1/beta) - 1

#Steady state values
c_star = (alpha * beta)^(alpha / (1 - alpha)) - (alpha * beta)^(1 / (1 - alpha)) # equation 14
k_star = (1/ (alpha * beta))^ (1/ (alpha - 1))

# Calculate derivatives

# U`(c*)
consumption_first_der = 1/c_star
# U``(c*)
consumption_second_der = - 1 / c_star^2
# F``(k*)
capital_second_der = alpha * (alpha - 1) * k_star^(alpha - 2)

matrix_A = zeros(2,2)
matrix_A[1,1] = 1 + ((beta * (consumption_first_der * capital_second_der)) / consumption_second_der)
matrix_A[2,1] = -1
matrix_A[1,2] = -((consumption_first_der*capital_second_der)/ consumption_second_der)
matrix_A[2,2] = 1 + theta

# f)
# Av = λv

# h) - use CRRA utility function
# U(c) = c¹-σ - 1 / 1 - σ
sigma = 2
delta = 1


# Compute stady state variables
k_star_1 = ( ( 1 - beta + delta * beta) / (alpha * beta))^(1/(alpha - 1)) # F'(k) = alpha*k^(alpha-1) = 1/beta + delta -1
c_star_1 = k_star_1^alpha - delta*k_star_1

# Compute eigenvalues and the linearly approximated consumption policy function
# first term: u`(c*) / u``(c*) = -c*/σ
# second term: F``(k*) = α(1 -α)k*^(α-2)
coeff1 = (-c_star_1/sigma) * (alpha*(alpha - 1)*k_star_1^(alpha - 2)) # u'(c)*F''(K)/u''(c)
matrix_A_1 = zeros(2,2)
matrix_A_1[1,1] = 1 + beta * coeff1
matrix_A_1[1,2] = -coeff1
matrix_A_1[2,1] = -1
matrix_A_1[2,2] = theta + 1

eigenValues, eigenVectors = eigen(matrix_A_1)
stableEigenIndex = 0
stableEigenValue = 0
for i =1:length(eigenValues)
    
    if(eigenValues[i] < 1)
        global stableEigenIndex = i
        global stableEigenValue = eigenValues[i]
    end
end

range_value = round(k_star_1*0.6, digits = 4)
step_value = k_star_1/30
capital_range = collect(-range_value:step_value:range_value) #-k_star_1*0.6/(30:k_star_1)*0.6

# The linearly approximated consumption policy function
eigenVectorDivision = (eigenVectors[1, stableEigenIndex] / eigenVectors[2, stableEigenIndex])
consumption_policy_function_evaluated = eigenVectorDivision*capital_range

linearly_appr_cons_fun_plot = plot(capital_range, consumption_policy_function_evaluated)
title!("Linearly approximated consumption policy function")
xlabel!("Capital today, deviation from st. st.")
ylabel!("Consumption today, deviation from st. st.")


# Calculate the transition path for consumption and capital starting from the same starting value
T = 50
simulation_storage = zeros(2, T)

# Starting values
simulation_storage[2,1] = -0.1*k_star_1
simulation_storage[1,1] = eigenVectorDivision*simulation_storage[2,1]

for t = 2:T
    simulation_storage[2, t] = dot(matrix_A_1[2, :],simulation_storage[:,t-1])
    simulation_storage[1, t] = eigenVectorDivision * simulation_storage[2, t]
end

# deviation_plot = plot(1:T, simulation_storage)
# display(deviation_plot)
capital_con_deviation_sim_plot = plot(1:T, simulation_storage[1, :], label="Consumption deviation")
plot!(simulation_storage[2, :], label="Capital deviation")
xlabel!("Time")

cap_con_dev_relative_plot = plot(simulation_storage[2, :], simulation_storage[1, :], markershape=:diamond, color=:inferno)
ylabel!("Consumption deviation")
xlabel!("Capital deviation")