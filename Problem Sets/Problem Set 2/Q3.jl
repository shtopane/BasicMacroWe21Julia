# a) calculate k*
# Parameters
sigma = 2
alpha = 0.5
beta = 0.95
delta = 0.1
A = 1
theta = (1/beta) - 1

# Functions
cobb_douglass(capital::Float64) = A*capital^alpha

k_star = ((delta + theta) / (alpha * A))^(1/(alpha - 1))
c_star = cobb_douglass(k_star) - (delta * k_star)

# b) evaluate derivatives for: u`(c*); u``(c*) F``(k*)
utility_first_der(consumption::Float64) = consumption^-sigma
utility_second_der(consumption::Float64) = -sigma*consumption^(-sigma-1)

consumption1 = utility_first_der(c_star)
consumption2 = utility_second_der(c_star)

cobb_douglass_first_der(capital::Float64) = alpha * A*capital^(alpha - 1)
cobb_douglass_second_der(capital::Float64) = alpha *(alpha - 1)*A*capital^(alpha - 2)

capital1 = cobb_douglass_first_der(k_star)
capital2 = cobb_douglass_second_der(k_star)

# c) - calculate deviations from steady-state values
deviation_period_start = reshape([1,2], :, 1)
matrixA = [[1.0073, -1] [-0.0078, 1.0526]]
# Setup for which period we want deviations
period_end = 19
period_start = 17
power = period_end - period_start

deviation_period_end = matrixA^power * deviation_period_start
println("deviation_period_end $deviation_period_end" )