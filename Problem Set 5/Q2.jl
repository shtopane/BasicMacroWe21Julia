using Plots; plotly()

# c) Find the rate of return on capital for households for 5 periods
# Parameters
alpha = 0.3
beta = 0.95
sigma = 2
theta = (1 / beta) - 1

k_initial = 2.6257
c_initial = 1.0733
delta_new = 1

# Capital 5 periods
get_capital_for_next_period(alpha, beta, capital_prev) =  capital_prev^alpha - ((1 - alpha*beta)*(capital_prev^alpha))

future_capital_wrapper = fill(NaN, 5)
future_capital_wrapper[1] = k_initial

for i = 1:4
    result = get_capital_for_next_period(alpha, beta, future_capital_wrapper[i])
    future_capital_wrapper[i + 1] = round(result, digits=4)
end

# Rate of return 5 periods
# ror - rate of return
get_ror_for_next_period(alpha, capital_prev) = (alpha * (capital_prev^(alpha - 1))) - 1
future_rate_of_return_wrapper = fill(NaN, 5)

for j = 1:5
    result = get_ror_for_next_period(alpha, future_capital_wrapper[j])
    future_rate_of_return_wrapper[j] = round(result, digits = 4)
end

# plot(future_capital_wrapper, future_rate_of_return_wrapper)

# d) generate wage rage for 5 periods
# wage_rate = f(kt) - f`(kt)*kt
# f(kt) = F(kt, 1) = kt^alpha
# f`(kt) = alphakt^alpha - 1
get_wage_for_next_period(alpha, capital_prev) = (capital_prev^alpha) - ((alpha * (capital_prev^(alpha - 1)) * capital_prev))
future_wage_wrapper = fill(NaN, 5)

for k = 1:5
    result = get_wage_for_next_period(alpha, future_capital_wrapper[k])
    future_wage_wrapper[k] = round(result, digits=4)
end
