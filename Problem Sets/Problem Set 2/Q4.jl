using Plots; plotly()
# Parameters
sigma = 2
alpha = 0.5
A = 1
delta = 0.1
beta = 0.95

# a)

function capital_next(capital_now, consumption_now)
    result = capital_now^alpha + (1 - delta)*capital_now - consumption_now
    return round(result, digits = 3)
end

function consumption_next(capital_next, consumption_now)
    capital_formation = (beta*((alpha * A)*(capital_next^(alpha - 1)) + (1 - delta)))^alpha
    result = consumption_now * capital_formation
    return round(result, digits=3)
end
# consumption_next(capital_next, consumption_now) = consumption_now*(beta*((alpha * A)*(capital_next^(alpha - 1)) + (1 - delta)))^alpha

# b)
# capital_1 = 1
# consumption_1 = 1

# capital_2 = capital_next(capital_1, consumption_1)
# consumption_2 = consumption_next(capital_2, consumption_1)

# struct NegativeCapitalException <: Exception end
# calculate for first 5 periods
function simulate_consumption_capital_path(capital_1, consumption_1, periods_to_simulate)
capital = zeros(1, periods_to_simulate)
consumption = zeros(1, periods_to_simulate)
capital[1] = capital_1
consumption[1] = consumption_1

for i = 1:(periods_to_simulate - 1)
    capital[i + 1] = capital_next(capital[i], consumption[i])

    # stop calculating when capital is negative
    if(capital[i + 1] < 0)
        println("Capital in period $(i + 1)is negative! Value: $(capital[i+1]) ")
        capital[i + 1] = NaN
        consumption[i + 1] = NaN
        #throw(NegativeCapitalException("Capital in period $i + 1 is negative!"))
    end

    consumption[i + 1] = consumption_next(capital[i + 1], consumption[i])
end

# println("capital $capital")
# println("consumption $consumption")

return capital, consumption
end

# capital_1 = 1
# consumption_1 = 1
# capital_period = 4
# consumption_periods = 3
simulate_consumption_capital_path(1,1, 5)

# c) simulate for 3 periods with different starting values
# capital_1 = 1
# consumption_1 = 0.627582223701021
# capital_period = 3
# consumption_periods = 3
simulate_consumption_capital_path(1,0.627582223701021, 3)

# d) - same system as c) but for 30 periods
# capital_1 = 1
# consumption_1 = 0.627582223701021
# capital_period = 30
# consumption_periods = 30
periods = 30
capital_result, consumption_result = simulate_consumption_capital_path(1,0.627582223701021, periods)


plot(capital_result[1:end], consumption_result[1:end], lw=2, label="System $periods periods", color=:coral, markershape=:diamond)
xlabel!("Capital over time")
ylabel!("Consumption over time")
