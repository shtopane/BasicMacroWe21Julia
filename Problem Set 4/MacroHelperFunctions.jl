module MacroHelperFunctions
function getSteadyStateConsumption(steady_state_capital, alpha, delta)
    (steady_state_capital^alpha) - (delta * steady_state_capital)
end

function getSteadyStateCapital(alpha, delta, theta)
    (alpha / (delta + theta))^(1 / (1 - alpha)) # F'(k) = alpha*k^(alpha-1) = 1/beta + delta -1 
end

export getSteadyStateCapital, getSteadyStateConsumption
end