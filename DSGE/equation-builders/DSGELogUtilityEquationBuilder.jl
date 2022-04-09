module DSGELogUtilityEquationBuilder

# k̄ = (1/β - 1 + δ) ^ 1/ α - 1
#           α
function getSteadyStateCapital(alpha, beta, delta)
    (((1 / beta) - 1 + delta) / alpha)^(1 / (alpha - 1))
end
# c̄ = s̄k̄ᵅ  δk̄
function getSteadyStateConsumption(k_bar, alpha, delta, s_bar)
    s_bar * (k_bar^alpha) - delta * k_bar
end

# constant: s̄ = 1
function getSteadyStateProductivity(s_bar_constant)
    s_bar_constant
end

# Formulas for a1,a2,b1,b2,b3 coefficients
function getLinearCoefficients(alpha, beta, delta, k_bar, c_bar)
    a1 = beta * alpha * (alpha - 1) * k_bar^(alpha - 1)
    a2 = beta * alpha * k_bar^(alpha - 1)
    b1 = 1 - delta + alpha * k_bar^(alpha - 1)
    b2 = k_bar^(alpha - 1)
    b3 = -c_bar/k_bar

    return (round(a1, digits=4), round(a2, digits=4), round(b1, digits=4), round(b2, digits=4), round(b3, digits=4))
end

export getSteadyStateCapital,getSteadyStateConsumption, getSteadyStateProductivity
end