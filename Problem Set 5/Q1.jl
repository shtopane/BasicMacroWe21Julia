# Problem set 5, Question 1
# Production function: F(kt, nt) = kₜᵅnₜᵅ-1
# Utility function: U(ct, lt) = log(ct) - ϕ(1 - lt), ϕ > 0

# b) find steady state values of k*, c*, n*, l*

# Parameters
alpha = 0.3
beta = 0.95
delta = 0.1
theta = (1 / beta) - 1
phi = 0.01 # TEST to derive steady state values

# Using hint that Fk = α yt / kt (8)
# and Fn = (1 - α) yt/nt (9)
# prod_capital_partial(alpha, output, capital) = alpha * (output / capital) # (8)
# prod_labor_partial(alpha, output, labor) = (1 - alpha) * (output / labor) # (9)

# 1. Impose steady state in the Euler equation:
# α*kstar^α-1 * nstar^1 - α = δ + Θ (10)
# 1.1 Using (8) with (10)
# α* ystar/kstar = δ + Θ (11)
output_to_capital_ratio(theta, delta, alpha) = (theta + delta) / alpha

# 2. Impose steady state on labor to find
# 2.1 ϕ = 1/cstar * (1-α)*kstar^α*nstar^-α (12)
# 2.2 Use (9) and (12) to find
# ϕ = 1/cstar * (1-α) * ystar/nstar

# 3. Imposing steady state on the dynamic resource constraint we get:
# 3.1 cstar = ystar - δ*kstar (14)
# 3.2 cstar / kstar = ystar/kstar - δ (15)
# 3.3 From (11) we use ystar/kstar in (15) to find
#     cstar/kstar = Θ + (1 - α)δ /α (16) !! This is computable
consumption_to_capital_ratio(theta, alpha, delta) = (theta + (1 - alpha) * delta)/alpha
consumption_to_capital = consumption_to_capital_ratio(theta, alpha, delta)

# 3.4 From (16) and (11) we can find
# cstar/ystar = Θ + (1 - α)δ / Θ + δ (17)
# ⇑ alternative consumption_to_capital_ratio / output_to_capital_ratio
consumption_to_output_ratio(theta, alpha, delta) = (theta + (1 - alpha) * delta) / (theta + delta)
consumption_to_output1 = consumption_to_output_ratio(theta, alpha, delta)
# consumption_to_output2 = consumption_to_capital_ratio(theta, alpha, delta) / output_to_capital_ratio(theta, delta, alpha) # Same as the function output

# 3.5 Obtain n*: from (17)
# nstar = 1/ϕ * ystar/cstar* (1 - α)
# nstar = (Θ + δ)(1 - α) / ϕ(Θ + (1 - α)δ) (18)
nstar = ((theta + delta)*(1 - alpha)) / (phi * (theta + (1 - alpha)*delta))
# 3.6 Obtain l*: l* = 1 - n*
# lstar = 1 - (Θ + δ)(1 - α) / ϕ(Θ + (1 - α)δ) (19)
lstar = 1 - nstar
# 3.7 Combine (18) with (10) to get k*
# kstar = (1 - α)(Θ + δ)^α/α-1 * α^1/1 - α / Θϕ + (1 - α)δϕ (20)
kstar = ((1 - alpha)*((theta + delta)^(alpha/(alpha-1)) * (alpha^(1/(1 - alpha))))) / ((theta*phi) + ((1 - alpha)*delta*phi))
# 3.8 Obtain c*
# cstar = 1/ϕ(1 - α)(Θ + δ)^1/α - 1 * α^\1/1 - α (21)
cstar = 1/phi * (1 - alpha) * ((theta + delta)^1/(alpha - 1)) * (alpha^1/(1 - alpha))