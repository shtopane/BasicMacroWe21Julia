# This is the code for Problem_Set_1, Question 5 code.
# Macroeconomics, Hintermaier, Winter 2021/2022, Bonn
# Production function: 
# (αktˢ + (1 - α)L̄ˢ)\^(1/s) 

using Plots; plotly()
# Parameters
alpha = 1/3
delta = 0.5
Lbar = 1 # = L - some constant
s = 0.1

T = 1000 # periods simulated
initial_multiple = 1/2 # Situation (b)

# Golden rule steady state values - from (a)
ksharp = ((delta^(s/(1-s)) - alpha^(1/(1-s)))/ ((1-alpha) * alpha^(s/(1-s)) * Lbar^s)) ^(-1/s) # F`(k) = δ, equation (19), page 4
csharp = (alpha*ksharp^s+(1-alpha)*Lbar^s)^(1/s) - delta * ksharp # equation (22), page 4: c_sharp = F(k_sharp) - δk_sharp

kgrid = LinRange(0, 2*ksharp, T)

prod_today = (alpha * kgrid.^s .+ (1-alpha)*Lbar^s).^(1/s)
kplus = prod_today .+ (1 - delta)*kgrid .- csharp # k_t+1 = F(k) + (1-δ)*k - c_sharp; equation (23)

plot(kgrid, kplus, lw=2)
plot!(kgrid,kgrid, color="red", lw=2)

min_kgrid = minimum(kgrid)
max_kgrid = maximum(kgrid)
kgrid_lims = (min_kgrid, max_kgrid)
plot!(xlims = kgrid_lims, ylims=kgrid_lims)
plot!(xlims = (minimum(kgrid), maximum(kgrid)), ylims=(minimum(kgrid), maximum(kgrid)))

scatter!((ksharp, ksharp), color="pink", markersize = 5)
xlabel!("Capital today")
ylabel!("Capital tomorrow")
title!("Dynamics under the Golden Rule")

# Simulating Dynamics
k_t = fill(NaN, T)
k_t[1] = ksharp*initial_multiple

for i = 1:(T - 1)
    k_t[i + 1] = (alpha* k_t[i]^s + (1-alpha)*Lbar^s)^(1/s)+ (1 - delta)*k_t[i] - csharp # k_t+1 = F(k) + (1-delta)*k - c_sharp

    if(k_t[i + 1] < 0)
        k_t[i + 1] = NaN
        break
    end
end

# Plotting
# plot!(k_t[1:end - 1], k_t[2:end], color="green")
scatter!(k_t[1:end - 1], k_t[2:end], markersize=8, markershape=:diamond, color=:blue )
xlabel!("Capital today")
ylabel!("Capital tomorrow")

# Filter out the NaNs and plot the capital over time
capital_with_values = filter(capital -> capital !== NaN, k_t)
time = 1:length(capital_with_values)

scatter(time, capital_with_values, markershape=:diamond, color=:blue, label="Capital")
xlabel!("Time")
ylabel!("Capital")


# Problem set 1, Question 5 with Cobb-Douglass Prod function
# F(kt) = Aktᵅ
cobb_alpha = 0.5
cobb_A = 1
cobb_beta = 0.95
cobb_delta = 0.1

cobb_ksharp = (cobb_delta / (cobb_alpha*cobb_A))^(1/cobb_alpha-1) # αAkᵅ-1 = δ, ksharp = (δ/ αA)^1/ α - 1
cobb_csharp = (cobb_A * cobb_ksharp^cobb_alpha) - (cobb_delta * cobb_ksharp) # c_shrap = F(k_shrap) - delta*k_sharp

cobb_kgrid = LinRange(0, 2*cobb_ksharp, T)

cobb_kplus = (cobb_A * cobb_kgrid.^cobb_alpha) .+ (1-cobb_delta).*cobb_kgrid .- cobb_csharp# k_t+1 = F(kgrid) + (1-δ)*kgrid - c_sharp;
plot(cobb_kgrid, cobb_kplus, lw=3)
plot!(cobb_kgrid, cobb_kgrid, color=:red, lw=3)

min_cobb_kgrid = minimum(cobb_kgrid)
max_cobb_kgrid = maximum(cobb_kgrid)
cobb_kgrid_lims = (min_cobb_kgrid, max_cobb_kgrid)
plot!(xlims = cobb_kgrid_lims, ylims=cobb_kgrid_lims)

scatter!((cobb_ksharp, cobb_ksharp), color="pink", markersize = 5)
xlabel!("Capital today(Cob)")
ylabel!("Capital tomorrow(Cob)")
title!("Dynamics under the Golden Rule(Cob)")

# Simulate
cobb_k_t = fill(NaN,T)
cobb_k_t[1] = cobb_ksharp * initial_multiple

for i = 1:(T-1)
    cobb_k_t[i+1] = (cobb_A * cobb_k_t[i].^cobb_alpha) .+ (1-cobb_delta).*cobb_k_t[i] .- cobb_csharp# k_t+1 = F(kgrid) + (1-δ)*kgrid - c_sharp;

    if(cobb_k_t[i + 1] < 0)
        cobb_k_t[i + 1] = NaN
        break
    end

end

scatter!(cobb_k_t[1:end - 1], cobb_k_t[2:end], markersize=8, markershape=:diamond, color=:blue )
# Filter out the NaNs and plot the capital over time
capital_with_values = filter(capital -> capital !== NaN, cobb_k_t)
time = 1:length(capital_with_values)

scatter(time, capital_with_values, markershape=:diamond, color=:blue, label="Capital")
xlabel!("Time")
ylabel!("Capital")