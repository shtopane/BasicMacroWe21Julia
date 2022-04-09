a = 1 + 1
alpha = 0.3
beta = 0.95
delta = 0.10
theta = (1 / beta) - 1

beta_t = 0.91
e1 = round((1/beta_t) - 1, digits=4)
k_star = ((1 - alpha) * (theta + delta^(alpha/(alpha - 1))) * alpha^(1/(1-alpha)))/theta + ((1-alpha)*delta)