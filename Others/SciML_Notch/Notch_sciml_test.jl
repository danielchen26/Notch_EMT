## =================================================================
using OrdinaryDiffEq
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, Optim
using DiffEqFlux, Flux
using Plots
using Latexify
gr()
using JLD2, FileIO
using Statistics
# Set a random seed for reproduceable behaviour
using Random
Random.seed!(1234)

## =================================================================
svname = "Notch_EMT_unknown_"
@parameters k0 k1 k2 d m p k pp kk δ α1 A t
# @variables R(t) NR(t) M(t) MR(t)
@variables KDM5A(t) H4(t) H0(t) PRC2(t) H27(t) KDM6A(t) KMT(t)
D = Differential(t)

eqs = [
    # D(R) ~ -A * R + NR - KDM5A * M * R + KDM6A * MR,
    # D(NR) ~ A * R - NR,
    # D(M) ~ -KDM5A * M * R + KDM6A * MR,
    # D(MR) ~ KDM5A * M * R - KDM6A * MR,
    D(KDM5A) ~ k * H27 - δ * KDM5A + α1, # k0 * MR +
    D(H4) ~ -d * H4 * KDM5A + H0 * KMT,
    D(H0) ~ d * H4 * KDM5A - m * H0 * PRC2 + H27 * KDM6A - H0 * KMT,
    D(PRC2) ~ p * H27 - δ * PRC2 + α1,
    D(H27) ~ - H27 * KDM6A, #  m * H0 * PRC2
    D(KDM6A) ~ kk * H4 - δ * KDM6A + α1,
    D(KMT) ~ pp * H4 - δ * KMT + α1]
@named sys = ODESystem(eqs)
latexify(sys)
parameters(sys)

## =================================================================
# Define the experimental parameter
tspan = (0.0f0, 50.0f0)
# u0 = [6.0, 0.0, 6.0, 40.0, 500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]
u0 = [500.0, 0.6876, 0.678, 500.0, 50.6344, 1.0, 2.0]
# p_ = [10.0, 1.0, 1.0, 0.2, 0.53, 1.8, 3.77, 19.08, 19.08, 1.0, 1.0, 0.0]
p_ = [1.0, 1.0, 3.77, 0.2, 0.53, 1.8, 19.08, 19.08]
prob = ODEProblem(sys, u0, tspan, p_, jac = true)
solution = solve(prob, Tsit5(), abstol = 1e-12, reltol = 1e-12, saveat = 0.1)


# Ideal data
X = Array(solution)
t = solution.t


# Add noise in terms of the mean
x̄ = mean(X, dims = 2)
noise_magnitude = Float32(5e-3)
Xₙ = X .+ (noise_magnitude * x̄) .* randn(eltype(X), size(X))

plot(solution, alpha = 0.75, color = :black, label = ["True Data" nothing])
scatter!(t, transpose(Xₙ), color = :red, label = ["Noisy Data" nothing])
plot(solution, vars = [H4, H27, KDM5A, PRC2])

## Define the network
# Gaussian RBF as activation
rbf(x) = exp.(-(x .^ 2))
# Multilayer FeedForward
U = FastChain(
    FastDense(7, 10, rbf), FastDense(10, 10, rbf), FastDense(10, 10, rbf), FastDense(10, 1)
)
# Get the initial parameters
p = initial_params(U)

# Define the hybrid model
function ude_dynamics!(du, u, p, t, p_)
    KDM5A, H4, H0, PRC2, H27, KDM6A, KMT = u
    û = U(u, p) # Network prediction
    du[1] = p_[3] * H27 - p_[2] * KDM5A + p_[1]
    du[2] = -p_[4] * H4 * KDM5A + H0 * KMT
    du[3] = p_[4] * H4 * KDM5A - û[1] + H27 * KDM6A - H0 * KMT
    du[4] = p_[6] * H27 - p_[2] * PRC2 + p_[1]
    du[5] = -H27 * KDM6A#û[2]
    du[6] = p_[7] * H4 - p_[2] * KDM6A + p_[1]
    du[7] = p_[8] * H4 - p_[2] * KMT + p_[1]
end

# Closure with the known parameter
nn_dynamics!(du, u, p, t) = ude_dynamics!(du, u, p, t, p_)
# Define the problem
prob_nn = ODEProblem(nn_dynamics!, Xₙ[:, 1], tspan, p)


## Function to train the network
# Define a predictor
function predict(θ, X = Xₙ[:, 1], T = t)
    Array(solve(prob_nn, Vern7(), u0 = X, p = θ,
        tspan = (T[1], T[end]), saveat = T,
        abstol = 1e-6, reltol = 1e-6,
        sensealg = ForwardDiffSensitivity()
    ))
end

# Simple L2 loss
function loss(θ)
    X̂ = predict(θ)
    sum(abs2, Xₙ .- X̂)
end

# Container to track the losses
losses = Float32[]

# Callback to show the loss during training
callback(θ, l) = begin
    push!(losses, l)
    if length(losses) % 1 == 0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    false
end

## Training

# First train with ADAM for better convergence -> move the parameters into a
# favourable starting positing for BFGS
res1 = DiffEqFlux.sciml_train(loss, p, ADAM(0.2f0), cb = callback, maxiters = 200)
println("Training loss after $(length(losses)) iterations: $(losses[end])")
# Train with BFGS
res2 = DiffEqFlux.sciml_train(loss, res1.minimizer, BFGS(initial_stepnorm = 0.01f0), cb = callback, maxiters = 10000)
println("Final training loss after $(length(losses)) iterations: $(losses[end])")
