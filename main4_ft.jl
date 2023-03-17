using Plots
using IntervalArithmetic
using IntervalRootFinding: roots

##

function phase_π(; n, γ, ω0, ω, t)
    # this defines the total phase is how many times of π
    phase_value = n * ω * t .+ atan(γ * n * ω / (ω0^2 - (n * ω)^2)^2) / π
    return phase_value
end


##
plt = plot()
for n ∈ 1:2:9
    tt = 0:0.1:10
    phase_value = n * ω * t .+ atan(γ * n * ω / (ω0^2 - (n * ω)^2)^2) / 2π
    @show phase_value
    phase_π(; n = 1, γ = 1, ω0 = 1, ω = 1, t = tt)
    plt = plot!(tt, phase_value, label = "$n th term", xlabel = "t", ylabel = "phase")
end
plt

# =======  test for t = 0.1 where, n= 1 term of x(t)_1 phase is π =================
t1 = phase_π(; n = 1, γ = 1, ω0 = 1, ω = 1, t = 0.5)
t3 = phase_π(; n = 3, γ = 1, ω0 = 1, ω = 1, t = 0.5)
t5 = phase_π(; n = 5, γ = 1, ω0 = 1, ω = 1, t = 0.5)
t7 = phase_π(; n = 7, γ = 1, ω0 = 1, ω = 1, t = 0.5)


##
function find_nterms_phase_distribution(; multi_π = 2, γ = 1, ω0 = 12, ω = 1, time_range = 0..100)
    phaset(t) = phase_π(; n = 1, γ = γ, ω0 = ω0, ω = ω, t) - multi_π
    N1_max_t = roots(phaset, time_range)
    N1_max_t = first.(mid.(interval.(N1_max_t)))
    @show N1_max_t
    n_term = []
    n_phase = []
    for n = 1:2:9
        @show n
        append!(n_term, n)
        N_phaset = phase_π(; n = n, γ = γ, ω0 = ω0, ω = ω, t = N1_max_t)
        @show N_phaset
        append!(n_phase, N_phaset)
    end
    return N1_max_t, n_term, n_phase
end




## ===================== when 1st term phase is 2π, the next few terms phase distribution =====================
anim = @animate for ω ∈ 0.1:0.1:10
    N1_max_t, n_term, n_phase = find_nterms_phase_distribution(multi_π = 2, γ = 0.6, ω0 = 3, ω = ω)
    plot(n_phase * π, n_term, proj = :polar,
        m = :red, label = "$ω")
end
gif(anim, fps = 2)





## ===================== when 1st term phase is π, the next few terms phase distribution =====================
