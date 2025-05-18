# ======= Collection of functions ========
using DifferentialEquations,LinearAlgebra


# =================== Basic Functions # ===================# ===================
# ===================# ===================# ===================# ===================

# Random numbers with fixed sum function
function randfixsum(n, dim, tot)
    sol = []
    for i in 1: n
        v = rand(1,dim)
        nv = tot*v./sum(v)
        push!(sol, nv)
    end
    return sol
end


# function make_cb2(var1, var1_val, var2, var2_val,ts_in)
#     ts = ts_in
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#       if integrator.t == ts[1]
#           integrator.p[var1] = var1_val
#           integrator.p[var2] = var2_val
#       elseif integrator.t == ts[2]
#           integrator.p[var1] = 0.0
#           integrator.p[var2] = 0.0
#       end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     return ts, cb
# end

# function make_cb2(ts_in, vars...)
#     ts = ts_in
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#       if integrator.t == ts[1]
#           integrator.p[vars[1]] = vars[2]
#           integrator.p[vars[3]] = vars[4]
#           integrator.p[vars[5]] = vars[6]
#       elseif integrator.t == ts[2]
#           integrator.p[vars[1]] = 0.0
#           integrator.p[vars[3]] = 0.0
#           integrator.p[vars[5]] = 0.0
#       end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     @show vars[1], vars[2], vars[3], vars[4], vars[5], vars[6]
#     return ts, cb
# end


function make_cb(ts_in, vars...)
    ts = ts_in
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
        for  i = 1:2: length(vars)
            if integrator.t == ts[1]
                integrator.p[vars[i]] = vars[i+1]
            elseif integrator.t == ts[2]
                integrator.p[vars[i]] = 0.0
            end
        end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    @show vars
    return ts, cb
end




# using Reduce
# function subs(old, new, Eq)
#     Algebra.sub(old => new, Eq) |> mat
# end
# function solve(Eq, var)
#     Algebra.solve(:(0  == $Eq),var)
# end






# =================== Stability Test # ===================# ===================
# ===================# ===================# ===================# ===================
function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,p,t)
    return (jac,eigen(jac).values)
end

function EigenVector(Demethy_TF_MC,p,ss)
    J=rand(length(Demethy_TF_MC.syms),length(Demethy_TF_MC.syms))
    t=0
    jacfun(Demethy_TF_MC)(J,ss,p,t)
    return eigvecs(J)
end

function CRN_Eig(rn,params,ss)
    J=rand(length(rn.syms),length(rn.syms))
    t=0
    jacfun(rn)(J,ss,params,t)
    return eigvals(J), eigvecs(J)
end
# ===================# ===================# ===================# ===================




# ===================# ===================# ===================# ===================
# Defineing the Statbility function for checking  SSS
"""
    stability_tianchi(ss, model, param)
# We classify the stability depends on the real and complex part of the Eigen spectrum of the steady state solution(SSS)
1. Caculate The `Eigen spectrum` of the SSS
2. Check the `Real` and `Complex` part of the `Eigen spectrum`
    - SSS is `Real`, check the sign
        - All `Real` are negative, then := `Stable Real`
        - At least one `Real` is positive , then:= `Unstable Real`
    - SSS is `Complex`, check both `real`, and `complex` part
        - All `Real` are negative, then := `Stable Complex`
        - At least one `Real` is positive, then:= `Unstable Complex`

# Color Coding Stability
- `Stable Real` := `blue`
- `Stable Complex` := `cyan`
- `Unstable Complex` := `orange`
- `Unstable Real` := `red`

"""
function stability_tianchi(ss, model, p, con_len)

    function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
        jac = zeros(length(rn.syms),length(rn.syms))
        rn.jac(jac,solution,p,t)
        return (jac,eigen(jac).values)
    end
    # Caculate The Eigen_spectrum of the Steady states solutions
    Eigen_spectrum = [JE_stability(i, model, p)[2] for i in ss]
    # Check the real part of the Eigen_spectrum of SSS with -1e-10 thred
    sb_ind =  [sum(real(i).< -1e-10) .+ con_len for i in Eigen_spectrum]
    # Show the Total number of dimensions that is stable + the conservation law
    # @show sb_ind
    return sb_ind .== length(ss[1])
end























######=======---------------- Macros ---------------------========#####
