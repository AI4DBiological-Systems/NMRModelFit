

using LinearAlgebra, FFTW
import BSON, Statistics, Random

import NMRSignalSimulator

include("../src/NMRModelFit.jl")
import .NMRModelFit

import PlotlyJS
using Plots; plotly()

# for loading something with Interpolations.jl
import OffsetArrays
import Interpolations

import MultistartOptimization, NLopt


include("helper.jl")

# run prep_full.jl first.

r = 1
T = Float64
N_starts = 500
local_optim_algorithm = NLopt.LN_BOBYQA
xtol_rel = 1e-9
maxeval = 200 # 2 # 50
maxtime = Inf
β_optim_algorithm = :GN_DIRECT_L
w_lb_default = 1e-1
w_ub_default = 100.0
β_max_iters = 500 # 2 # 500
β_xtol_rel = 1e-9
β_ftol_rel = 1e-9
β_maxtime = Inf


N_d = sum( NMRModelFit.getNd(Bs[n]) for n = 1:length(Bs) )
#N_β = sum( NMRModelFit.getNβ(κs_β_DOFs[n], Bs[n]) for n = 1:length(Bs) )

shift_lb = -ones(N_d)
shift_ub = ones(N_d)

a_setp, b_setp, minxs,
    rets = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
#

y_cost = y[cost_inds_set[r]]
U_cost = U_y[cost_inds_set[r]]
P_cost = P_y[cost_inds_set[r]]
LS_inds = 1:length(U_cost)

# prepare.
#N_d = sum( NMRModelFit.getNd(Bs[n]) for n = 1:length(Bs) )
@assert length(shift_ub) == length(shift_lb) == N_d

U_rad_cost = U_cost .* (2*π)

# setup inner optim over β.
q0, updatedfunc, getshiftfunc, N_vars_set,
run_optim, obj_func_β, E_BLS, w_BLS, b_BLS, updateβfunc, updatewfunc,
q_β = NMRModelFit.setupcostnesteddwarpw(Bs, As, fs, SW, LS_inds, U_rad_cost,
    y_cost, Δsys_cs, a_setp, b_setp;
    β_optim_algorithm = β_optim_algorithm,
    w_lb_default = w_lb_default,
    w_ub_default = w_ub_default,
    β_max_iters = β_max_iters,
    β_xtol_rel = β_xtol_rel,
    β_ftol_rel = β_ftol_rel,
    β_maxtime = β_maxtime)

# set up outer optim over shifts.
N_β = sum( NMRModelFit.getNβ(Bs[n]) for n = 1:length(Bs) )
p_β = zeros(T, N_β) # persistant buffer.

obj_func = pp->NMRModelFit.costnesteddw(U_rad_cost, y_cost, updatedfunc,
updatewfunc, pp,  Bs, run_optim,
E_BLS, w_BLS, b_BLS, p_β)

p = ones(N_d) .* 0.02
#p = [0.02;] # Serine 700 MHz BMRB.
#p = [0.13513301735885624;]
updatedfunc(p)

fill!(p_β, zero(T)) # always start from 0-phase?
minf, minx, ret, N_evals = run_optim(p_β)
p_β[:] = minx # take out?

#p_β[2] = p_β[1]

#NMRModelFit.updateβ!(Bs, κs_β_orderings, κs_β_DOFs, p_β, 1)
NMRModelFit.updateβ!(Bs, p_β, 1)
updatewfunc(1.0)

cost = norm(reinterpret(T, E_BLS)*w_BLS - b_BLS)^2

#### visualize.
N_viz = 50000

U = LinRange(u_min, u_max, N_viz)
P = hz2ppmfunc.(U)
U_rad = U .* (2*π)

q2 = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs; w = w_BLS )
q_U = q2.(U_rad)

plotregion(P, U, q_U, P_y, y, P_cost, y_cost, 0.0, 1,
    "", "test", "";
    save_plot_flag = false,
    display_plot_flag = true,
    canvas_size = (1000,400))

##revisit cost again.

function evalcost(q2, U_cost_rad, y_cost)
    cost2 = 0.0
    for m = 1:length(U_cost_rad)
        cost2 += abs2(q2(U_cost_rad[m]) - y_cost[m])
    end

    return cost2
end

U_cost_rad = U_cost .* (2*π)
cost2 = evalcost(q2, U_cost_rad, y_cost)

#norm(b_BLS - reinterpret(T, y_cost))
@assert 12==3

β_lb = ones(T, N_β) .* (-π)
β_ub = ones(T, N_β) .* (π)

p_lb = β_lb
p_ub = β_ub

q3, updateβfunc3, updatewfunc3, E_BLS3, w_BLS3, b_BLS3,
getβfunc3 = NMRModelFit.setupcostβLSw(Bs, As, LS_inds, U_rad_cost, y_cost;
    w_lb_default = w_lb_default,
    w_ub_default = w_ub_default)
#
updatewfunc3(1.0)
# E_BLS3*w_BLS3 - E_BLS*w_BLS
# norm(E_BLS3- E_BLS)
#
# (E_BLS3*w_BLS3)[1]
# u_rad = U_cost[1]*2*pi
# println("q0(u_rad) = ", q0(u_rad))
# println("q3(u_rad) = ", q3(u_rad))

# q is simulated spectra. f is cost function.
f = pp->NMRModelFit.costβLSw(U_rad_cost, y_cost, updateβfunc, updatewfunc, pp,
E_BLS, w_BLS, b_BLS, q3)
