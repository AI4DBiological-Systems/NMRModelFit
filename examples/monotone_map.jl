
### try to analytically or semi-analytically get a, b, such that f(target) = target whilst preserving the domain wrap.

#a_setp, b_setp, _ = NMRModelFit.setupitpab(0.1, 10, 0.7; optim_algorithm = :LN_BOBYQA)
a_setp, b_setp, _ = NMRModelFit.setupitpab(0.1, 10, 0.7;
optim_algorithm = :LN_BOBYQA,
p_ub = [0.9; 5.0])

#
T = Float64
import MonotoneMaps

target = 0.7
#target_01 = target
target_01 = NMRModelFit.convertcompactdomain(target, -one(T), one(T), zero(T), one(T))
a = a_setp(target_01)
b = b_setp(target_01)

function getb(fa, x)#, lb, ub)
    y = (x+1)/2
    return y/(1-y) - (2/(x+1) - 1)^fa(x)
end

b = getb(a_setp, target_01)

targets = LinRange(0,1, 20)
as = a_setp.(targets)
bs = b_setp.(targets)

f = xx->MonotoneMaps.evalcompositelogisticprobit(xx, a, b, -one(T), one(T))
f_last = xx->MonotoneMaps.evalcompositelogisticprobit(xx, as[end], bs[end], -one(T), one(T))
f_first = xx->MonotoneMaps.evalcompositelogisticprobit(xx, as[1], bs[1], -one(T), one(T))
t = LinRange(-1, 1, 500)

# f = xx->MonotoneMaps.evalcompositelogisticprobit(xx, a, b, zero(T), one(T))
# f_last = xx->MonotoneMaps.evalcompositelogisticprobit(xx, as[end], bs[end], zero(T), one(T))
# f_first = xx->MonotoneMaps.evalcompositelogisticprobit(xx, as[1], bs[1], zero(T), one(T))
# t = LinRange(0, 1, 500)




PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

Random.seed!(25)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t, f.(t), label = "f")
PyPlot.plot(t, f_last.(t), label = "f_last")
PyPlot.plot(t, f_first.(t), label = "f_first")

PyPlot.legend()
PyPlot.title("fs vs t")

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(bs, as)
PyPlot.plot(bs, as, "x")

PyPlot.legend()
PyPlot.title("as vs bs")


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(targets, as)
PyPlot.plot(targets, as, "x")

PyPlot.legend()
PyPlot.title("as vs targets")


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(targets, bs)
PyPlot.plot(targets, bs, "x")

PyPlot.legend()
PyPlot.title("bs vs targets")


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(targets, f.(targets))
PyPlot.plot(targets, f.(targets), "x")

PyPlot.legend()
PyPlot.title("f vs targets")
