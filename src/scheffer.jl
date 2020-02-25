using DifferentialEquations
using Plots

# Parameters
const ib = 2e-5; # g⋅m⁻²⋅day⁻¹ immigration rate of bream
const ip = 2e-5; # g⋅m⁻²⋅day⁻¹ immigration rate of pike
const r = 7.5e-3; # day⁻¹ maximum growth rate of bream
const H₁ = 0.5; # half saturation constant
const H₂ = 10; # % half saturation constant
const H₃ = 20; # g⋅m⁻² half saturation constant
const H₄ = 15; # g⋅m⁻² half saturation constant
const cb = 7.5e-5; # m⁻²⋅g⁻¹⋅day⁻¹ intraspecific competition constant for bream
const cp = 2.75e-4; # m⁻²⋅g⁻¹⋅day⁻¹ intraspecific competition constant for pike
const prmax = 5e-2; # day⁻¹ maximum predation rate of pike
const ce = 0.1; # pike food conversion efficiency to growth
const mp = 2.25e-3; # day⁻¹ mortality rate of pike
const K = 100; # % maximum vetetation coverage

function scheffer!(du,u,p,t)
    B, P = u # g⋅m⁻² bream/pike density
    ib, r, nutr, H₁, H₂, H₃, H₄, cb, prmax, ip, ce, mp, cp, K = p

    V = K*(H₃^2/(H₃^2 + B^2)); # % of lake covered vetetation
    FR = B^2/(B^2 + H₄^2); # fuctional response of pike

    du[1] = dB = ib + r*(nutr/(nutr + H₁))*B - cb*B^2 - prmax*FR*P
    du[2] = dP = ip + ce*prmax*FR*P*(V/(V + H₂)) - mp*P - cp*P^2
end

function run_ss_solver(;
             nutr = 0.5, #nutrient level
             u0 = [40.0, 2.6], #Initial populations
             p = [ib, r, nutr, H₁, H₂, H₃, H₄, cb, prmax, ip, ce, mp, cp, K])
    prob = SteadyStateProblem(scheffer!,u0,p);
    solve(prob)
end

equilibria = Set(); 
for i in range(0.0, stop=100.0, length=10)
    for j in range(0.0, stop=10.0, length=10)
        for k in [0.5, 1.5, 2.5];
            eq = run_ss_solver(u0=[i,j], nutr=k);
            push!(equilibria, (eq[1], eq[2]));
        end
    end
end

function dB(B, P)
    FR = B^2/(B^2 + H₄^2); # fuctional response of pike

    ib + r*(nutr/(nutr + H₁))*B - cb*B^2 - prmax*FR*P
end

function dP(B, P)
    V = K*(H₃^2/(H₃^2 + B^2)); # % of lake covered vetetation
    FR = B^2/(B^2 + H₄^2); # fuctional response of pike

    ip + ce*prmax*FR*P*(V/(V + H₂)) - mp*P - cp*P^2
end

bream = range(0.0, stop=100.0, length=100)
pike = range(0.0, stop=10.0, length=100)

nutr = 2.5;
contour(bream, pike, dB, levels=[0], color=:blue, cbar=false, legend=true, label="bream'=0, nutr=2.5");
nutr = 1.5;
contour!(bream, pike, dB, levels=[0], color=:red, cbar=false, label="bream'=0, nutr=1.5");
nutr = 0.5;
contour!(bream, pike, dB, levels=[0], color=:green, cbar=false, label="bream'=0, nutr=0.5");
contour!(bream, pike, dP, levels=[0], color=:black, cbar=false, label="pike'=0");

eq = collect(equilibria);
scatter!(first.(eq), last.(eq), label="equilibria");

gui()

