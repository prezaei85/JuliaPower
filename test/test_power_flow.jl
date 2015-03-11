using DataFrames, Makeps, ACPowerFlow, Optim
C = Constps()

# user input to the code
verbose = true
PartFactFlag = false
PolarFlag = true
require("case6ww_ps")

LoadFactor = 1.0
ps = case6ww_ps()
ps.shunt[:,:P] = ps.shunt[:,:P] * LoadFactor
ps.shunt[:,:Q] = ps.shunt[:,:Q] * LoadFactor

# extract stuff from ps
nBus = size(ps.bus,1)
G = ps.bus_i[ps.gen[:,:bus]]
D = ps.bus_i[ps.shunt[:,:bus]]
Ybus = getYbus(ps)[1]
Vmag = ps.bus[:,:Vmag]
theta = zeros(nBus)
Sbus, Sd, Sg = getSbus(ps)

# Convert loads to simple constant power
Sd = sum(Sd,2)
# Create bus versions of Sd,Sg
Sd_bus = zeros(Complex{Float64},nBus); Sd_bus[D] = Sd
Sg_bus = zeros(Complex{Float64},nBus); Sg_bus[G] = Sg
# Find the bus types
pq  = ps.bus[:,:type].==C.PQ
ref = ps.bus[:,:type].==C.REF
pv  = ps.bus[:,:type].==C.PV
if PartFactFlag
    PartFact = pv | ref
    PartFact = convert(Vector{Float64}, PartFact)
else
    PartFact = Float64[]
end

npq = sum(pq)
if PolarFlag
	# set up x for polar power flow
    ix_theta = 1:(nBus-1);
    ix_Vmag  = (1:npq) + maximum(ix_theta)
    if PartFactFlag
        ix_rho = 1 + maximum(ix_Vmag); # generator ramping variable
        nx = maximum(ix_rho);
    else
        nx = maximum(ix_Vmag);
    end
    x0 = zeros(nx)
    x0[ix_theta] = theta[!ref]
    x0[ix_Vmag]  = Vmag[pq]
    if PartFactFlag
        x0[ix_rho] = 0
    end
    mismatch_func = mismatch_polar
else
	# set up x for rectangular power flow
end
mismatch(x;args...) = mismatch_func(x,Ybus,Vmag,Sg_bus,Sd_bus,pq,pv,ref,PartFact;args...)

## Solve power flow
# newton-raphson
x, flag = nrsolve(mismatch, x0, verbose=true, linesearch="exact")

# solve the unconstrained optimization problem 1/2 * g(x)' * g(x) using levenberg-marquardt 
# algorithm, where g(x) is the mismatch in power injections
g(x) = mismatch(x)[1]
Jac(x) = mismatch(x, NeedMismatch=false, NeedJac=true)[2]
results = Optim.levenberg_marquardt(g, Jac, x0)

