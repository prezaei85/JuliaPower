module ACSIMSEP

using Makeps, ACPowerFlow

export get_left_eigv,
	   find_Sm,
       find_sens_fact

function get_left_eigv(J::SparseMatrixCSC)
    # gives the left eigen vector of the zero eigen value for the jacobian and checks for some errors
    EPS = 0.03
    if size(J)[1] > 2
        # right now, eigs works only when number of rows (cols) is greater than 2. 
        # Issue reported. Will be fixed in future versions!
        w_= eigs(J', nev=1, which=:SM) # gives the eigen value/vector for smallest 
                                          # (magnitude) eigen value
        if abs(imag(w_[1])[1]) > EPS || abs(real(w_[1])[1]) > EPS
            println("Check this case. The smallest eigen value has a large real/imaginary value.")
        end
        if maximum(abs(imag(w_[2]))) > EPS
            println("Check this case. Eigen vector has imaginary values.")
        else
            w = real(w_[2])
        end
    else
        w_ = eig(full(J'))
        idx = abs(w_[1]) .< EPS
        if sum(idx) == 0
            println("Check this case. No zero eigen values found.")
        elseif sum(idx) > 1
            println("Check this case. More than one zero eigen values found.")
        end
        w = w_[2][:,idx]
    end
    return w
end

function find_Sm(mismatch::Function, x0::Vector{Float64}; verbose::Bool = true)
    # This function finds the minimum distance in parameter space from 
    # an unsolvable point to the solvable region (based on Overbye's paper)
    x, converged = nrsolve(mismatch, x0, verbose=verbose, linesearch="exact")
    if converged == true
        dist = 0
        println("The initial power flow converged.")
    end
    g(x) = mismatch(x)[1]
    Jac(x) = mismatch(x, NeedJac=true)[2]
    S = mismatch(x0)[3]
    # find the left eigen vector of jacobian
    w = get_left_eigv(Jac(x))
    F(x) = 0.5*dot(g(x),g(x))
    F_new = F(x)
    if verbose
        println("F(x): $F_new")
    end
    F_prev = 0
    while abs(F_new - F_prev) > 1e-4
        S_i = S + sum(g(x).*w)*w
        mismatch_new(x;args...) = mismatch(x, Si=vec(S_i);args...)
        x, converged = nrsolve(mismatch_new, x0, verbose=false, linesearch="exact")
        F_prev = copy(F_new)
        F_new = F(x)
        if F_new > F_prev
            println("WARNING: The objective function has increased.")
        end
        w = get_left_eigv(Jac(x))
        if verbose
            println("F(x): $F_new")
        end
    end
    # find the minimum distance
    beta = sqrt(2*F_new)
    return x, Jac(x), beta
end

function find_sens_fact(ps::Caseps, PartFactFlag; verbose = true)
    mismatch, x0, Vmag, pq, ref, Yf, Yt, ix_Vmag, ix_theta = init_mismatch(ps, true, PartFactFlag, true);
    # find the minimum distance to solvability
    xm, Jm, beta = find_Sm(mismatch,x0,verbose = verbose);
    # get the left eigen vector of Jacobian associated with zero eigen value 
    w_m = get_left_eigv(Jm)
    # resolve the direction of w_m (it needs to be outward the solvable region)
    g(x) = mismatch(x)[1]
    if dot(g(xm),w_m) < 0
        w_m = -w_m
    end
    nbus = size(ps.bus,1)
    Sbus, = getSbus(ps)
    if PartFactFlag
        # beta_P: the sensitivity factors when reducing real power at buses
        beta_P = w_m[1:nbus]
        
        # beta_Q: the sensitivity factors when reducing reactive power at buses
        beta_Q = zeros(nbus)
        beta_Q[pq] = w_m[nbus+1:end]
    else
        # beta_P: the sensitivity factors when reducing real power at buses
        beta_P = zeros(nbus)
        beta_P[!ref] = w_m[1:nbus-1]
        
        # beta_Q: the sensitivity factors when reducing reactive power at buses
        beta_Q = zeros(nbus)
        beta_Q[pq] = w_m[nbus:end]
    end
    # beta_S: the sensitivity factors when reducing complex power at buses while maintaining power factor
    # only considers buses with positive load
    Pbus = real(Sbus)
    Qbus = imag(Sbus)
    LoadBusID = (Pbus.<=0) & (Qbus.<=0) & !((Pbus.==0) & (Qbus.==0)) # finds buses with at least one non-zero P/Q
    LoadBusID = find(LoadBusID)
    power_factor = zeros(nbus)
    power_factor[LoadBusID] = -Pbus[LoadBusID]./sqrt(Pbus[LoadBusID].^2 + Qbus[LoadBusID].^2)
    beta_S = zeros(nbus)
    beta_S[LoadBusID] = beta_P[LoadBusID].*power_factor[LoadBusID] + beta_Q[LoadBusID].*sqrt(1-power_factor[LoadBusID].^2)
    
    # beta_C: the sensitivity factors when switching a capacitor at buses assuming there was no capacitor before
    beta_C = zeros(nbus)
    Vmag[pq] = xm[ix_Vmag]
    beta_C[pq] = Vmag[pq].^2
    beta_C = beta_C.*beta_Q

    return beta, beta_P, beta_Q, beta_S, beta_C
end

end

