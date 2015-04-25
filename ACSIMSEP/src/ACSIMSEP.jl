module ACSIMSEP

using Makeps, ACPowerFlow

export get_left_eigv,
       find_sens_fact,
       find_beta

# find the sensitivity factors to laod shedding
include("find_sens_fact.jl")
# include("find_sens_fact_old.jl")

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

function find_beta(mismatch::Function, x0::Vector{Float64}; verbose::Bool = true)
    # This function finds the minimum distance in parameter space (beta) from 
    # an unsolvable point to the solvable region (based on Overbye's paper)
    x, converged = nrsolve(mismatch, x0, verbose=verbose, linesearch="exact")
    if converged == true
        dist = 0
        println("Warning: The initial power flow converged.")
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


end

