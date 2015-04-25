function acpf!(ps::Caseps;
               verbose::Bool = true,
               PartFactFlag::Bool = false,
               PartFact::Vector{Float64} = Float64[],
               PolarFlag::Bool = true)
    # if PartFact is empty (the default) and PartFactFlag is true, all the PV buses (and the reference bus)
    # will be added to the list of participation factor list
    nBus = size(ps.bus,1)
    nodes = ps.bus[:,:ID]
    br_status = ps.branch[:,:status]
    links = ps.branch[br_status, [:from, :to]]
    graphNos, linkNos = FindSubGraphs(nodes,links)
    Vmag = ps.bus[:Vmag] # note that Vmag and theta are passed by reference, so they are in fact the same as ps.bus[:Vmag] and ps.bus[:Vang]
    ps.bus[:,:Vang] = 0.0
    theta = ps.bus[:Vang]
    n_subgraph = maximum(graphNos)
    Ybus, Yf, Yt = getYbus(ps)
    Sbus, Sd, Sg = getSbus(ps)
    Sd = sum(Sd,2)
    D = ps.bus_i[ps.shunt[:,:bus]]
    Sd_bus = zeros(Complex{Float64},nBus); Sd_bus[D] = Sd
    F = ps.bus_i[ps.branch[:,:from]]
    T = ps.bus_i[ps.branch[:,:to]]
    not_converged = Int64[] # list of the islands that did not converge
    for i = 1:n_subgraph
        # find buses in this island
        busi_sub = find(graphNos .== i)
        # find branches in this island
        bri_sub = find(linkNos .== i)
        mismatch, x0, ix_Vmag_sub, ix_theta_sub, pq_sub, ref_sub = 
            init_mismatch(ps, busi_sub, verbose, PartFactFlag, PartFact, PolarFlag)
        ## Solve power flow
        # newton-raphson
        x, converged = nrsolve(mismatch, x0, verbose=verbose, linesearch="exact")
        
        if converged
            # save the results
            if PolarFlag
                Vmag_sub = Vmag[busi_sub]
                Vmag_sub[pq_sub] = x[ix_Vmag_sub]
                Vmag[busi_sub] = Vmag_sub # this updates ps.bus[:,:Vmag]
                theta_sub = theta[busi_sub]
                theta_sub[!ref_sub] =  x[ix_theta_sub]
                theta_sub[ref_sub]  = 0.
                theta[busi_sub] = theta_sub # this updates ps.bus[:,:Vang]
            else
                # ignore this for now
                #=
                e = Vmag
                f = zeros(n)
                e[!ref] = x[ix_e]
                f[!ref] = x[ix_f]
                Vmag = sqrt(e.^2 + f.^2)
                theta = atan(f./e)
                =#
            end
            # update results in ps only for this island
            V_sub = Vmag_sub.*exp(im*theta_sub)
            V = zeros(Complex{Float64},nBus)
            V[busi_sub] = V_sub
            V_sub = Vmag_sub.*exp(im*theta_sub)
            Yf_sub = Yf[bri_sub,busi_sub]
            Yt_sub = Yt[bri_sub,busi_sub]
            If_sub = Yf_sub*V_sub
            It_sub = Yt_sub*V_sub
            F_sub = F[bri_sub]
            Sf_sub = V[F_sub] .* conj(If_sub)
            T_sub = T[bri_sub]
            St_sub = V[T_sub] .* conj(It_sub)
            ps.branch[bri_sub,:ImagF] = abs(If_sub)
            ps.branch[bri_sub,:ImagT] = abs(It_sub)
            ps.branch[bri_sub,:Pf] = real(Sf_sub) * ps.baseMVA
            ps.branch[bri_sub,:Qf] = imag(Sf_sub) * ps.baseMVA
            ps.branch[bri_sub,:Pt] = real(St_sub) * ps.baseMVA
            ps.branch[bri_sub,:Qt] = imag(St_sub) * ps.baseMVA
            # update ps.gen
            Ybus_sub = Ybus[busi_sub,busi_sub]
            Sg_bus_sub = V_sub.*conj(Ybus_sub*V_sub) + Sd_bus[busi_sub]
            Pg_bus_sub = real(Sg_bus_sub)
            Qg_bus_sub = imag(Sg_bus_sub)
            gen_busID = ps.gen[:,:bus]
            gen_busi = ps.bus_i[gen_busID]
            for (j,this_bus_i) = enumerate(busi_sub)
                gi = (gen_busi .== this_bus_i)
                if V[this_bus_i] == 0.
                    ps.gen[gi,:P] = 0.
                    ps.gen[gi,:Q] = 0.
                else
                    if sum(gi) == 1
                        Pg_ratio = 1
                    else
                        Pg_ratio = ps.gen[gi,:Pmax]/sum(ps.gen[gi,:Pmax])
                    end
                    ps.gen[gi,:P] = Pg_bus_sub[j] * Pg_ratio * ps.baseMVA
                    ps.gen[gi,:Q] = Qg_bus_sub[j] * Pg_ratio * ps.baseMVA
                end
            end
        else
            if verbose
                println("Power flow did not converge on island $i of $n_subgraph.")
            end
            not_converged = [not_converged, i]
            # update ps with 0 values for non-converged islands
            Vmag[busi_sub] = 0.
            theta[busi_sub] = 0.
            ps.branch[bri_sub,:ImagF] = 0.
            ps.branch[bri_sub,:ImagT] = 0.
            ps.branch[bri_sub,:Pf] = 0.
            ps.branch[bri_sub,:Qf] = 0.
            ps.branch[bri_sub,:Pt] = 0.
            ps.branch[bri_sub,:Qt] = 0.
            gen_busID = ps.gen[:,:bus]
            gen_busi = ps.bus_i[gen_busID]
            for j = busi_sub
                gi = (gen_busi .== j)
                ps.gen[gi,:P] = 0.
                ps.gen[gi,:Q] = 0.
            end
        end
    end

    return not_converged

    #=
    # solve the unconstrained optimization problem 1/2 * g(x)' * g(x) using levenberg-marquardt 
    # algorithm, where g(x) is the mismatch in power injections
    g(x) = mismatch(x)[1]
    Jac(x) = mismatch(x, NeedJac=true)[2]
    results = Optim.levenberg_marquardt(g, Jac, x0)
    F(x) = g(x)'*g(x)
    =#
end

function init_mismatch(ps::Caseps,
                       busi_sub::Vector{Int64},
                       verbose::Bool,
                       PartFactFlag::Bool,
                       PartFact::Vector{Float64},
                       PolarFlag::Bool)
    # initializes mismatch with all the necessary inputs and only accepting x0 as input
    C = Constps()
    # extract stuff from ps
    nBus = size(ps.bus,1)
    if isempty(busi_sub)
        busi_sub = 1:nBus
    end
    nBus_sub = length(busi_sub)
    G = ps.bus_i[ps.gen[:,:bus]]
    D = ps.bus_i[ps.shunt[:,:bus]]
    Ybus, Yf, Yt = getYbus(ps)
    Vmag = ps.bus[:,:Vmag]
    theta0_sub = zeros(nBus_sub)
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
    if PartFactFlag && isempty(PartFact)
        PartFact = pv | ref
        PartFact = convert(Vector{Float64}, PartFact)
    end
    if sum(ref) == 0 || sum(ref) > 1
        error("No reference bus specified.")
    end
    ref_sub = ref[busi_sub]
    if sum(ref_sub) == 0
        # assign the bus with the biggest generation as the reference
        Pg_bus = real(Sg_bus)
        ref_sub_busi = find(Pg_bus[busi_sub] .== maximum(Pg_bus[busi_sub]))[1]
        ref_sub[ref_sub_busi] = true
    end
    Vmag_sub = Vmag[busi_sub]
    pq_sub = pq[busi_sub]
    npq_sub = sum(pq_sub)
    if PolarFlag
        # set up x for polar power flow
        ix_theta_sub = 1:(nBus_sub-1)
        ix_Vmag_sub  = (1:npq_sub) + maximum(ix_theta_sub)
        if PartFactFlag
            ix_rho_sub = 1 + maximum(ix_Vmag_sub) # generator ramping variable
            nx = maximum(ix_rho_sub)
        else
            nx = maximum(ix_Vmag_sub)
        end
        x0 = zeros(nx)
        x0[ix_theta_sub] = theta0_sub[!ref_sub]
        x0[ix_Vmag_sub]  = Vmag_sub[pq_sub]
        if PartFactFlag
            x0[ix_rho_sub] = 0
        end
        mismatch_func = mismatch_polar
    else
        # set up x for rectangular power flow
    end
    if isempty(busi_sub)
        Ybus_sub, Vmag_sub, Sg_bus_sub, Sd_bus_sub, pq_sub, pv_sub, ref_sub, PartFact_sub = 
        Ybus, Vmag, Sg_bus, Sd_bus, pq, pv, ref, PartFact
    else
        Ybus_sub = Ybus[busi_sub,busi_sub]
        Yf_sub = Yf[:,busi_sub] # because status is taken into account in making Ybus, Yf and Yt, this should work
        Yt_sub = Yt[:,busi_sub]
        Vmag_sub = Vmag[busi_sub]
        Sg_bus_sub   = Sg_bus[busi_sub]
        Sd_bus_sub   = Sd_bus[busi_sub]
        pv_sub   = pv[busi_sub]
        if PartFactFlag
            PartFact_sub = PartFact[busi_sub]
        else
            PartFact_sub = Float64[]
        end
    end

    mismatch(x;args...) = mismatch_func(x,Ybus_sub,Vmag_sub,Sg_bus_sub,Sd_bus_sub,pq_sub,pv_sub,ref_sub,PartFact_sub;args...)

    return mismatch, x0, ix_Vmag_sub, ix_theta_sub, pq_sub, ref_sub

end

