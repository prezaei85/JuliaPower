function find_sens_fact(ps::Caseps, PartFactFlag, not_converged; verbose = true)
    nBus = size(ps.bus,1)
    nodes = ps.bus[:,:ID]
    br_status = ps.branch[:,:status]
    links = ps.branch[br_status, [:from, :to]]
    graphNos, linkNos = FindSubGraphs(nodes,links)
    n_subgraph = maximum(graphNos)
    beta_P = zeros(nBus); beta_Q = zeros(nBus)
    beta_S = zeros(nBus); beta_C = zeros(nBus)
    Sbus, = getSbus(ps)
    Vmag = ps.bus[:,:Vmag]
    C = Constps()
    for i = not_converged
        # find buses in this island
        busi_sub = find(graphNos .== i)
        nBus_sub = length(busi_sub)
        beta_P_sub = zeros(nBus_sub); beta_Q_sub = zeros(nBus_sub)
        beta_S_sub = zeros(nBus_sub); beta_C_sub = zeros(nBus_sub)
        # find branches in this island
        bri_sub = find(linkNos .== i)
        mismatch, x0, ix_Vmag_sub, ix_theta_sub, pq_sub, ref_sub = 
            init_mismatch(ps, busi_sub, verbose, PartFactFlag, Float64[], true)
        # find the minimum distance to solvability
        xm, Jm, beta = find_beta(mismatch,x0,verbose = verbose);

        # get the left eigen vector of Jacobian associated with zero eigen value 
        w_m = get_left_eigv(Jm)
        # resolve the direction of w_m (it needs to be outward the solvable region)
        g(x) = mismatch(x)[1]
        if dot(g(xm),w_m) < 0 # resolve the direction of the eigen vector
            w_m = -w_m
        end
        
        if PartFactFlag
            # beta_P: the sensitivity factors when reducing real power at buses
            beta_P_sub = w_m[1:nBus_sub]
            
            # beta_Q: the sensitivity factors when reducing reactive power at buses
            beta_Q_sub[pq_sub] = w_m[nBus_sub+1:end]
        else
            # beta_P: the sensitivity factors when reducing real power at buses
            beta_P_sub[!ref_sub] = w_m[1:nBus_sub-1]
            
            # beta_Q: the sensitivity factors when reducing reactive power at buses
            beta_Q_sub[pq_sub] = w_m[nBus_sub:end]
        end
        # beta_S: the sensitivity factors when reducing complex power at buses while maintaining power factor
        # only considers buses with positive load
        Sbus_sub = Sbus[busi_sub]
        Pbus_sub = real(Sbus_sub)
        Qbus_sub = imag(Sbus_sub)
        load_busi_sub = (Pbus_sub.<=0) & (Qbus_sub.<=0) & !((Pbus_sub.==0) & (Qbus_sub.==0)) # finds buses with at least one non-zero P/Q
        load_busi_sub = find(load_busi_sub)
        power_factor = zeros(nBus_sub)
        power_factor[load_busi_sub] = -Pbus_sub[load_busi_sub]./sqrt(Pbus_sub[load_busi_sub].^2 + Qbus_sub[load_busi_sub].^2)
        beta_S_sub[load_busi_sub] = beta_P_sub[load_busi_sub].*power_factor[load_busi_sub] + 
                                    beta_Q_sub[load_busi_sub].*sqrt(1-power_factor[load_busi_sub].^2)
        
        # beta_C: the sensitivity factors when switching a capacitor at buses assuming there was no capacitor before
        Vmag_sub = Vmag[busi_sub]
        Vmag_sub[pq_sub] = xm[ix_Vmag_sub]
        beta_C_sub[pq_sub] = Vmag_sub[pq_sub].^2
        beta_C_sub = beta_C_sub.*beta_Q_sub

        # update the original list
        beta_P[busi_sub] = beta_P_sub
        beta_Q[busi_sub] = beta_Q_sub
        beta_S[busi_sub] = beta_S_sub
        beta_C[busi_sub] = beta_C_sub
    end

    return beta, beta_P, beta_Q, beta_S, beta_C
end


