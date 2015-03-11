function mismatch_polar(x::Vector{Float64},
						Ybus::SparseMatrixCSC{Complex{Float64},Int64},
						Vmag::Vector{Float64},
						Sg_bus::Vector{Complex{Float64}},
						Sd_bus::Vector{Complex{Float64}},
						pq::BitArray{1},
						pv::BitArray{1},
						ref::BitArray{1},
						PartFact::Vector{Float64};
						NeedMismatch::Bool = true,
						NeedJac::Bool = false)
	if isempty(PartFact)
		PartFactFlag = false
	else
		PartFactFlag = true
	end
	nBus = size(Ybus,1)
	if PartFactFlag && all(PartFact .== 0)
		error("No ramping generators provided.")
	end
	if sum(ref) != 1 
		error("Must have only one ref bus"); 
	end
	# build an index
	npq = sum(pq)
	ix_theta = 1:(nBus-1)
	ix_Vmag  = (1:npq) + maximum(ix_theta)
	if PartFactFlag
	    # add in a variable for the generator ramping
	    ix_rho = 1 + maximum(ix_Vmag)
	    nx = nBus + npq
	else
	    nx = nBus - 1 + npq
	end
	# extract things from x
	theta = zeros(nBus)
	theta[!ref] = x[ix_theta]

	Vmag[pq]    = x[ix_Vmag]
	if PartFactFlag
	    rho = x[ix_rho]
	end
	# calculate the voltage
	V = Vmag.*exp(im*theta)
	# calculate the total load according to ZIPE matrix
	if size(Sd_bus,2) == 1
	    SZipe = Sd_bus
	elseif size(Sd_bus,2) == 5
	    S_Z = Sd_bus[:,1] .* Vmag.^2
	    S_I = Sd_bus[:,2] .* Vmag
	    S_P = Sd_bus[:,3]
	    S_E = Sd_bus[:,4] .* Vmag.^Sd_bus[:,5]
	    SZipe = S_Z + S_I + S_P + S_E
	else
	    error("zipe load model matrix is not the right size.")
	end
	if NeedMismatch
		Sbus = Sg_bus - SZipe
		
		if PartFactFlag
		    # make adjustments to Sg for gen ramping
		    RampGen        = ref | pv
		    Sbus[RampGen]  = Sbus[RampGen] + PartFact[RampGen]*rho
		end
		# compute the final mismatch
		miscx = (V .* conj(Ybus * V)) - Sbus
		if PartFactFlag
		    g = [real(miscx); imag(miscx[pq])]
		else
		    g = [real(miscx[!ref]); imag(miscx[pq])]
		end
	else
		g = Float64[]
	end
	if NeedJac
		# make the Jacobian	    
	    # do some matrix algebra (borrowed from MATPOWER)
	    Ibus = Ybus * V
	    diagV     = sparse(1:nBus, 1:nBus, V, nBus, nBus)
	    diagIbus  = sparse(1:nBus, 1:nBus, Ibus, nBus, nBus)
	    diagVnorm = sparse(1:nBus, 1:nBus, exp(im*theta), nBus, nBus)
	    dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
	    dSbus_dVa = im * diagV * conj(diagIbus - Ybus * diagV)
	    
	    if PartFactFlag
	        # now put these into dg_dx
	        # dP_dtheta
	        rows,cols,values = findnz(real(dSbus_dVa[:,!ref]))
	        dg_dx = sparse(rows,cols,values,nx,nx)
	        # dQ_dtheta
	        rows,cols,values = findnz(imag(dSbus_dVa[pq,!ref]))
	        dg_dx = dg_dx + sparse(rows+nBus,cols,values,nx,nx)
	        # dP_dVmag
	        rows,cols,values = findnz(real(dSbus_dVm[:,pq]))
	        dg_dx = dg_dx + sparse(rows,cols+(nBus-1),values,nx,nx)
	        # dQ_dVmag
	        rows,cols,values = findnz(imag(dSbus_dVm[pq,pq]))
	        dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),values,nx,nx)
	        # dP_drho
	        dg_dx = dg_dx + sparse((1:nBus)[RampGen],ix.rho,-PartFact[RampGen],nx,nx)
	    else
	        # no participation factor (standard power flow)
	        # dP_dtheta
	        rows,cols,values = findnz(real(dSbus_dVa[!ref,!ref]))
	        dg_dx = sparse(rows,cols,values,nx,nx)
	        # dQ_dtheta
	        rows,cols,values = findnz(imag(dSbus_dVa[pq,!ref]))
	        dg_dx = dg_dx + sparse(rows+(nBus-1),cols,values,nx,nx)
	        # dP_dVmag
	        rows,cols,values = findnz(real(dSbus_dVm[!ref,pq]))
	        dg_dx = dg_dx + sparse(rows,cols+(nBus-1),values,nx,nx)
	        # dQ_dVmag
	        rows,cols,values = findnz(imag(dSbus_dVm[pq,pq]))
	        dg_dx = dg_dx + sparse(rows+(nBus-1),cols+(nBus-1),values,nx,nx)
	    end
	    
	    if size(Sd_bus,2) == 5 # zipe load model 
	        # fix the derivatives with ZIP[E] contributions
	        dP_E_dVmag = Sd_bus[:,5].*real(Sd_bus[:,4]).*Vmag.^(Sd_bus[:,5]-1)
	        dQ_E_dVmag = Sd_bus[:,5].*imag(Sd_bus[:,4]).*Vmag.^(Sd_bus[:,5]-1)
	        # fix the derivatives with [Z]IPE contributions
	        dP_Z_dVmag = 2.0*real(Sd_bus[:,1]).*Vmag
	        dQ_Z_dVmag = 2.0*imag(Sd_bus[:,1]).*Vmag
	        # fix the derivatives with Z[I]PE contributions
	        dP_I_dVmag = real(Sd_bus[:,2])
	        dQ_I_dVmag = imag(Sd_bus[:,2])
	        
	        
	        SS_zipe = S_zipe
	        SS_zipe[ref | pv] = 0     # assume that exponential loads are not located in ref/pv buses
	        
	        rows = find(SS_zipe); cols = find(SS_zipe[pq])
	        dg_dx = dg_dx + sparse(rows,cols+(nBus-1),dP_E_dVmag[rows],nx,nx)
	        dg_dx = dg_dx + sparse(rows,cols+(nBus-1),dP_Z_dVmag[rows],nx,nx)
	        dg_dx = dg_dx + sparse(rows,cols+(nBus-1),dP_I_dVmag[rows],nx,nx)
	        
	        rows = find(SS_zipe[pq]); cols = rows
	        dQ_E_dVmag = dQ_E_dVmag[pq]
	        dQ_Z_dVmag = dQ_Z_dVmag[pq]
	        dQ_I_dVmag = dQ_I_dVmag[pq]
	        dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),dQ_E_dVmag[rows],nx,nx)
	        dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),dQ_Z_dVmag[rows],nx,nx)
	        dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),dQ_I_dVmag[rows],nx,nx)
	    end
	else
		dg_dx = spzeros(1,1)
	end

	return g, dg_dx
end