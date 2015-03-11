function getYbus(ps::Caseps,includeShunts = true)
	# Extract the Ybus matrix from power system data
	# Usage: Ybus, Yf, Yt = getYbus(ps,includeShunts,useMATPOWER)
	# Inputs:
	#  ps is the ps structure
	#  inclueShunts is a Bool indicating whether to put the impdeance shunts
	#   into the Ybus. (Default = true)
	# Outputs:
	#  Ybus - is the matrix
	#  Yf is a matrix that calculates complex currents on the from end 
	#    of each branch: If = Yf*V
	#  Yt is a matrix that calculates complex currents on the to end
	#    of each branch: It = Yt*V

	# collect some data
	n = size(ps.bus,1)
	nbr = size(ps.branch,1)
	eps = 1e-9

	# extract data from the branch matrix
	F      = ps.bus_i[ps.branch[:,:from]]
	T      = ps.bus_i[ps.branch[:,:to]]
	status = ps.branch[:,:status]
	R      = ps.branch[:,:R]
	X      = ps.branch[:,:X]
	B      = ps.branch[:,:B]
	shift  = ps.branch[:,:shift]
	tap    = ps.branch[:,:tap]
	tap[abs(tap).<eps] = 1.0

	# calculate the branch impedance values
	y_series = 1./(R+im*X)
	tap_shift = tap .* exp(-im*pi/180 * shift)
	y_tt = status.*( y_series + im*B/2)
	y_ff = status.*( y_tt ./ tap.^2)
	y_ft = status.*(-y_series ./ conj(tap_shift))
	y_tf = status.*(-y_series ./ tap_shift)
	# build these values into a sparse ybus matrix
	Ybus = sparse(F,F,y_ff,n,n) + 
	       sparse(T,T,y_tt,n,n) + 
	       sparse(F,T,y_ft,n,n) +
	       sparse(T,F,y_tf,n,n)
	# add the shunts if requested
	if includeShunts && isdefined(symbol("ps.shunt"))
	    factor = ps.shunt[:,:fracZ].*ps.shunt[:,:factor]
	    y_shunt = factor.*(ps.shunt[:,:P] + im*ps.shunt[:,:Q])/ps.baseMVA
	    sh_bus_i = ps.bus_i[ps.shunt[:,:bus]]
	    y_sh_bus = sparsevec(sh_bus_i,y_shunt,n)
	    Ybus = Ybus + sparse([1:n],[1:n],y_sh_bus,n,n)
	end
	# calculate Yf and Yt
	i = [1:nbr, 1:nbr]     # double set of row indices    
	Yf = sparse(i, [F, T], [y_ff, y_ft], nbr, n)
	Yt = sparse(i, [F, T], [y_tf, y_tt], nbr, n)

	return Ybus, Yf, Yt

end