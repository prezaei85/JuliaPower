function  getSbus(ps::Caseps,only_PQ_loads::Bool = false)
	# extract the Sbus injection from a power system structure
	# usage: Sbus, Sbus_load, Sbus_gen = getSbus(ps,only_PQ_loads)
	# Sbus will be a matrix with the components from the ZIPE model, as defined in shunt

	n = size(ps.bus,1)

	# generators
	status = ps.gen[:,:status]
	Sg = status.*(ps.gen[:,:P] + im*ps.gen[:,:Q])/ps.baseMVA
	ge_bus_i = ps.bus_i[ps.gen[:,:bus]]
	Sg_bus = sparsevec(ge_bus_i,Sg,n)

	# loads
	sh_bus_i = ps.bus_i[ps.shunt[:,:bus]]
	# figure out the ZIPE load types
	f_Z = ps.shunt[:,:fracZ]
	f_P = ps.shunt[:,:fracS]
	f_E = ps.shunt[:,:fracE]
	f_I = 1 - f_Z - f_P - f_E
	# get the actual load and locations
	Sd = ps.shunt[:,:factor].*(ps.shunt[:,:P] + im*ps.shunt[:,:Q])/ps.baseMVA
	# check the ZIPE load model
	if any(f_I.<0) || any(f_I.>1.0)
	    error("The load model in ps.shunt is not valid")
	end
	# get the exponent for the E portion
	gamma = ps.shunt[:,:gamma]
	# build the load matrix
	Sd_model = [Sd.*f_Z Sd.*f_I Sd.*f_P Sd.*f_E gamma]

	#build the final Sd_bus
	if only_PQ_loads
	    Sd_bus = sparsevec(sh_bus_i,Sd.*f_P,n)
	else
	    Sd_bus = sparsevec(sh_bus_i,Sd,n)
	end

	#combine to get Sbus
	Sbus = Sg_bus - Sd_bus

	return Sbus, Sd_model, Sg

end