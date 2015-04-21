module Makeps

using DataFrames
import Base.deepcopy
export 
	Caseps, 
	addcolnames, 
	updateps, 
	Constps

type Caseps
    baseMVA::Float64
    bus::DataFrame
    branch::DataFrame 
    gen::DataFrame 
    shunt::DataFrame 
    areas::DataFrame
    gencost::DataFrame
    bus_i::Array{Int64,1}
    Caseps() = new()
end

function addcolnames(ps::Caseps)
	# add column names for ps.bus
	BusColnames = ["ID","type","Pd","Qd","Gs","Bs","area","Vmag","Vang","basekV",
		"zone","Vmax","Vmin","lamP","lamQ","muVx","muVn","locX","locY","comm_status",
		"status","power_from_sh","grid_comm"];
	nColBus = size(ps.bus,2)
	names!(ps.bus,[convert(Symbol,i) for i in BusColnames[1:nColBus]])

	# add column names for ps.branch
	BrColnames = ["from","to","R","X","B","rateA","rateB","rateC","tap","shift",
		"status","Pf","Qf","Pt","Qt","muF","muT","ImagF","ImagT","switchable",
		"fail_rate","type","contg_st"];
	nColBr = size(ps.branch,2)
	names!(ps.branch,[convert(Symbol,i) for i in BrColnames[1:nColBr]])

	# add column names for ps.gen
	GenColnames = ["bus","P","Q","Qmax","Qmin","Vsp","mBase","status","Pmax","Pmin",
        "muPx","muPn","muQx","muQn","type","cost","part_fact","RRU","RRD","NA","NA2"];
	nColGen = size(ps.gen,2)
	names!(ps.gen,[convert(Symbol,i) for i in GenColnames[1:nColGen]])

	# add column names for ps.shunt
	ShuntColnames = ["bus","P","Q","fracS","fracZ","factor","type","value","fracE",
		"gamma","nearGen","ID"];
	nColShunt = size(ps.shunt,2)
	names!(ps.shunt,[convert(Symbol,i) for i in ShuntColnames[1:nColShunt]])

	return ps
end

function updateps(ps::Caseps)
	# add column labels to ps
	ps = addcolnames(ps)

	# make bus ID numbers an integer
	ps.bus[:ID] = int(ps.bus[:ID])
	ps.branch[:from] = int(ps.branch[:from])
	ps.branch[:to] = int(ps.branch[:to])
	ps.gen[:bus] = int(ps.gen[:bus])
	ps.shunt[:bus] = int(ps.shunt[:bus])
	# make branch status a bool variable
	ps.branch[:status] = bool(ps.branch[:status])

	# add ps.bus_i to ps
	n = size(ps.bus,1)
	max_bus_no = maximum(ps.bus[:,:ID])
	ps.bus_i = zeros(Int64,max_bus_no)
	ps.bus_i[ps.bus[:,:ID]] = 1:n

	return ps
end


# define a type and function for constants. 
immutable Constps
	PQ::Int64
	PV::Int64
	REF::Int64
	Constps() = new(
		1, #PQ
		2, #PV
		3  #REF
		)
end

# add a method to copy ps cases

function Base.deepcopy(ps::Caseps) 
	newps = Caseps()
	for i = 1:length(names(ps))
		newps.(names(ps)[i]) = deepcopy(ps.(names(ps)[i]))
	end
	return newps
end

#=
# define a type and function for options. 
type Optps
	sim_StopThreshold::Float64
	Optps() = new()
end
function psoptions()
	opt = Optps()
	opt.sim_StopThreshold = 0

	return opt
end
=#

end
