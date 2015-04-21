function FindSubsGraphs(nodes,links)
	# find connected components in a graph
	nodes = vec(nodes) # make sure nodes are in vector format and not matrix
	N = length(nodes)
    bus_i = zeros(Int64,maximum(nodes))
    bus_i[nodes] = 1:N
	F = bus_i[links[:,1]]
	T = bus_i[links[:,2]]
    A = sparse([F,T],[T,F],ones(Int64,2*size(links,1)),N,N)

	grNo = 1
	graphNos = zeros(Int64,N)
	next = 1
	while !isempty(next)
        included = falses(N)
        included[next] = true
        oldLen = 0
        while sum(included) != oldLen
            oldLen = sum(included)
            Ai = findnz(A[:,included])[1]
            included[Ai] = true
        end
        graphNos[included] = grNo
        grNo = grNo+1
        next = find(graphNos.==0)
	end
    nSubGraphs = grNo-1
    return graphNos
end