function FindSubGraphs(nodes,links)
	# find connected components in a graph
	nodes = vec(nodes) # make sure nodes are in vector format and not matrix
	N = length(nodes)
    nbr = size(links,1)
    bus_i = zeros(Int64,maximum(nodes))
    bus_i[nodes] = 1:N
	F = bus_i[links[:,1]]
	T = bus_i[links[:,2]]
    A = sparse([F,T],[T,F],ones(Int64,2*size(links,1)),N,N)

	grNo = 1
	graphNos = zeros(Int64,N)
    linkNos = zeros(Int64,nbr)
	next = 1
	while next != 0
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
        next = findfirst(graphNos.==0,1)
	end
    # find where each link is
    for i = 1:nbr
        this_busi = F[i]
        linkNos[i] = graphNos[this_busi]
    end

    return graphNos, linkNos
end

