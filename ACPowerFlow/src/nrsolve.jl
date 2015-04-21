import Optim
function nrsolve(mismatch::Function,
				 x0::Vector{Float64};
				 max_iter::Integer = 20,
				 gtol::Real = 1e-9,
				 verbose::Bool = true,
				 linesearch::String = "none")

	# solve power flow with the mismatch function
	converged = false
	iter = 0 
	x = deepcopy(x0)
	if verbose
		println("Iter   Max(|g|)   |g|_2      max(Jac'*g)      alpha")
	end
	norm2(x) = sum(x.^2)

	for iter = 1:max_iter
		# evaluate the function
		g,Jac = mismatch(x,NeedJac=true)
		max_mismatch = maximum(abs(g))
		if max_mismatch < gtol
			converged = true
            if verbose 
    			println("Solution found.")
            end
			break
		end

		# check first order optimality condition
        Jg = Jac.'*g;
        max_Jg = maximum(Jg)
        g2 = sqrt(sum(g.^2))

        # choose the search direction
        p = -(Jac\g)
        
        # save x for the plot below with different alpha's
        x_previous = deepcopy(x)
        
        # do some sort of line search to select the step size and the new x
        if linesearch == "none"
        	alpha = 1
        	x = x_previous + alpha*p
        elseif linesearch == "exact"
    		obj_fun(a) = norm2(mismatch(x_previous + a*p)[1])
        	results = Optim.optimize(obj_fun, -5.0, 5.0)
        	if results.converged
        		alpha = results.minimum
        	else
        		alpha = 0
        	end
        	x = x_previous + alpha*p
        end
        # print something
        if verbose
            @printf("%4d %10.7f %10.7f %10.7f %10.7f \n", iter, max_mismatch, g2, max_Jg, alpha)
        end
        if alpha < 1e-6
        	break
        end
    end

    if  verbose && !converged
    	println(" Did not find a solution to g(x) = 0.")
	end

	return x, converged
end
