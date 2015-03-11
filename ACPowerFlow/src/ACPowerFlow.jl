module ACPowerFlow
# This module contains all the types and functions for power flow
using Makeps
import Base.real, Base.imag # extend these two to compute real and imag of a sparse matrix
export 
	getYbus, 
	getSbus,
	mismatch_polar,
	mismatch_rectangular,
	nrsolve
	
# Extend real() and imag() for sparce matrices (they currently don't take a sparce matrix as input)
real(A::SparseMatrixCSC) = SparseMatrixCSC(A.m,A.n,A.colptr,A.rowval,real(A.nzval))
imag(A::SparseMatrixCSC) = SparseMatrixCSC(A.m,A.n,A.colptr,A.rowval,imag(A.nzval))

# Make Ybus
include("getYbus.jl")

# Make Sbus
include("getSbus.jl")

# Mismatch function for polar voltage
include("mismatch_polar.jl")

# Mismatch function for rectangular voltage
include("mismatch_rectangular.jl")

# Newton-Raphson solver
include("nrsolve.jl")

end


