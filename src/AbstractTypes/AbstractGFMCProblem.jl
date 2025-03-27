"""
    Abstract type for GFMC problems.
# Interface:
- `runGFMC!(problem::AbstractGFMCProblem,args...;kwargs...)`: Run the GFMC algorithm for the given problem.
"""
abstract type AbstractGFMCProblem end
