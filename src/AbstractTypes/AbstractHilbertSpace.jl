abstract type AbstractHilbertSpace end
abstract type AbstractConstraint end

constraints(::AbstractHilbertSpace) = error("No constraints defined for this Hilbert space.")
size(::AbstractHilbertSpace) = error("No size defined for this Hilbert space.")

