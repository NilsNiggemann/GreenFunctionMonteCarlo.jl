abstract type AbstractPropagator end

propagateWalkers!(X::AbstractWalkerEnsemble, dX::AbstractArray, P::AbstractPropagator,params) = throw(MethodError(propagateWalkers!, (X, dX, P, params)))