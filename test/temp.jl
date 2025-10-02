using GreenFunctionMonteCarlo

σz(n::Bool) = (1 - 2 * n)
σz(i, conf::AbstractArray) = σz(conf[i])

# Define parameters
J = 1.0       # Interaction strength
h = 1.5       # Transverse field strength
lattice_size = 6  # Number of spins
periodic = false  # No periodic boundary conditions

using SparseArrays, LinearAlgebra
⊗(x,y) = kron(x,y)
function TransverseFieldIsing_sparse(;N,h)
    id = [1 0; 0 1] |> sparse
    σˣ = [0 1; 1 0] |> sparse
    σᶻ = [1 0; 0 -1] |> sparse
    
    first_term_ops = fill(id, N)
    first_term_ops[1] = σᶻ
    first_term_ops[2] = σᶻ
    
    second_term_ops = fill(id, N)
    second_term_ops[1] = σˣ
    
    H = spzeros(Int, 2^N, 2^N) # note the spzeros instead of zeros here
    for i in 1:N-1
        H -= foldl(⊗, first_term_ops)
        first_term_ops = circshift(first_term_ops,1)
    end
    
    for i in 1:N
        H -= h*foldl(⊗, second_term_ops)
        second_term_ops = circshift(second_term_ops,1)
    end
    H
end
bit_rep(num::Integer, N::Integer) = BitArray(parse(Bool, i) for i in string(num, base=2, pad=N))
H_TFI = Array(TransverseFieldIsing_sparse(N=lattice_size, h=h))

function corrs(state)
    N = Int(log2(length(state)))
    SiSj = zeros(N,N)
    for i in eachindex(state)
        bstate = bit_rep(i-1,N)
        # println(bstate,state[i])
        for (i_site,spin1) in enumerate(bstate)
            for (j_site,spin2) in enumerate(bstate)
                s1 = σz(spin1)
                s2 = σz(spin2)
                SS = s1 * s2
                SiSj[i_site,j_site] += state[i]^2 * SS
            end
        end
    end
    return SiSj
end
Energy_exact = minimum(eigen(H_TFI).values)
corrs_exact = corrs(eigen(H_TFI).vectors[:,1])

##
function transverse_field_ising(Nsites, h, J; periodic = false)
    Hilbert = BosonHilbertSpace(Nsites, HardCoreConstraint())

    moves = [Bool[0 for _ in 1:Nsites] for _ in 1:Nsites]
    offdiagElements = zeros(Float64, Nsites)

    for i in eachindex(moves)
        moves[i][i] = true
        offdiagElements[i] = -h
    end

    function Hxx(conf)
        E = -J * sum(σz(i, conf) * σz(i + 1, conf) for i in eachindex(conf)[1:end-1])
        if periodic
            E += -J * σz(Nsites, conf) * σz(1, conf)
        end
        return E
    end

    H = localOperator(moves, offdiagElements, DiagOperator(Hxx), Hilbert)
    return (; Hilbert, H)
end

# Define the Hamiltonian
(;Hilbert,H) = transverse_field_ising(lattice_size, h, J; periodic)
##
NSteps_total = 2^17  # Number of Monte Carlo steps
num_bins = 20
NSteps = NSteps_total ÷ num_bins  # Number of steps per parallel run
NWalkers = 1  # Number of walkers
dtau = 0.1    # imaginary time propagation step
mProj = 100
startConfig = BosonConfig(Hilbert) # creates an initial configuration where all occupation numbers are 0
##

using Random, CairoMakie, Statistics
function E_critPoint_exact(L, h=1, periodic=false)
    (!periodic && h==1) || return NaN 
    
    return 1 - csc(pi / (2 * (2 * L + 1)))
end

vij_jastrow = zeros(Float32,lattice_size,lattice_size)
for i in axes(vij_jastrow, 1)[1:end-1]
    vij_jastrow[i,i+1]=vij_jastrow[i+1,i] = 0.4
end

logψ = Jastrow(
    zeros(Float32,lattice_size),
    vij_jastrow,
)
P = ProblemEnsemble([GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert) for i in 1:num_bins])
Observers_jastrow = [ConfigObserver(rand!(copy(startConfig)), NSteps, NWalkers) for _ in 1:num_bins]

runGFMC!(P, NoObserver(),1000,logger=nothing) #equilibrate
runGFMC!(P, Observers_jastrow,NSteps,logger=nothing)

energies_jastrow = [getEnergies(Observer, mProj) for Observer in Observers_jastrow]
tau = (0:mProj-1) * dtau
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"Energy$$")

    lines!(ax, tau, mean(energies_jastrow), label=L"shortrange Jastrow $N_w=%$NWalkers$")
    band!(ax, tau, mean(energies_jastrow) - std(energies_jastrow), mean(energies_jastrow) + std(energies_jastrow), color=(:green, 0.2))

    hlines!(ax, [Energy_exact], color=:black, label=L"Exact$$", linestyle=:dash)
    # ylims!(ax,-14.95/12*lattice_size, -14.68/12*lattice_size)
    axislegend(ax, position=:rt)
    fig
end
##
Observable = OccupationNumber(lattice_size)

n_avg = [stack(getObs_diagonal(O,Observable,1:mProj)) for O in Observers_jastrow]

n_avg_mean = mean(n_avg)
n_avg_std = std(n_avg)

let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"$\langle S^z_i\rangle$")
    tau =( 0:mProj-1) * dtau
    for i in (1:lattice_size)[1:2:end]
        lines!(ax, tau, n_avg_mean[i,:])
        band!(ax, tau, n_avg_mean[i,:] - n_avg_std[i,:], n_avg_mean[i,:] + n_avg_std[i,:], alpha = 0.5)
    end
    fig
end

##
struct SpinCorrelationsVec{T<:Real} <: AbstractObservable
    ObservableBuffer::Vector{T}
end
SpinCorrelationsVec(Nsites) = SpinCorrelationsVec(zeros(Nsites*Nsites))

Base.copy(O::SpinCorrelationsVec) = SpinCorrelationsVec(copy(O.ObservableBuffer))
GreenFunctionMonteCarlo.obs(O::SpinCorrelationsVec) = O.ObservableBuffer
function (O::SpinCorrelationsVec)(out,config)
    outmat = reshape(out, length(config), length(config))
    for i in eachindex(config)
        for j in eachindex(config)
            outmat[i,j] = σz(i,config) * σz(j,config)
        end
    end
    return out
end

Observable = SpinCorrelationsVec(lattice_size)
Corrs = [stack(getObs_diagonal(O,Observable,1:mProj)) for O in Observers_jastrow]

Corr_mean = reshape(mean(Corrs), lattice_size, lattice_size, mProj)
Corr_std = reshape(std(Corrs), lattice_size, lattice_size, mProj)

let 
    fig = Figure()
    i = 2
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"$\langle S^z_{%$i}S^z_j\rangle$")
    tau =( 0:mProj-1) * dtau
    
    for j in (1:lattice_size)[1:2:end]
        lines!(ax, tau, Corr_mean[i,j,:])
        band!(ax, tau, Corr_mean[i,j,:] - Corr_std[i,j,:], Corr_mean[i,j,:] + Corr_std[i,j,:], alpha = 0.5)
    end
    fig
end
##
NSteps_total = 2^18  # Number of Monte Carlo steps
NWalkers = 1
import GreenFunctionMonteCarlo as GFMC
problem = GFMCProblem(startConfig, NWalkers, ContinuousTimePropagator(dtau); logψ, H, Hilbert)

mean_TotalWeights, w_avg_estimate = GFMC.estimate_weights_continuousTime!(problem;Nepochs=20,Nsamples=800)

CT = ContinuousTimePropagator(dtau, w_avg_estimate)
# num_bins = NSteps_total ÷2^10
num_bins = 100
BAcc = BasicAccumulator(mProj, NWalkers;weight_normalization = mean_TotalWeights, num_bins = num_bins , bin_elements = NSteps_total ÷ num_bins)
conf_averages = WalkerAVGObserver(startConfig, NSteps_total)
ObsAcc = ObservableAccumulator(SpinCorrelationsVec(lattice_size),BAcc, mProj, NWalkers, 1)
Observer = CombinedObserver((BAcc, ObsAcc,conf_averages))

# Logging.with_logger(Logging.ConsoleLogger()) do
@time runGFMC!(problem, Observer, NSteps_total; Propagator = CT,parallelization = SingleThreaded())
# @time runGFMC!(problem, Observer, NSteps_total÷2+1:NSteps_total-1; Propagator = CT,parallelization = SingleThreaded())
# end
##
using MakieHelpers
bunchErrs_all = GFMC.bunching_errors_recursive.(eachrow(conf_averages.average_configs))

with_theme(theme_SimpleTicks()) do
    fig = Figure(size = (800, 500))
    bunchErrs = mean(bunchErrs_all)
    lens = 2 .^eachindex(bunchErrs)
    ax = Axis(fig[1, 1], xlabel=L"bunching interval$$", ylabel=L"apparent error$$";xscale = log2,xminorticksvisible = true,xticks = (lens[1:4:end],string.(lens[1:4:end])),
    xminorticks = IntervalsBetween(3))
    # scatter!(ax,lens,bunchErrs,markersize=10)
    errlines!(ax, lens,bunchErrs,std(bunchErrs_all),color=:blue,markersize=10)
    # err_est = bunchErrs[end-5]
    # hlines!(ax, [err_est],label = L"estimated error$$")
    # axislegend(ax, position=:rt)
    # errlines(2 .^(1:floor(Int,log2(N))),reverse(dropmean(bunchErrs,dims=2)),reverse(dropstd(bunchErrs,dims=2)),markersize=10,axis = (;xscale = log2))
    fig
end
##
Energies = GFMC.get_energy_from_accumulator_bunching(BAcc,30)
Energies_mean = mean(Energies)
Energies_std = std(Energies)

let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"Energy$$")
    tau =( 0:mProj-1) * dtau
    
    err0 = Energies_std[begin]
    err1 = Energies_std[end]
    ll = lines!(ax, tau, Energies_mean)
    band!(ax, tau, Energies_mean - Energies_std, Energies_mean + Energies_std, alpha = 0.5)
    hlines!(ax, Energy_exact/lattice_size, color=:black, label=L"Exact$$", linestyle=:dash)
    text!(ax, Point(0.3,0.7),text = L"\Delta E(0) = %$err0,",space = :relative, color=:black)
    text!(ax, Point(0.3,0.65),text = L"\Delta E(τ) = %$err1 ",space = :relative, color=:black)
    fig
end
##
All_en_err_bunch = stack([std(GFMC.get_energy_from_accumulator_bunching(BAcc,i)) for i in 1:50])
let 
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=L"$n_\textrm{bunch}$", ylabel=L"\Delta Energy$$",yscale = log10,)
    tau =( 0:mProj-1) * dtau
    
    tau_inds = axes(All_en_err_bunch, 1)
    for i_tau in (tau_inds[1], tau_inds[end])
        ll = scatterlines!(ax, All_en_err_bunch[i_tau,:], label=L"τ = %$tau[i_tau]$$",markersize=5)
    end

    err_mean = mean(All_en_err_bunch, dims=1)[:]
    err_std = std(All_en_err_bunch, dims=1)[:]

    lines!(ax, err_mean, color=:black, label=L"mean$$")
    band!(ax,eachindex(err_mean), err_mean - err_std, err_mean + err_std, alpha = 0.5, color=:black)
    fig
    
end
##

# @time Corrs = GFMC.get_obs_from_accumulator(ObsAcc)
@time Corrs = GFMC.get_obs_from_accumulator_bunching(ObsAcc,4)
Corr_mean = reshape(mean(Corrs), lattice_size, lattice_size, mProj)
Corr_std = reshape(std(Corrs), lattice_size, lattice_size, mProj)

let 
    fig = Figure()
    i = 2
    ax = Axis(fig[1, 1], xlabel=L"$τ$ (imaginary time)", ylabel=L"$\langle S^z_{%$i}S^z_j\rangle$")
    tau =( 0:mProj-1) * dtau
    
    for j in (1:lattice_size)[1:1:end]
        ll = lines!(ax, tau, Corr_mean[i,j,:])
        band!(ax, tau, Corr_mean[i,j,:] - Corr_std[i,j,:], Corr_mean[i,j,:] + Corr_std[i,j,:], alpha = 0.5)
        for corr_i in Corrs
            lines!(ax, tau, reshape(corr_i,lattice_size, lattice_size, mProj)[i,j,:],alpha = 0.2,linewidth = 0.5,linestyle = :solid, color = ll.color)
        end
        hlines!(ax, [corrs_exact[i,j]], color=ll.color, label=L"Exact$$", linestyle=:dash)
    end
    fig
end

##
Bunchs = ([
    std(GFMC.get_obs_from_accumulator_bunching(ObsAcc,i))[2] for i in 1:300
])