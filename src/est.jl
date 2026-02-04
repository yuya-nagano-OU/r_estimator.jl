using FITSIO
using Plots
using LaTeXStrings
#using Optim


struct ClBBModel
    lmax::Int
    ell::Vector{Int}
    cl_lens::Vector{Float64}
    cl_tens::Vector{Float64}
    cl_sys::Vector{Float64}
end

function ClBBModel(lmax::Int, cl_sys::Vector{Float64}, reffpath_r0::String, reffpath_r1::String)
    # Load fiducial power spectra
    cl_r0 = read_planck_cl_as_TT_EE_BB_TE(reffpath_r0)
    cl_r1 = read_planck_cl_as_TT_EE_BB_TE(reffpath_r1)
    # Extract BB component (index 3 for BB)
    cl_lens = cl_r0.BB
    cl_tens = cl_r1.BB
    ell = collect(2:lmax)
    return ClBBModel(lmax, ell, cl_lens[ell.+1], cl_tens[ell.+1], cl_sys[ell.+1])
end

function ClBBModel(lmax::Int, cl_sys::Vector{Float64}, cl_lens::Vector{Float64}, cl_tens::Vector{Float64})
    ell = collect(2:lmax)
    return ClBBModel(lmax, ell, cl_lens[ell.+1], cl_tens[ell.+1], cl_sys[ell.+1])
end

function ClBBModel(;lmax::Int, cl_sys::Vector{Float64}, cl_lens::Vector{Float64}, cl_tens::Vector{Float64})
    ell = collect(2:lmax)
    return ClBBModel(lmax, ell, cl_lens[ell.+1], cl_tens[ell.+1], cl_sys[ell.+1])
end

function read_planck_cl_as_TT_EE_BB_TE(fname; hdu_index=2)
    FITS(fname) do f
        hdu = f[hdu_index]  # TableHDU
        TT = read(hdu, "TEMPERATURE")
        EE = read(hdu, "GRADIENT")
        BB = read(hdu, "CURL")
        TE = read(hdu, "G-T")
        N = length(TT)
        ell = collect(0:N-1)
        return (ell=ell, TT=TT, EE=EE, BB=BB, TE=TE)
    end
end

function r_iterative_estimator(params::ClBBModel,; rmin, rmax, rresol, itr = 1)
    r_grid = range(rmin, rmax; length = Int(rresol))
    logL = zeros(length(r_grid))
    cl_obs = params.cl_lens .+ params.cl_sys
    ell = params.ell
    Δr = 0.0
    for itr_idx in 1:itr
        for (idx, r) in enumerate(r_grid)
            cl_th = params.cl_lens .+ r .* params.cl_tens
            logL[idx] = _logL(ell, cl_obs, cl_th)
        end
        maxid = argmax(logL)
        Δr = r_grid[maxid]
        r_grid = range(Δr - Δr*(0.5/(itr_idx)), Δr + Δr*(0.5/(itr_idx)); length = Int(rresol))
    end
    return Δr
end


@inline function _logL(ell, cl_obs, cl_th)
return sum(-1/2 .* (2 .* ell .+ 1) .* (cl_obs ./ cl_th .+ log.(cl_th) .-((2 .* ell .- 1) ./ (2 .* ell .+ 1)) .* log.(cl_obs)))
end

function δr_estimator(params::ClBBModel,; rmin, rmax, rresol)
    r_grid = range(rmin, rmax; length = Int(rresol))
    logL = zeros(length(r_grid))
    cl_obs = params.cl_lens .+ params.cl_sys
    ell = params.ell
    for (idx, r) in enumerate(r_grid)
        cl_th = params.cl_lens .+ r .* params.cl_tens
        logL[idx] = _logL(ell, cl_obs, cl_th)
    end
    L = exp.(logL .- maximum(logL))
    δr, idx = delta_r_at_level(r_grid, L; level = 0.68)
    return δr, L, r_grid, idx
end

function delta_r_at_level(r_grid, L; level = 0.68)
    # Use normalized likelihood for numerical stability.
    dr = step(r_grid)
    cdf = cumsum((L[1:end-1] .+ L[2:end]) .* (dr / 2))
    target = level * cdf[end]
    idx = findfirst(>=(target), cdf)
    return r_grid[idx + 1], idx+1
end