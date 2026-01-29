using FITSIO
using Plots
using LaTeXStrings
using Optim


struct BBClModel
    lmax::Int
    ell::Vector{Int}
    cl_lens::Vector{Float64}
    cl_tens::Vector{Float64}
    cl_sys::Vector{Float64}
end

function BBClModel(lmax::Int, cl_sys::Vector{Float64}, reffpath_r0::String, reffpath_r1::String)
    # Load fiducial power spectra
    cl_r0 = read_planck_cl_as_TT_EE_BB_TE(reffpath_r0)
    cl_r1 = read_planck_cl_as_TT_EE_BB_TE(reffpath_r1)
    # Extract BB component (index 3 for BB)
    cl_lens = cl_r0.BB
    cl_tens = cl_r1.BB
    ell = collect(2:lmax)
    return BBClModel(lmax, ell, cl_lens[ell.+1], cl_tens[ell.+1], cl_sys[ell.+1])
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

function r_iterative_estimator(params::BBClModel, rmin, rmax, rresol,; itr = 1)
    r_ = range(rmin, rmax; length = Int(rresol))
    likelihoods = zeros(length(r_))
    cl_obs = params.cl_lens .+ params.cl_sys
    r_result = 0.0
    ell = params.ell
    for itr_idx in 1:itr
        for (idx, r) in enumerate(r_)
            cl_th = params.cl_lens .+ r .* params.cl_tens
            likelihoods[idx] = _likelihood(ell, cl_obs, cl_th)
        end
        maxid = argmax(likelihoods)
        delta_r = r_[maxid]
        r_result = delta_r
        r_ = range(delta_r - delta_r*(0.5/(itr_idx)), delta_r + delta_r*(0.5/(itr_idx)); length = Int(rresol))
    end
    return r_result
end


@inline function _likelihood(ell, cl_obs, cl_th)
return sum(-1/2 .* (2 .* ell .+ 1) .* (cl_obs ./ cl_th .+ log.(cl_th) .-((2 .* ell .- 1) ./ (2 .* ell .+ 1)) .* log.(cl_obs)))
end