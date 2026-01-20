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

