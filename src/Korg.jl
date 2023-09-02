module Korg
    _data_dir = joinpath(@__DIR__, "../data") 

    include("CubicSplines.jl")             # 1D cubic Splines
    include("constants.jl")                # physical constants
    include("atomic_data.jl")              # symbols and atomic weights
    include("isotopic_abundances.jl")      # self-explanatory
    include("species.jl")                  # types for chemical formulae and species
    include("linelist.jl")                 # parse linelists, define Line type
    include("line_absorption.jl")          # opacity, line profile, voigt function
    include("hydrogen_line_absorption.jl") # hydrogen lines get special treatment
    include("autodiffable_conv.jl")        # wrap DSP.conv to be autodiffable
    include("read_statmech_quantities.jl") # approximate Us, Ks, chis
    include("statmech.jl")                 # statistical mechanics, molecular equilibrium
    include("atmosphere.jl")               # parse model atmospheres
    include("RadiativeTransfer/RadiativeTransfer.jl")        # radiative transfer formal solution
    include("utils.jl")                                # functions to apply LSF, vac<->air wls, etc.
    include("ContinuumAbsorption/ContinuumAbsorption.jl")   # Define continuum absorption functions.
    include("synthesize.jl")                                # top-level API
    include("fit_EWs.jl")                  # fit EWs
    include("fit.jl")                      # routines to infer stellar params from data

    using .Fit: find_best_fit_params
    export synthesize, read_linelist, read_model_atmosphere, interpolate_marcs, format_A_X, 
        find_best_fit_params
end # module
