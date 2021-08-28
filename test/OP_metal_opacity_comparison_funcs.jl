# This file defines some functions that are used to compare our implementation of bound-free
# opacity calculations using data from the Opacity Project against calculations of the same
# opacities using the hydrogenic aproximation. This nominally means that we're comparing
# opacities from H I and He II
#
# In the future, it might also be nice to include comparisons for multi-electron data against
# experimental data.
#
# Tests are explicitly NOT defined in this file so that functions defined in this file can be used
# in scripts outside of this framework (to create plots in order to qualitatively assess
# disagreement)
#
# in the future, we should really consolidate this file with opacity_comparison_funcs.jl

function calc_hydrogenic_bf_absorption_coef(λ_vals,  T, ndens_species, species_name;
                                            use_OP_data = false)
    @assert species_name in ["H_I", "He_II"]

    if use_OP_data
        elec_conf_file, cross_sec_file = if species_name == "H_I"
            (joinpath(@__DIR__, "data/TOPbase_electron_config_H.txt"),
             joinpath(@__DIR__, "data/TOPbase_cross_section_H_I.txt"))
        else
            (joinpath(@__DIR__, "data/TOPbase_electron_config_He.txt"),
             joinpath(@__DIR__, "data/TOPbase_cross_section_He_II.txt"))
        end
        Korg.ContinuumOpacity.absorption_coef_bf_TOPBase(λ_vals, T, ndens_species, species_name;
                                                         extrapolation_bc=0.0,
                                                         elec_conf_file = elec_conf_file,
                                                         cross_sec_file = cross_sec_file)
    else
        ν_vals = (Korg.c_cgs*1e8)./λ_vals # Hz
        ndens_div_partition = ndens_species/Korg.partition_funcs[species_name](T)
        ρ = 1.0 # this is unphysical, but it works out fine for a hydrogenic atom
        κ_dflt_approach = if species_name == "H_I"
            H_I_ion_energy = 13.598
            Korg.ContinuumOpacity.H_I_bf.(ndens_div_partition, ν_vals, ρ, T, H_I_ion_energy)
        else
            He_II_ion_energy = 54.418
            Korg.ContinuumOpacity.He_II_bf.(ndens_div_partition, ν_vals, ρ, T, He_II_ion_energy)
        end
        κ_dflt_approach/ρ
    end
end

function _OP_hydrogenic_dflt_λ_vals(species_name)
    first, last = (species_name == "H_I") ? (80,80000) : (25, 40000)
    map((λ) -> Float32(λ), first:1.0:last)
end

function _OP_hydrogenic_λ_comp_intervals(species_name)
    # because the energy state data used by the Opacity Project is somewhat inaccurate (note the
    # disagreement gets worse at higher energy levels), we can only make meaningful comparisons
    # over specific wavelength intervals.

    if species_name == "H_I"
        [80.0       911.0;
         1000.0     3645.0;
         4000.0     8203.0;
         9300.0     14579.0;
         18000.0    22789.0]
    else
        [25.0       10.0^2.35;
         10.0^2.37  10.0^2.95;
         10.0^2.99  10.0^3.3;
         10.0^3.37  10.0^3.55;
         10.0^3.65  10.0^3.75]
    end
end

_OP_hydrogenic_rtol(species_name) = (species_name == "H_I") ? 0.11 : 0.06
