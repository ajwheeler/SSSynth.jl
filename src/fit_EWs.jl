function fit_EWs(atm::ModelAtmosphere, linelist, A_X, EWs; kwargs...)
    example_ind = 35

    example_line = linelist[example_ind]
    A_X, itps, example_sol = fit_EWs_exactly(atm, [example_line], A_X, EWs[example_ind]; 
                                             abundance_grid=-5:0.1:3,kwargs...)
    curve_of_growth = itps[1]

    line_center_ind = argmin(abs.(example_sol.wavelengths .- example_line.wl*1e8))
    example_sol.alpha[line_center_ind:line_center_ind]
    τ = similar(example_sol.alpha[:, line_center_ind])
    RadiativeTransfer.BezierTransfer.compute_tau_bezier!(τ, get_zs(atm), example_sol.alpha[:, line_center_ind:line_center_ind])
    contribution = blackbody.(get_temps(atm), example_line.wl) .* exp.(-τ)
    #formation_layer_ind = argmax(contribution)
    formation_layer_ind = 41
    formation_layer = atm.layers[formation_layer_ind]

    #display([τ contribution])
    #println(formation_layer_ind)
    #println(formation_layer)

    α_eachline = reverse(ContinuumAbsorption.total_continuum_absorption(
        reverse([c_cgs/line.wl for line in linelist]),
        formation_layer.temp,
        example_sol.electron_number_density[formation_layer_ind],
        Dict([k=>example_sol.number_densities[k][formation_layer_ind] 
              for k in keys(example_sol.number_densities)]),
        Korg.default_partition_funcs))

    α_cntm_example = example_sol.alpha_cntm_itps[formation_layer_ind](example_line.wl)
    θ = log10(ℯ)/(kboltz_eV * formation_layer.temp) 

    # Gray equation 16.6 
    A_Xs = curve_of_growth.(EWs)
    ΔA_Xs = map(linelist, α_eachline) do line, α
        line.log_gf - example_line.log_gf 
        + log10(line.wl/example_line.wl) 
        #- log10(α/α_cntm_example)
        #- θ * (line.E_lower - example_line.E_lower)
    end
    A_Xs .- ΔA_Xs
end

"""
TODO
"""
function fit_EWs_exactly(atm::ModelAtmosphere, linelist, A_X, EWs; abundance_grid=-2:0.1:1, 
                        solar_abundances=default_solar_abundances, max_window_size=5.0,
                        wavelength_resolution=0.01, synthesize_kwargs...)
    @assert issorted(linelist, by=l->l.wl)

    atoms = map(linelist) do line
        get_atoms(line.species)
    end
    if !all(length.(atoms) .== 1)
        throw(ArgumentError("Can't fit the EWs of a molecular line"))
    end
    element = atoms[1][1]
    if !all(first.(atoms) .== [element])
        throw(ArgumentError("Can't fit the EWs of lines from different elements"))
    end

    windows = map(enumerate(linelist)) do (i, line)
        lambda0 = line.wl*1e8

        prev_line = i == 1 ? lambda0 - max_window_size : linelist[i-1].wl*1e8
        ll = max(lambda0 - max_window_size/2, (prev_line+lambda0)/2)

        next_line = i == length(linelist) ? lambda0+max_window_size : linelist[i+1].wl*1e8
        ul = min(lambda0 + max_window_size/2, (next_line+lambda0)/2) - wavelength_resolution

        ll : wavelength_resolution : ul
    end

    # do an initial synthesis to get chemical equilibrium, continuum alpha, and subspectra
    sol = synthesize(atm, linelist, A_X, windows; electron_number_density_warn_threshold=1e10,
                     line_buffer=0.0, hydrogen_lines=false, synthesize_kwargs...)
    nₑs = sol.electron_number_density
    number_densities = sol.number_densities
    subspectra = sol.subspectra
    cntm_alpha = sol.alpha_cntm_itps
    α5 = sol.alpha_5000

    cntm_sol = synthesize(atm, [], A_X, windows; hydrogen_lines=false, line_buffer=0.0,
                         electron_number_density_warn_threshold=1e10,
                         _alpha_cntm_itps=cntm_alpha, _nes=nₑs, _number_densities=number_densities, 
                         _alpha_5000=α5, synthesize_kwargs...)

    absorptions = map(abundance_grid) do X_H
        # naively rescale the number densities for all species involving element.
        # these aren't fully physically self-consistent.
        ns = copy(number_densities)
        for spec in keys(ns)
            ns[spec] = ns[spec] * 10^((X_H+solar_abundances[element]-A_X[element]) * sum(get_atoms(spec) .== element))
        end

        sol = synthesize(atm, linelist, A_X, windows; hydrogen_lines=false, line_buffer=0.0,
                         electron_number_density_warn_threshold=1e10,
                         _alpha_cntm_itps=cntm_alpha, _nes=nₑs, _number_densities=ns, 
                         _alpha_5000=α5, synthesize_kwargs...)
        1 .- sol.flux ./ cntm_sol.flux
    end

    ps = map(subspectra, EWs) do s, EW
        synthEWs = map(absorptions) do absorption 
            sum(absorption[s]) * wavelength_resolution
        end
        itp = CubicSplines.CubicSpline(synthEWs, abundance_grid, extrapolate=true)
        itp(EW), itp
    end

    first.(ps), last.(ps), sol
end
