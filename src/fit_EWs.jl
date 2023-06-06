"""
TODO
"""
function fit_EWs(atm::ModelAtmosphere, linelist, A_X, EWs; abundance_grid=-2:0.1:1, 
                 solar_abundances=default_solar_abundances, max_window_size=5.0,
                 wavelength_resolution=0.01, synthesize_kwargs...)
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

    # sort the linelist by wavelength
    perm = sortperm(linelist, by=l->l.wl)
    linelist = linelist[perm]
    EWs = EWs[perm]

    windows = map(enumerate(linelist)) do (i, line)
        lambda0 = line.wl*1e8
        lambda0 - max_window_size/2: wavelength_resolution : lambda0 + max_window_size/2
    end

    # do an initial synthesis to get chemical equilibrium, continuum alpha, and subspectra
    sol = synthesize(atm, linelist, A_X, windows, hydrogen_lines=false, synthesize_kwargs...)
    nₑs = sol.electron_number_density
    number_densities = sol.number_densities
    subspectra = sol.subspectra
    cntm_alpha = sol.alpha_cntm_itps
    α5 = sol.alpha_5000

    cntm_sol = synthesize(atm, [], A_X, windows, hydrogen_lines=false, 
                         _alpha_cntm_itps=cntm_alpha, _nes=nₑs, _number_densities=number_densities, 
                         _alpha_5000=α5, synthesize_kwargs...)

    absorptions = map(abundance_grid) do X_H
        A_X[element] = solar_abundances[element] + X_H

        # naively rescale the number densities for all species involving element.
        # these aren't fully physically self-consistent.
        ns = copy(number_densities)
        for spec in keys(ns)
            ns[spec] = ns[spec] * 10^(X_H * sum(get_atoms(spec) .== element))
        end

        sol = synthesize(atm, linelist, A_X, windows, hydrogen_lines=false, line_buffer=0.0,
                         _alpha_cntm_itps=cntm_alpha, _nes=nₑs, _number_densities=ns, 
                         _alpha_5000=α5, synthesize_kwargs...)
        1 .- sol.flux ./ cntm_sol.flux
    end

    map(subspectra, EWs) do s, EW
        synthEWs = map(absorptions) do absorption 
            sum(absorption[s]) * wavelength_resolution
        end
        itp = CubicSplines.CubicSpline(synthEWs, abundance_grid, extrapolate=true)
        itp(EW), itp
    end
end
