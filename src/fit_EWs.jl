"""
TODO
"""
function fit_EWs(atm::ModelAtmosphere, linelist, A_X, EWs; abundance_grid=-2:0.1:1, 
                 solar_abundances=default_solar_abundances, window_size=5.0,
                 wavelength_resolution=0.01, synthesize_kwargs...)
    atoms = map(linelist) do line
        get_atoms(line.species)
    end
    if !all(length.(atoms) .== 1)
        throw(ArgumentError("Can't fit the EWs of a molecular line"))
    end
    element = atoms[1][1]
    if !all(atoms .== [element])
        throw(ArgumentError("Can't fit the EWs of lines from different elements"))
    end

    window_size *= 1e-8
    wavelength_resolution *= 1e-8
    windows = map(linelist) do line
        line.wl - window_size/2: wavelength_resolution : line.wl + window_size/2
    end

    wl_lb_ind = 1 # the index into α of the lowest λ in the current wavelength range
    subspectra = []
    for λs in wl_ranges
        wl_inds = wl_lb_ind : wl_lb_ind + length(λs) - 1
        push!(subspectra, wl_inds)
        wl_lb_ind += length(λs)
    end

    absorptions = map(abundance_grid) do X_H
        A_X[element] = solar_abundances[element] + X_H
        sol = synthesize(atm, linelist, A_X, windows, synthesize_kwargs...)
        cntm_sol = synthesize(atm,[], A_X, windows, synthesize_kwargs...)
        1 .- sol.flux ./ cntm_sol.flux
    end

    map(sol.subspectra, EWs) do s, EW
        synthEWs = map(absorptions) do absorption 
            sum(absorption[s]) * wavelength_resolution
        end
        itp = CubicSplines.CubicSpline(synthEWs, abundance_grid)
        itp(EW), itp
    end
end
