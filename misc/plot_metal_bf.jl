using Plots
using Printf: @sprintf
include("../src/Korg.jl")
include("../test/OP_metal_opacity_comparison_funcs.jl")

function plot_general(; T = 5000.0)
    λ_vals = map((λ) -> Float32(λ), 100:1.0:40000)

    H_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "H_I")

    He_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "He_I")
    He_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "He_II")

    C_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "C_I")
    C_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "C_II")

    N_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "N_I")
    N_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "N_II")

    O_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "O_I")
    O_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "O_II")

    Mg_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Mg_I")
    Mg_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Mg_II")

    Al_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Al_I")
    Al_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Al_II")

    Si_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Si_I")
    Si_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Si_II")

    Ca_Icross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Ca_I")
    Ca_IIcross_sec = Korg.ContinuumOpacity._OP_bf_tabulated_cross_section(λ_vals, T, "Ca_II")

    p = plot(λ_vals, [H_Icross_sec,
                      He_Icross_sec, He_IIcross_sec,
                      C_Icross_sec, C_IIcross_sec,
                      N_Icross_sec, N_IIcross_sec,
                      O_Icross_sec, O_IIcross_sec,
                      Mg_Icross_sec, Mg_IIcross_sec,
                      Al_Icross_sec, Al_IIcross_sec,
                      Si_Icross_sec, Si_IIcross_sec,
                      Ca_Icross_sec, Ca_IIcross_sec],
             labels = ["H I" "He I" "He II" "C I" "C II" "N I" "N II" "O I" "O II" "Mg I" "Mg II" "Al I" "Al II" "Si I" "Si II" "Ca I" "Ca II"],
             ylim = [0.0001, 80], yscale = :log10, xlim = [0, 1e4])

    gui()
end

function _pos_value_extrema(arr)
    min_val = floatmax(eltype(arr))
    max_val = 0.0

    for elem in arr
        if !(elem > 0)
            continue
        end
        min_val = min(min_val, elem)
        max_val = max(max_val, elem)
    end
    (min_val,max_val)
end



function _plot_x_intervals!(p, x_intervals)
    y_lim = ylims(p)
    for i in 1:size(x_intervals)[1]
        cur_x = [x_intervals[i,1], x_intervals[i,2]]
        cur_y1 = [y_lim[1], y_lim[1]]
        cur_y2 = [y_lim[2], y_lim[2]]
        plot!(p, cur_x, cur_y1, fillrange = cur_y2, fillalpha = 0.10,
              ylim = y_lim, fillcolor = :grey,
              label = "")
    end
end



function plot_H_I_bf(T = 7800.0, nH_I = 3.0e16)

    λ_vals = map((λ) -> Float32(λ), 80:1.0:80000)

    α_H_I_bf_OP = calc_hydrogenic_bf_absorption_coef(λ_vals, T, nH_I, "H_I"; use_OP_data = true)
    α_H_I_bf_dflt = calc_hydrogenic_bf_absorption_coef(λ_vals, T, nH_I, "H_I"; use_OP_data = false)

    p = plot()

    ylims = extrema((_pos_value_extrema(α_H_I_bf_OP)..., _pos_value_extrema(α_H_I_bf_dflt)...))
    p = plot!(p, λ_vals, [α_H_I_bf_OP, α_H_I_bf_dflt],
              labels = ["Opacity Project" "Reference (Hydrogenic)"],
              yscale = :log10, xscale = :log10, ylims = ylims,
              ylabel = "absorption coef, α (cm⁻¹)",
              xlabel = "λ (Å)")

    vline!(p, [911.75, 3646, 8204, 14580, 22790, 32820], label = "known binding energies")

    λ_comp_intervals = _OP_hydrogenic_λ_comp_intervals("H_I")
    _plot_x_intervals!(p, λ_comp_intervals)
    gui()
end



function plot_He_II_bf(T = 7800.0, nHe_II = 3.0e16)

    λ_vals = map((λ) -> Float32(λ), 25:1.0:40000)

    α_He_II_bf_OP = calc_hydrogenic_bf_absorption_coef(λ_vals,  T, nHe_II, "He_II";
                                                       use_OP_data = true)
    α_He_II_bf_dflt = calc_hydrogenic_bf_absorption_coef(λ_vals,  T, nHe_II, "He_II";
                                                         use_OP_data = false)

    λ_comp_intervals = _OP_hydrogenic_λ_comp_intervals("He_II")

    p = plot()

    ylims = extrema((_pos_value_extrema(α_He_II_bf_OP)..., _pos_value_extrema(α_He_II_bf_dflt)...))
    p = plot!(p, λ_vals, [α_He_II_bf_OP, α_He_II_bf_dflt],
              labels = ["Opacity Project" "Reference (Hydrogenic)"],
              yscale = :log10, xscale = :log10, ylims = ylims,
              ylabel = "absorption coef, α (cm⁻¹)",
              xlabel = "λ (Å)")

    _plot_x_intervals!(p, λ_comp_intervals)
    gui()
end



