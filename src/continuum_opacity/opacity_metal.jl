# the only reason that we NEED to load in the electron_configurations data is that it's not
# completely clear to me exactly how I would compute certain quantities (i.e. the statistical weight and the ground state energy)

include("../Korg.jl")
using Base
using Interpolations

struct StateID
    iSLP::UInt16
    iLV::UInt8
    function StateID(iSLP, iLV)
        @assert 99 < iSLP < 992 # this upper limit is set based on the data source
        @assert (iSLP % 10) < 2 # trailing digit must be 0 or 1
        @assert 0 < iLV # the upper limit is set based on the data source
        new(iSLP, iLV)
    end
end

# I'm not sure any of this is necessary, but I wanted to document this
Base.:(<)(a::StateID, b::StateID) = (a.iSLP < b.iSLP) || ((a.iSLP == b.iSLP) && (a.iLV < b.iLV))
Base.:(==)(a::StateID, b::StateID) = (a.iSLP == b.iSLP) && (a.iLV == b.iLV)
Base.hash(obj::StateID) = Base.hash((UInt32(obj.iSLP) << 16) + obj.iLV)
Base.show(io::IO, obj::StateID) = print(io, "StateID{iSLP = ", obj.iSLP, ", iLV = ", obj.iLV, "}")

# this isn't strictly necessary, but I wanted to document it
function _quantum_prop_type(state_id::StateID)
    str = String(iSLP)
    S = div(Integer(str[1])-1,2)
    L = Integer(str[2])
    parity = ["even", "odd"][Integer(str[3])+1]
    (S, L, parity)
end
get_S_quantum_num(state_id::StateID) = _quantum_prop_type(state_id)[1]
get_L_quantum_num(state_id::StateID) = _quantum_prop_type(state_id)[2]
get_parity(state_id::StateID) = _quantum_prop_type(state_id)[3]

struct StateProp
    ion_energy_ryd::Float64 # energy with respect to ionization potential in Rydbergs
    excitation_potential_ryd::Float64 # energy above ground state in Rydbergs
    # the following properties can technically be computed from the species name and StateID
    # but the procedure to do so isn't obvious to me
    statistical_weight::UInt8 # technically this can be directly computed
    electron_config::String
end

# parses the atomic_num, number of electrons, iSLP, and iLV entries of the table
function _parse_Z_numelectrons_stateid(str::AbstractString, start_Z::Integer,
                                       start_numElectrons::Integer, start_iSLP::Integer,
                                       start_iLV::Integer, last_iLV::Integer)
    atomic_num = parse(UInt8, str[start_Z:start_numElectrons-1])
    num_electrons = parse(UInt8, str[start_numElectrons:start_iSLP-1])
    iSLP = parse(UInt16, str[start_iSLP:start_iLV-1])
    iLV = parse(UInt16, str[start_iLV:last_iLV])
    (atomic_num, num_electrons, StateID(iSLP, iLV))
end

_get_species_string(atomic_num::Integer, num_electrons::Integer) =
    Korg.atomic_symbols[atomic_num] * "_" * Korg.numerals[atomic_num-num_electrons+1]


const _ELECTRON_CONFIG_TABLE_BREAK =
    " ===================================================================="

"""
    read_electron_configurations(fname)

Parse the provided table produced by an Energy levels query on
[TOPbase](http://cdsweb.u-strasbg.fr/topbase/energy.html).

This expects that the query:
- only includes a single atomic number is included in the query
- orders the output levels in the Level order in each spectroscopic series
- the output data for each level should only include the "Electron
  configuration," "Energy (Ryd) wrt ionization potential," and "Statistical
  weight."
"""
function read_electron_configurations(fname)

    # for simplicity, we're going to make a dict. But the values are ordered
    #data_pairs = Vector{Tuple{StateID,StateProp}}()

    outputs = Vector{Tuple{String,Dict{StateID,StateProp}}}()

    linelist = open(fname) do f
        file_atomic_num = nothing
        current_num_electrons = nothing
        current_dict = nothing

        # confirm that standard file header is present and move the position of the stream past
        # this point
        @assert readline(f) == _ELECTRON_CONFIG_TABLE_BREAK
        @assert readline(f) ==
            "      i NZ NE iSLP iLV iCONF                 E(RYD)      TE(RYD)   gi"
        @assert readline(f) == _ELECTRON_CONFIG_TABLE_BREAK

        for line in eachline(f)
            if line == _ELECTRON_CONFIG_TABLE_BREAK # the last line
                break
            end

            atomic_num, num_electrons, state_id =
                _parse_Z_numelectrons_stateid(line, 9, 11, 14, 19, 22)
            species = _get_species_string(atomic_num, num_electrons)

            if file_atomic_num === nothing
                file_atomic_num = atomic_num
            else
                @assert file_atomic_num == atomic_num
            end

            if current_num_electrons != num_electrons
                current_dict = Dict{StateID,StateProp}()
                push!(outputs, (species, current_dict))
                current_num_electrons = num_electrons
            end

            electron_config = strip(line[24:39])
            ion_energy_ryd = parse(Float32, line[40:51])
            excitation_potential_ryd = parse(Float32, line[53:64])

            # finally parse statistical weight. This is always listed as 1 or 2 digits followed by
            # decimal point followed by a zero. Since statistical weight is always an int, we'll
            # only parse the leading 2 characters.

            if line[66:69] == "****"
                # I think that this is only an issue for 1 highly excited state of Fe
                continue
            end
            statistical_weight = parse(UInt8, line[66:67])

            current_dict[state_id] = StateProp(ion_energy_ryd, excitation_potential_ryd,
                                               statistical_weight, electron_config)
        end
    end
    outputs
end


# The files holding the photoionization cross sections can be massive. With that in mind, I'm
# writing an iterator to handle this. I'm following the design pattern used in the Julia standard
# library for implementing the ``EachLine`` iterable object. We need to acknowledge the Julia Base
# library's MIT license

const _TABLE_BREAK = " ================================================"
const _HEADER_LINE_LENGTH = length(_TABLE_BREAK)

struct EachPhotoIonSubtable{IOT <: IO}
    stream::IOT
    ondone::Function
end

"""
    each_photo_ion_subtable(filename::AbstractString)

Create an iterable `EachPhotoIonSubtable` object that will yield each photoionization subtable from
the input filename. The underlying file is only closed after the `EachPhotoIonSubtable` is garbage
collected. If this is a problem, we probably need to do some light refactoring.
"""
function each_photo_ion_subtable(filename::AbstractString)
    stream = open(filename, "r")
    # confirm that standard file header is present and move the position of the
    # stream past this point
    @assert readline(stream) == _TABLE_BREAK
    @assert readline(stream) == "       I  NZ  NE  ISLP  ILV        E(RYD)      NP"
    @assert readline(stream) == _TABLE_BREAK
    EachPhotoIonSubtable(stream, ()->close(stream))::EachPhotoIonSubtable
end

function Base.iterate(itr::EachPhotoIonSubtable, state=nothing)
    next_line = readline(itr.stream)

    if next_line == _TABLE_BREAK
        return (itr.ondone(); nothing)
    end

    # first parse the subtable's header
    atomic_num, num_electrons, state_id =
        _parse_Z_numelectrons_stateid(next_line, 9, 13, 17, 23, 27)
    species = _get_species_string(atomic_num, num_electrons)
    # ion_energy_ryd is duplicated from the other table. Technically, we don't need the other
    # table. we just need to implement a function to compute the statistical weight from state_id
    # and the species
    ion_energy_ryd = parse(Float32, next_line[28:41])

    num_entries = parse(UInt32, next_line[42:49])

    # this could probably be 32-bit.
    photon_energy_ryd = Array{Float32,1}(undef, num_entries)
    cross_section_MegaBarn = Array{Float32,1}(undef, num_entries)
    for i = 1:num_entries # this properly handles the case when num_entries == 0
        cur_line = readline(itr.stream)
        photon_energy_ryd[i] = parse(Float32, cur_line[1:14])
        cross_section_MegaBarn[i] = parse(Float32, cur_line[15:24])
    end

    # parse the species, and the level
    ((species, state_id, photon_energy_ryd, cross_section_MegaBarn), nothing)
end

Base.eltype(::Type{<:EachPhotoIonSubtable}) =
    Tuple{String, StateID, Array{Float32,1}, Array{Float32,1}}

Base.IteratorSize(::Type{<:EachPhotoIonSubtable}) = SizeUnknown()




"""
    _OP_bf_tabulated_cross_section(λ_vals, T, species_name)

Calculates the LTE cross section (in Mbarn) of an ion species from the Opacity Project tables. This
basically computes the weighted average of all of the energy levels (weighted by the LTE energy 
levels). This is corrected for stimulated emission.

To get the opacity, this needs to be multiplied by the number density of the species and divided by
the total density.
"""
function _OP_bf_tabulated_cross_section(λ_vals, T, species_name)

    split_name = split(species_name, "_")
    @assert length(split_name) == 2
    element_name = split_name[1]
    ion_state = split_name[2]
    @assert ion_state in Korg.numerals

    electron_config_file = joinpath(Korg._data_dir,
                                    string("TOPbase/electron_config/", element_name, ".txt"))
    photo_ion_file = joinpath(Korg._data_dir,
                              string("TOPbase/cross_sections/", species_name, ".txt"))
    println(electron_config_file)
    println(photo_ion_file)

    result = read_electron_configurations(electron_config_file)

    rslt_index = 0
    for i = 1:length(result)
        if result[i][1] == species_name
            rslt_index = i
        end
    end
    @assert rslt_index != 0

    data_dict = result[rslt_index][2]

    photon_energies = (Korg.hplanck_eV * Korg.c_cgs * 1.0e8) ./ λ_vals ./ Korg.RydbergH_eV
    inv_partition_func_val = 1.0/Korg.partition_funcs[species_name](T)
    weighted_average = 0.0 .* photon_energies

    β_Ryd = Korg.RydbergH_eV/(Korg.kboltz_eV * T)

    total_weight = 0.0

    itr = each_photo_ion_subtable(photo_ion_file)
    i = 0
    for (species, state_id, table_photon_energy_ryd, table_cross_section_MBarn) in itr
        i+=1
        if length(table_photon_energy_ryd) == 0
            continue
        end

        if species_name in ["He_I", "C_I", "C_II", "Al_II", "O_II"]
            # In a couple tables, there's a place where the photon energy is listed out of order
            # This seems to be okay given that this is in a section of the table where the
            # cross-section is zero.
            for j in 1:length(table_photon_energy_ryd)-1
                if table_photon_energy_ryd[j+1] < table_photon_energy_ryd[j]
                    @assert j > 1
                    @assert (j + 2) <= length(table_photon_energy_ryd)
                    @assert table_photon_energy_ryd[j+2] >= table_photon_energy_ryd[j+1]
                    @assert table_photon_energy_ryd[j+2] >= table_photon_energy_ryd[j]
                    @assert table_photon_energy_ryd[j+1] >= table_photon_energy_ryd[j-1]
                    @assert table_photon_energy_ryd[j] >= table_photon_energy_ryd[j-1]

                    @assert table_cross_section_MBarn[j+2] == 0.0
                    @assert table_cross_section_MBarn[j+1] == 0.0
                    @assert table_cross_section_MBarn[j] == 0.0
                    @assert table_cross_section_MBarn[j-1] == 0.0

                    # swap the order of the 2 values
                    temp = table_photon_energy_ryd[j]
                    table_photon_energy_ryd[j] = table_photon_energy_ryd[j+1]
                    table_photon_energy_ryd[j+1] = temp
                end
            end
        end

        func = LinearInterpolation(table_photon_energy_ryd, table_cross_section_MBarn,
                                   extrapolation_bc=0.0)
        tmp = data_dict[state_id]
        statistical_weight = tmp.statistical_weight
        energy_ryd = abs(tmp.excitation_potential_ryd)

        weight = inv_partition_func_val*statistical_weight * exp(-energy_ryd * β_Ryd)

        total_weight += weight

        current_cross_section = func.(photon_energies) .* (1.0 .- exp.(-photon_energies.*β_Ryd))
        weighted_average .+= (current_cross_section .* weight)
    end
    weighted_average
end

# Need to standardize some tests and delete all of the following:


using Plots

T = 5000.0
λ_vals = map((λ) -> Float32(λ), 100:1.0:40000)

H_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "H_I")

He_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "He_I")
He_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "He_II")

C_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "C_I")
C_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "C_II")

N_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "N_I")
N_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "N_II")

O_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "O_I")
O_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "O_II")

Mg_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Mg_I")
Mg_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Mg_II")

Al_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Al_I")
Al_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Al_II")

Si_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Si_I")
Si_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Si_II")

Ca_I_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Ca_I")
Ca_II_opacity = _OP_bf_tabulated_cross_section(λ_vals, T, "Ca_II")

p = plot(λ_vals, [H_I_opacity,
                  He_I_opacity, He_II_opacity,
                  C_I_opacity, C_II_opacity,
                  N_I_opacity, N_II_opacity,
                  O_I_opacity, O_II_opacity,
                  Mg_I_opacity, Mg_II_opacity,
                  Al_I_opacity, Al_II_opacity,
                  Si_I_opacity, Si_II_opacity,
                  Ca_I_opacity, Ca_II_opacity],
         labels = ["H I" "He I" "He II" "C I" "C II" "N I" "N II" "O I" "O II" "Mg I" "Mg II" "Al I" "Al II" "Si I" "Si II" "Ca I" "Ca II"],
         ylim = [0.0001, 80], yscale = :log10, xlim = [0, 1e4])

gui()


