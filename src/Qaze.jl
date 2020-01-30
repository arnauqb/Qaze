module Qaze
export initialize_line, start_lines
#ENV["PYCALL_JL_RUNTIME_PYTHON"]="/home/arnau/Documents/qwind/env/bin/python"
using TOML
include("constants.jl")
include("structures.jl")
include("black_hole.jl")
include("grids.jl")
include("radiation.jl")
include("integrate.jl")
include("streamline.jl")
include("utils.jl")
include("config.jl")


function initialize_line(i, r, wind::WindStruct, is_first_iter)
    T = wind.sed.disk_nt_temperature4(r)^(0.25)
    v_th = thermal_velocity(T)
    if wind.config["wind"]["v_0"] == "thermal"
        v_0 = v_th
    end
    n_0 = wind.config["wind"]["n_0"]
    z_0 = wind.config["wind"]["z_0"]
    line_integ = initialize_line(i, r, z_0, v_0, n_0, v_th, wind, is_first_iter)
    return line_integ
end

function start_lines(wind::WindStruct)
    is_first_iter = true
    for it_num in 1:wind.config["wind"]["iterations"]
        for (i, r) in enumerate(wind.lines_range)
            line = initialize_line(i, r, wind, is_first_iter)
            println("Solving line $i of $(length(wind.lines))")
            solve!(line)
            write_line(wind.config["general"]["save_path"], line.p, it_num)
        end
        write_properties_and_grids(wind.config["general"]["save_path"], wind, it_num)
        is_first_iter = false
    end
end
#wind = initialize(config)
#start_lines(wind)
#T = wind.sed.disk_nt_temperature4(100)^(0.25)
#v_th = thermal_velocity(T)
#println("T: $T, v_th: $v_th")
#line = initialize_line(1, 100, 1, v_th, 2e8, v_th, wind)
#write_line(json_file)

end #module
