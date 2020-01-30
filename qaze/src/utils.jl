using JSON
using Printf
#include("structures.jl")
#using ..structures
function drawline(x1::Int64, y1::Int64, x2::Int64, y2::Int64)
    x=x1
    y=y1
    dx=abs(x2-x1)
    dy=abs(y2-y1)
    s1=sign(x2-x1)
    s2=sign(y2-y1)
    swap=false
    length = max(dx,dy) + 1
    results = zeros(Int64, length, 2)
    results[1,1] = x1
    results[1,2] = y1
    if(dy>dx)
        temp=dx
        dx=dy
        dy=temp
        swap=true
    end
    p=2*dy-dx
    for i in 2:dx
        while(p>=0)
            p=p-2*dx
            if(swap)
                x+=s1
            else
                y+=s2
            end
        end
        p=p+2*dy;
        if(swap)
            y+=s2
        else
            x+=s1
        end
        results[i, 1] = x
        results[i, 2] = y
    end
    results[length,1] = x2
    results[length,2] = y2
    return results
end

function get_index(array, value)
    return max(searchsortedlast(array, value), 1)
end
function initialize_json(json_file, wind::Wind)
    data = Dict()
    data["metadata"] = Dict()
    data["metadata"]["M"] = wind.bh.M
    data["metadata"]["mdot"] = wind.bh.mdot
    data["metadata"]["spin"] = wind.bh.spin
    data["metadata"]["eta"] = wind.bh.eta

    data["metadata"]["disk_r_min"] = wind.config["disk"]["inner_radius"]
    data["metadata"]["disk_r_max"] = wind.config["disk"]["outer_radius"]

    data["metadata"]["grid_r_min"] = wind.grids.r_range[1]
    data["metadata"]["grid_z_min"] = wind.grids.z_range[1]
    data["metadata"]["grid_r_max"] = wind.grids.r_range[end]
    data["metadata"]["grid_z_max"] = wind.grids.z_range[end]
    data["metadata"]["n_grid_r"] = wind.grids.n_r
    data["metadata"]["n_grid_z"] = wind.grids.n_z
    data["metadata"]["n_grid_disk"] = wind.grids.n_disk
    data["metadata"]["n_lines"] = wind.grids.n_lines
    data["metadata"]["r_range"] = wind.grids.r_range
    data["metadata"]["z_range"] = wind.grids.z_range

    data["metadata"]["f_uv"] = wind.radiation.f_uv
    data["metadata"]["f_x"] = wind.radiation.f_x
    data["metadata"]["rho_shielding"] = wind.config["wind"]["n_shielding"]
    data["metadata"]["r_in"] = wind.config["wind"]["initial_radius"]
    data["metadata"]["r_out"] = wind.config["wind"]["final_radius"]
    for it in 1:wind.config["wind"]["iterations"]
        data[@sprintf("it_%02d", it)] = Dict(
                "lines" => [],
                "lines_escaped" => [],
                "grids" => Dict()
        )
    end
    open(json_file, "w") do f
        JSON.print(f, data)
    end

end

function write_line(json_file, line::Streamline, it_num) 
    line_data = Dict(
        "r" => line.u_hist[:,1],
        "z" => line.u_hist[:,2],
        "v_r" => line.u_hist[:,3],
        "v_z" => line.u_hist[:,4],
        "v_T" => sqrt.(line.u_hist[:,3].^2 + line.u_hist[:,4].^2),
        "n" => line.n_hist,
        "fm" => line.fm_hist,
        "xi" => line.xi_hist
    )
    data = open(json_file, "r") do f
        data = JSON.parse(f)
    end
    open(json_file, "w") do f
        push!(data[@sprintf("it_%02d", it_num)]["lines"], line_data)
        JSON.print(f, data)
    end
end

function write_properties_and_grids(json_file, wind::Wind, it_num)
    #=
    properties = Dict(
        "mdot_w_gs" => wind.mdot_w,
        "mdot_w_msunyr" => wind.mdot_w / constants.M_SUN * constants.YEAR_TO_SEC,
        "kin_lumin" => wind.kinetic_luminosity,
        "angle" => wind.angle,
        "terminal_velocity" => wind.v_terminal,
    )

    print("updating tau_x grid...")
    wind.grid.tau_x_grid.update_all()
    print("updating ionization grid...")
    wind.grid.ionization_grid.update_all()
    =#
    grids = Dict(
        "density" => wind.grids.density,
        "force_multiplier" => wind.grids.fm,
    )
    data = open(json_file, "r") do f
        data = JSON.parse(f)
    end
    open(json_file, "w") do f
        #data[@sprintf("it_%02d", it_num)]["properties"] = properties
        data[@sprintf("it_%02d", it_num)]["grids"] = grids 
        JSON.print(f,data)
    end
    #with open(os.path.join(folder_name, "wind_properties.csv"), "w") as f:
    #    f.write(json.dumps(properties))
end