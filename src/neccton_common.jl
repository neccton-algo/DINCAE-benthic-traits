# Common functions to all steps
# in the preparing, training and validation of the reconstruction

using NCDatasets
using Interpolations

# Resolution
dlon = dlat = 1/16

# lon/lat grid
gridlon = 27.4:dlon:32
gridlat = 41.8:dlat:46.8

# File location
basedir = expanduser("~/Data/NECCTON")
datadir = joinpath(basedir,"General","Datasets")
auxdir = joinpath(basedir,"Auxdata_$(1/dlon)")


split_fname = joinpath(basedir,"split.nc")

moddir = joinpath(basedir,"PourSeverine")
moddir = "/home/abarth/Data/NECCTON/PourSeverine"

station_fname = joinpath(datadir,"Stations","stations.csv")
CWM_response_fname = joinpath(datadir,"Community Weighted Mean (sites x traits)","CWM_response.csv")
env_matrix_fname = joinpath(datadir,"Environnement","matrice environnement.txt")


function neccton_load_mask()

    #fnames0 = sort(glob("Cl*_1d_*_*_btrc_T_*-*.nc",moddir))
    fnames = sort(glob("Cl1992_1d_*_*_grid_T_*-*.nc",moddir))
    NCDataset(fnames[1]) do ds_T
        sbtemper = nomissing(ds_T["sbtemper"][:],NaN)
        mask = .!(isnan.(sbtemper[:,:,1]) .|| sbtemper[:,:,1] .== 0)
        return mask
    end
end


"""
    lon,lat = neccton_load_coord(fname)

Load longitude and latitude from NEMO file `fname`.
"""
function neccton_load_coord(fname)
    NCDataset(fname) do ds
        lon = ds["nav_lon"][:]
        lat = ds["nav_lat"][:]

        lon = allowmissing(lon)
        lon[lon .== -1] .= missing

        lat = allowmissing(lat)
        lat[lat .== -1] .= missing

        lon1 = mapslices(s -> first(skipmissing(s)),lon,dims=2)[:,1]
        lat1 = mapslices(s -> first(skipmissing(s)),lat,dims=1)[1,:]

        return lon1,lat1
    end
end

"""
    df = df_load(station_fname,CWM_response_fname)

Load dataframe from file `station_fname` and `CWM_response_fname` and join
on station name.
"""
function df_load(station_fname,CWM_response_fname)
    station = CSV.read(station_fname,DataFrame);
    resp = CSV.read(CWM_response_fname,DataFrame,decimal=',');
    resp = rename(resp, :Column1 => :sta)
    df = innerjoin(resp, station, on = :sta)
    return df
end


"""
    saveinterp((lon,lat),field,(gridlon,gridlat),varname,interp_fname)

Interpolate the field `field` defined over the grid `(lon,lat)`
on the grid `(gridlon,gridlat)`. The result is saved in the NetCDF files
`interp_fname` under the same `varname`.
`lon`, `lat`, `gridlon`, `gridlat` are all vectors.
"""
function saveinterp((lon,lat),SS2,(gridlon,gridlat),varname,interp_fname)
    @info "interpolate"
    itp = interpolate((lon,lat), SS2, Gridded(Linear()));
    SSi = itp(gridlon,gridlat);

    @show extrema(SSi)

    ds = Dataset(interp_fname,"c")
    # Dimensions

    ds.dim["lon"] = length(gridlon)
    ds.dim["lat"] = length(gridlat)

    # Declare variables

    nclon = defVar(ds,"lon", Float64, ("lon",))
    nclon.attrib["units"] = "degrees_east"
    nclon.attrib["standard_name"] = "longitude"
    nclon.attrib["long_name"] = "longitude"

    nclat = defVar(ds,"lat", Float64, ("lat",))
    nclat.attrib["units"] = "degrees_north"
    nclat.attrib["standard_name"] = "latitude"
    nclat.attrib["long_name"] = "latitude"

    ncvar = defVar(ds,lowercase(varname), Float32, ("lon", "lat"))
    ncvar.attrib["_FillValue"] = Float32(9.96921e36)
    ncvar.attrib["missing_value"] = Float32(9.96921e36)
    ncvar.attrib["long_name"] = varname


    # Define variables

    nclon[:] = gridlon
    nclat[:] = gridlat
    ncvar[:,:] = SSi

    close(ds)
end



function plotmap(bathname = joinpath(ENV["HOME"],"projects","Julia","DIVAnd-example-data","Global","Bathymetry","gebco_30sec_4.nc");
                  patchcolor = [.8,.8,.8], coastlinecolor = nothing)

    xl = xlim()
    yl = ylim()
    # work-around
    xl = xl[1]:0.1:xl[2]
    yl = yl[1]:0.1:yl[2]

    bx,by,b = DIVAnd.extract_bath(bathname,true,xl,yl)
    if patchcolor !== nothing
        contourf(bx,by,b', levels = [-1e5,0],colors = [patchcolor])
    end

    if coastlinecolor !== nothing
        contour(bx,by,b', levels = [-1e5,0],colors = coastlinecolor, linestyles = "-")
    end
end

"""
    set_aspect_ratio()

Fixes the aspect ratio of a plot.
"""
function set_aspect_ratio()
    ax = gca()
    as = cosd(mean(ylim()))
    ax.set_aspect(1/as)
end
