using DataStructures
using OrderedCollections
using DataFrames
using HTTP
using Dates
using TimeZones
using Printf
using JSON3
using SQLite

using PhysConst.UnitConst

const SISREMOROOT = dirname(@__DIR__)
const DATAROOT    = joinpath(SISREMOROOT, "data")
const FIGDIR      = joinpath(SISREMOROOT, "figures")
const JSONROOT    = joinpath(DATAROOT, "json")
const DBPATH      = joinpath(DATAROOT, "ise_data.sqlite")

include("init_logging.jl")
include("energy_data/include_energy_data.jl")
