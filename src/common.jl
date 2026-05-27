using Revise, Infiltrator
using Statistics, CSV, DataFrames, YAXArrays, TOML
using ADRIA

const ROOT_PATH = dirname(@__DIR__)
_config_path = joinpath(ROOT_PATH, "config.toml")
_config = TOML.parsefile(_config_path)
pd_config = _config["PathwayDiversityAnalysis"]
