#!/bin/julia
# make.jl
# 2023-08-05
# curtis
# Make documentation for SimpleSonar

using Documenter, SimpleSonar, DocumenterLaTeX

makedocs(sitename="SimpleSonar.jl",
         format=LaTeX(),
         pages=["index.md"])

