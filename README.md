# DICE.jl-notebooks
Jupyter notebooks for working through [DICE.jl](https://github.com/Libbum/DICE.jl) examples.

# Usage

To get things installed and ready, you'll need to set up the project.
Do this via the following commands:

```bash
$ git clone git@github.com:Libbum/DICE.jl-notebooks.git
$ cd DICE.jl-notebooks
$ julia
julia> ]
(v1.1) pkg> activate .
(DICE.jl-notebooks) pkg> instantiate
(DICE.jl-notebooks) pkg> precompile
$ jupyter lab
```

The final command can also be `jupyter notebook` if you don't have `lab` installed.

# Status

At the moment, the notebooks are mostly comparisons of DICE.jl to the GAMS versions of DICE they correspond to.
Source code and output of each GAMS run can be found in the `GAMS` folder for you to run and compare yourself if you have a license to do so.

✔️   v2013R Vanilla and v2016R beta both see a 1:1 correspondence with the GAMS output, thus can be used interchangeably at present.

✔️  The v2013R Rocky Road version is also verified. There are a few deviations from 1:1 here, since GAMS is using truncated 32bit floats, whereas DICE.jl used 64bit floats. Additionally, there are bugs in the GAMS version of the `Stern` and `SternCalibrated` scenarios. We consider patched versions of both runs here, they can be found in the `GAMS` directory.

❌  v2016R2 is still in development, files here are solely here to aid in that process.
