# python-bertini

python-bertini is a Python library for numerical algebraic geometry.
Mostly, it's just a wrapper around [Bertini](https://bertini.nd.edu) (sold separately).

## Installing python-bertini

Don't install python-bertini.
It's not ready yet.
If you want something right now, and you're into Julia, you might try [Bertini.jl](https://github.com/PBrdng/Bertini.jl).

If you absolutely must, then clone this repository, cd to the base directory, and run either

```shell
$ conda env create -n bertini -f environment.yml # if you have anaconda
```
or (with your virtualenv activated)
```shell
$ pip install -U -r requirements.txt # if you just want to use a virtualenv
```

Then in either case run
```shell
$ pip install -e .
```