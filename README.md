# JSON Likelihoods for ATLAS SUSY 3L eRJR analysis

The JSON likelihoods are serialized in this folder. This is done by providing a background-only workspace containing the signal/control channels at `BkgOnly.json` as well as patch files for each mass point on the signal phase-space explored in the analysis.

Each [jsonpatch](http://jsonpatch.com/) file follows the format `patch.ERJR_mn2_mn1.json` where `mn2` is the mass of the chargino/second lightest-neutralino and `mn1` is the mass of the lightest supersymmetric particle (LSP).

## Producing signal workspaces

As an example, we use [python jsonpatch](https://python-json-patch.readthedocs.io/en/latest/) here:

```
jsonpatch BkgOnly.json patch.ERJR_500p0_0p0.json > ERJR_500p0_0p0.json
```

## Computing signal workspaces

For example, with [pyhf](https://diana-hep.org/pyhf/), you can do any of the following:

```
pyhf cls BkgOnly.json -p patch.ERJR_500p0_0p0.json

jsonpatch BkgOnly.json patch.ERJR_500p0_0p0.json | pyhf cls

pyhf cls ERJR_500p0_0p0.json
```
# 3L-RJ-mimic-likelihood-validation
