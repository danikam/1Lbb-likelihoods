# JSON Likelihoods for 1Lbb Analysis

The JSON likelihoods are serialized in this folder. This is done by providing a background-only workspace containing the signal/control channels at `BkgOnly.json` as well as patch files for each mass point on the signal phase-space explored in the analysis.

Each [jsonpatch](http://jsonpatch.com/) file follows the format `patch_C1N2_Wh_hbb_m1_m2.json` where [zzz understand what's going on zzz].

## Producing signal workspaces

As an example, we use [python jsonpatch](https://python-json-patch.readthedocs.io/en/latest/) here:

```
jsonpatch BkgOnly.json patch_C1N2_Wh_hbb_550_200.json > C1N2_Wh_hbb_550_200.json
```

## Computing signal workspaces

For example, with [pyhf](https://diana-hep.org/pyhf/), you can do any of the following:

```
pyhf cls BkgOnly.json -p patch_C1N2_Wh_hbb_550_200.json

jsonpatch BkgOnly.json patch_C1N2_Wh_hbb_550_200.json | pyhf cls

pyhf cls C1N2_Wh_hbb_550_200.json
```

