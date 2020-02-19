pyhf prune --sample MGPy8EG_A14N23LO_C1N2_WZ_500_50_3L_2L7 ERJR_500_50.json > BkgOnly.json
for f in ERJR*.json
do
  echo $f
  jsondiff BkgOnly.json $f > patch.$f
done
