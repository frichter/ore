

git add -u && git commit -m "Version update accepting lower bounds to AF input" && git push

## final commands run: 
## need to update version EVERYTIME in version.py
# confirm readme rst is parsed correctly:
python3 setup.py check --restructuredtext
# confirm previous builds are archived
ls dist/ore*
mv dist/ore* dist/archive/
# create new builds
python3 setup.py sdist
python3 setup.py bdist_wheel
twine upload dist/ore*
