rm -rf dist
echo "The currrent pkg version is: "
pip show edge_pydb | grep Version
python3 setup.py sdist bdist_wheel
twine upload dist/*
