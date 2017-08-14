MY_NOSE_FLAGS?=-v -s

doctest:
	py.test --doctest-modules debarcode/

pylint:
	pylint --errors-only debarcode/*.py

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
