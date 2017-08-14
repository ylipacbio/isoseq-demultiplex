MY_NOSE_FLAGS?=-v -s

doctest:
	py.test --doctest-modules debarcode/

pylint:
	pylint --errors-only debarcode/*.py

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel

test:
	cluster-to-consensus-primer /pbi/dept/secondary/siv/yli/isoseq/lima/data/flnc.fasta /pbi/dept/secondary/siv/yli/isoseq/lima/data/nfl.fasta /pbi/dept/secondary/siv/smrtlink/smrtlink-alpha/jobs-root/020/020643/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.cluster_report.csv out_dir
