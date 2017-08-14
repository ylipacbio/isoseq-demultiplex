MY_NOSE_FLAGS?=-v -s

doctest:
	py.test --doctest-modules debarcode/

pylint:
	pylint --errors-only debarcode/*.py

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel

c2cp:
	mkdir -p out_dir
	cluster-to-consensus-primer /pbi/dept/secondary/siv/yli/isoseq/lima/data/flnc.fasta /pbi/dept/secondary/siv/yli/isoseq/lima/data/nfl.fasta /pbi/dept/secondary/siv/smrtlink/smrtlink-alpha/jobs-root/020/020643/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.cluster_report.csv out_dir

z2cp:
	mkdir -p out_dir
	zmw-to-consensus-primer out_dir/flnc_z2c.csv out_dir/nfl_z2c.csv out_dir/cluster_dict.csv out_dir
