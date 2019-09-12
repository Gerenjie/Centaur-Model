src = obscode.py obs80.py obs80hc.py leapsec.py
rst = obscode.rst obs80.rst obs80hc.rst leapsec.rst

doc: conf.py index.rst $(rst) $(src)
	sphinx-build-3 . doc
