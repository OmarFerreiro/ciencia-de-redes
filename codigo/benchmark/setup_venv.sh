python3 -m venv $HOME/bench/venv 2>&1 | tail -3
$HOME/bench/venv/bin/pip install --quiet numpy scipy 2>&1 | tail -5
$HOME/bench/venv/bin/python -c "import numpy, scipy; print('numpy', numpy.__version__, 'scipy', scipy.__version__)"
du -sh $HOME/bench/venv
