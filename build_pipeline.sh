rm -r build
rm -r dist
python setup.py bdist_wheel
pip uninstall AlphaPeel -y
pip install dist/AlphaPeel-0.0.1-py3-none-any.whl

