rm -r build
rm -r dist
python setup.py bdist_wheel
# pip uninstall AlphaPeel -y
# pip install dist/AlphaPeel-0.0.1-py3-none-any.whl


#Compile manual
 ( cd docs; make latexpdf )


target=AlphaPeel
rm -rf $target
mkdir $target
cp dist/* $target
cp docs/build/latex/AlphaPeel.pdf $target
cp -r example $target
zip -r $target.zip $target

