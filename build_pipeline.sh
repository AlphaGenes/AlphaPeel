command=$1
if [ $# -eq 0 ] ; then 
    command=none
fi

rm -r build
rm -r dist
python setup.py bdist_wheel

if [[ ! -f src/tinypeel/tinyhouse/Pedigree.py ]] ; then
    echo Pedigree.py file not found. Check that the tinyhouse submodule is up to date
    exit 
fi

# Create python wheel.
rm -r build
rm -r dist
python setup.py bdist_wheel

if [ $command == "install" ] ; then
    pip uninstall AlphaPeel -y
    pip install dist/AlphaPeel-0.0.1-py3-none-any.whl
fi

#Compile manual
 ( cd docs; make latexpdf )


target=AlphaPeel
rm -rf $target
mkdir $target
cp dist/* $target
cp docs/build/latex/AlphaPeel.pdf $target
cp -r example $target
zip -r $target.zip $target

