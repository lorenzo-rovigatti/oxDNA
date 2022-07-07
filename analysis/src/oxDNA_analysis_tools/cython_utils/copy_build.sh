#cd cython_utils
python3 setup.py build_ext --inplace
cp *.so ../
#cd ..
#python3 m_mean.py
