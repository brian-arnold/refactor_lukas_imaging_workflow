

to create conda environment:

mamba create --name ZZ_calcium
mamba activate ZZ_calcium
mamba install conda-forge::pynrrd conda-forge::matplotlib-base conda-forge::h5py conda-forge::joblib conda-forge::opencv conda-forge::scipy anaconda::scikit-image
pip install tifffile