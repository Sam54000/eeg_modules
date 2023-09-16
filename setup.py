import setuptools
setuptools.setup(
    name="eeg_modules",
    author="Samuel Louviot",
    author_email="sam.louviot@gmail.com",
    url="https://github.com/Sam54000/eeg_modules",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "mne",
        "pyprep",
        "mne_bids"
        "simple_term_menu",
    ]
)
