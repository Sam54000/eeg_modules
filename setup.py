from setuptools import setup
set(
    name="eeg_modules",
    author="Samuel Louviot",
    author_email="sam.louviot@gmail.com",
    packages=["eeg_modules"],
    install_requires=[
        "numpy",
        "pandas",
        "mne",
        "pyprep",
        "mne_bids"
    ]
)
