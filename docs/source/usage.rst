Usage
=====
.. _Installation:

Installation
____________
The package relies on MNE-python_ which is quite heavy and can take a very long time to install.
I would suggest to install MNE-python_ first and then install this package.
To `install MNE-python`_ you have several choices that are very well documented on their website:
https://mne.tools/stable/install/index.html

To install the package:

1. Clone the repository on your computer

.. code-block:: console

   $ git clone https://github.com/Sam54000/eeg_modules.git

2. Move to the newly cloned repository ``eeg_modules`` on your computer

.. code-block:: console

   $ cd eeg_modules

3. Install the package using ``pip``

.. code-block:: console

   $ python -m pip install -e .

.. _MNE-python: https://mne.tools/stable/index.html
.. _install MNE-python: https://mne.tools/stable/install/index.html

.. _Usage:

Usage
_____
The biggest part of the package is the ``bids_conversion`` module.

