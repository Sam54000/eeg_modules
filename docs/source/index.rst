************
Introduction
************

Here are useful modules for EEG analysis using MNE-python_, mne_bids_, and pyprep_.
These modules have been developped in the specific case of the adult_TBI_attention
project at `Sudhin Shah's lab`_ at Weill Cornell Medicine. This package derives
from the analysis of EEG data during an attention network task (ANT) [#]_ [#]_ recorded in
a large cohort of TBI patients and controls. The EEG data were recorded useg a 128-channel
EGI system (HydroCel Geodesic Sensor Net). Here are also referenced modules for
converting the data from EGI to BIDS format. Feel free to take example from these modules
to adapt them to your own needs.

.. _MNE-python: https://mne.tools/stable/index.html
.. _mne_bids: https://mne.tools/mne-bids/stable/index.html
.. _pyprep: https://pyprep.readthedocs.io/en/latest/
.. _Sudhin Shah's lab: https://radiology.weill.cornell.edu/research/brain-health-imaging-institute/sudhin-shah-laboratory
   
.. note::

   This project is under active development.

.. [#] Fan, J., McCandliss, B. D., Sommer, T., Raz, A., & Posner, M. I. (2002). Testing the efficiency and independence of attentional networks. Journal of cognitive neuroscience, 14(3), 340-347. https://www.mitpressjournals.org/doi/abs/10.1162/089892902317361886
.. [#] Fan, J., McCandliss, B. D., Fossella, J., Flombaum, J. I., & Posner, M. I. (2005). The activation of attentional networks. Neuroimage, 26(2), 471-479. https://www.sciencedirect.com/science/article/pii/S1053811905000416

.. toctree::
   :maxdepth: 1
   
   index
   usage
   APIs
