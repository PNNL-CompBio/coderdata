API Reference
=============
.. toctree::
   :maxdepth: 2
 

CoderData Object
^^^^^^^^^^^^^^^^

.. automodule:: coderdata.download.downloader
   :members: download
   :no-undoc-members:


.. automodule:: coderdata.utils.utils
      :members: version, list_datasets
      :no-undoc-members:
      :show-inheritance:


.. automodule:: coderdata.cli
   :members: info, check_folder
   :no-undoc-members:
   :show-inheritance:


Dataset Object
^^^^^^^^^^^^^^

.. autoclass:: coderdata.dataset.dataset.Dataset

.. autofunction:: coderdata.dataset.dataset.Dataset.save

.. automodule:: coderdata.dataset.dataset
   :members: load, format, train_test_validate, split_train_other, split_train_test_validate
   :show-inheritance:

.. automodule:: coderdata.utils.stats
   :members: summarize_response_metric, plot_response_metric, plot_2D_respones_metric
   :no-undoc-members:
   :show-inheritance:
