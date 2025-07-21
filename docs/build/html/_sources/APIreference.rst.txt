API Reference
=============
.. toctree::
   :maxdepth: 2
 


CoderData Object
^^^^^^^^^^^^^^^^

.. automodule:: coderdata.download.downloader
   :members: download
   :undoc-members:


.. automodule:: coderdata.utils.utils
      :members: version, list_datasets

.. automodule:: coderdata.cli
   :members: info, check_folder
   :undoc-members:
   :show-inheritance:


Dataset Object
^^^^^^^^^^^^^^

.. autoclass:: coderdata.dataset.dataset.Dataset

.. automodule:: coderdata.dataset.dataset
   :members: load, format, train_test_validate, split_train_other, split_train_test_validate
   :undoc-members:
   :show-inheritance:

.. Hidden
.. ^^^^^^
.. .. automodule:: coderdata.dataset.dataset
   :members: _split_two_way, _create_classes, _balance_data, _filter, _determine_delimiter, _load_file
   :undoc-members:
   :show-inheritance:

.. Build
.. ^^^^^

.. .. automodule:: build.build_dataset
   :members: process_drugs, process_samples, process_omics, process_experiments, process_misc, process_genes, run_docker_cmd, process_docker, run_docker_validate_cmd, run_schema_checker
   :undoc-members:
   :show-inheritance:

.. .. autoclass:: coderdata.build.build_dataset.process_genes

.. .. automodule:: build.beatAML.GetBeatAML
   :members: generate_samples_file
   :undoc-members:
   :show-inheritance: