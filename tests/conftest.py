import os


def pytest_addoption(parser):
    group = parser.getgroup('coderdata build data')
    group.addoption('--build-dir', default=os.environ.get('CODERDATA_BUILD_DIR'),
                    help='folder of CoderData build or release files to check with '
                         'tests/test_build_data.py (or set CODERDATA_BUILD_DIR)')
    group.addoption('--previous-summary', default=os.environ.get('CODERDATA_PREVIOUS_SUMMARY'),
                    help='validation_summary.json of an earlier build, to report large changes')
    group.addoption('--build-datasets', default=os.environ.get('CODERDATA_BUILD_DATASETS'),
                    help='comma-separated datasets to check (default: all found in --build-dir)')
