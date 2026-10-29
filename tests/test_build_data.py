"""
Data checks on a real CoderData build or release, run as tests.

These are the checks `coderbuild/build_all.py --validate` runs, exposed to
pytest so any folder of CoderData files can be checked on its own: a local
build (local/all_files_dir) or a downloaded release.

    pytest tests/test_build_data.py --build-dir local/all_files_dir
    pytest tests/test_build_data.py --build-dir data/ --build-datasets beataml,hcmi
    CODERDATA_BUILD_DIR=local/all_files_dir pytest tests/test_build_data.py

Each dataset is one test, plus one for cross-dataset checks. Errors fail the
test; warnings are reported as pytest warnings. Without --build-dir the tests
are skipped. The checks themselves are unit tested in test_validate_build.py.
"""
import os
import sys
import warnings

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'coderbuild', 'utils'))

import validate_build as vb  # noqa: E402


class BuildDataWarning(UserWarning):
    """A validation warning: worth a look, not a failure."""


def _options(config):
    build_dir = config.getoption('--build-dir')
    datasets = config.getoption('--build-datasets')
    return build_dir, (datasets.split(',') if datasets else None)


def pytest_generate_tests(metafunc):
    if 'dataset' not in metafunc.fixturenames:
        return
    build_dir, datasets = _options(metafunc.config)
    if not build_dir:
        metafunc.parametrize('dataset', [pytest.param(None, marks=pytest.mark.skip(
            reason='no build to check: pass --build-dir or set CODERDATA_BUILD_DIR'))])
        return
    metafunc.parametrize('dataset', datasets or vb.discover_datasets(build_dir))


@pytest.fixture(scope='session')
def build_report(pytestconfig):
    build_dir, datasets = _options(pytestconfig)
    if not build_dir:
        pytest.skip('no build to check: pass --build-dir or set CODERDATA_BUILD_DIR')
    if not os.path.isdir(build_dir):
        pytest.fail(f'--build-dir {build_dir} is not a folder')
    summary = pytestconfig.getoption('--previous-summary')
    previous = vb.load_summary(summary) if summary else None
    report, _ = vb.validate(build_dir, datasets, previous)
    return report


def _check(report, scope):
    for check, _, message in [w for w in report.warnings if w[1] == scope]:
        warnings.warn(f'{check}: {message}', BuildDataWarning)
    errors = [f'{check}: {message}' for check, ds, message in report.errors if ds == scope]
    assert not errors, '\n'.join(errors)


def test_dataset(build_report, dataset):
    _check(build_report, dataset)


def test_cross_dataset(build_report):
    _check(build_report, 'all')
