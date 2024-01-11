from setuptools import setup, find_packages

setup(
    name='coderdata',  
    version='0.1.19',   
    author='Jeremy Jacobson',
    author_email='jeremy.jacobson@pnnl.gov',
    description='A package to download, load, and process multiple benchmark multi-omic drug response datasets',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/PNNL-CompBio/candleDataProcessing',
    packages = find_packages(),
    install_requires=[
        'pandas', 
        'requests',
        'numpy'
        ],
    entry_points={
        'console_scripts': [
            'coderdata=coderdata.cli:main',  # CLI command for downloading datasets
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6'
)
