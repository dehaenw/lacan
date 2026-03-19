from setuptools import setup
from pathlib import Path

setup(
    name='lacan',
    version='1.0.1',    
    description='molecular filter for adjacent fragments',
    long_description=Path("README.md").read_text(),
    long_description_content_type="text/markdown",
    url='https://github.com/dehaenw/lacan',
    author='Wim Dehaen',
    packages=['lacan','lacan.data'],
    install_requires=['rdkit>=2022.03'],
    package_data = {"lacan/data": ["*.pickle"]},
    include_package_data = True,
    extras_require={
        'dev': ['pytest'],
        'docs': ['sphinx', 'sphinx-rtd-theme'],
        'notebooks': ['scikit-learn','numpy','py3Dmol'],
    },
)


