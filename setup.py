from setuptools import setup

setup(
    name='lacan',
    version='0.0.1alpha',    
    description='filter for adjacent fragments',
    url='https://github.com/dehaenw/lacan',
    author='Wim Dehaen',
    packages=['lacan','lacan.data'],
    install_requires=['rdkit>=2022.03'],
    package_data = {"lacan/data": ["*.pickle"]},
    include_package_data = True,
)

