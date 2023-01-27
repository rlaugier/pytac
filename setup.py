import sys

from setuptools import setup

setup(name='pytac',
      version='0.1', # defined in the __init__ module
      description='A python tool to work with ESO TAC',
      url='--',
      author='Romain Laugier',
      author_email='romain.laugier@kuleuven.be',
      license='',
      classifiers=[
          'Development Status :: 3 - pre-Alpha',
          'Intended Audience :: Professional Astronomers',
          'Topic :: Astronomical instrumentation engineering :: System control',
          'Programming Language :: Python :: 3.7'
      ],
      packages=['pytac'],
      install_requires=[
          'numpy', 'sympy', 'scipy', 'matplotlib', 'control', 'graphviz'
      ],
      include_package_data=True,
      zip_safe=False)