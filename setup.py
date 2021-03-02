from setuptools import setup

with open('README.rst', 'r', encoding='utf-8') as f:
    long_description = f.read()
with open('version.txt', 'r') as f:
    version = f.read().rstrip()

setup(
   name='isosolve',
   version=version,
   description='integrative framework for measurements of isotope labeling',
   keywords='isotope labeling, least squares, symbolic solution',
   license='GNU General Public License v2 or later (GPLv2+)',
   long_description=long_description,
   author='Serguei Sokol',
   author_email='sokol@insa-toulouse.fr',
   url='https://github.com/sgsokol/isosolve',
   packages=['isosolve'],
   py_modules=['isosolve'],
   package_dir={'isosolve': '.'},
   package_data={
        'isosolve': ['version.txt', 'doc/*', 'example/*'],
   },
   install_requires=['sympy', 'markdown', 'numpy', 'scipy', 'pandas', 'nlsic'],
   dependency_links=[
      'https://github.com/didix21/mdutils/tarball/master',
   ],
   entry_points={
      'console_scripts': [
      'isosolve = isosolve:main_arg',
      ],
   },
   classifiers=[
      'Environment :: Console',
      'Intended Audience :: End Users/Desktop',
      'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3 :: Only',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
   ],
   project_urls={
      'Documentation': 'https://readthedocs.org/isosolve',
      'Source': 'https://github.com/sgsokol/isosolve',
      'Tracker': 'https://github.com/sgsokol/isosolve/issues',
   },
)
