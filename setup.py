import os
import os.path
from setuptools import setup

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

version = '1.0.1'
# For travis/pypi-test
if os.environ.get('TRAVIS') == 'true' and os.environ.get('TRAVIS_TAG') == '':
      N, M = os.environ['TRAVIS_JOB_NUMBER'].split('.')
      version = "{v}-a{N}.dev{M}".format(v=version, N=N, M=M)

setup(name='codonharmonizer',
      version=version,
      description='A different take on codon optimization',
      url='http://github.com/smsaladi/codonpair',
      author='Shyam Saladi',
      author_email='saladi@caltech.edu',
      license='GPLv3',
      install_requires=['pandas', 'biopython'],
      py_modules=['codonharmonizer'],
      entry_points={'console_scripts': ['codonharmonizer=codonharmonizer:main']},
      long_description=long_description,
      long_description_content_type='text/markdown',
      tests_require = [
            'pytest',
      ],
      test_suite="pytest",
      zip_safe=True,
      include_package_data=True
)
