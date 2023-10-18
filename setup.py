from setuptools import setup,find_packages



setup(name='VISOR',      
  version=1.2,
  description='VarIant SimulatOR',
  url='https://github.com/davidebolo1993/VISOR.git',
  requires=['python (>= 3.6)'],
  author='Davide Bolognini',
  author_email='davidebolognini7@gmail.com',
  license='LICENSE.txt',
  dependency_links=['https://github.com/rrwick/Badread/tarball/master#egg=Badread'],
  install_requires=['pyfaidx>=0.7.2.2', 'pysam>=0.22.0', 'pywgsim>=0.5.2', 'pybedtools>=0.9.1', 'mappy>=2.26', 'plotly==5.17.0', 'numpy>=1.26.1'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['VISOR=VISOR.VISOR:main']}          
)
