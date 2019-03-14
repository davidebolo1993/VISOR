from setuptools import setup,find_packages



setup(name='VISOR',      
  version=1.0,
  description='VarIant SimulatOR',
  url='https://github.com/davidebolo1993/VISOR.git',
  requires=['python (>= 3.6)'],
  author='Davide Bolognini',
  author_email='davidebolognini7@gmail.com',
  license='LICENSE.txt',
  install_requires=['pyfaidx >= 0.5.5.2', 'pybedtools >= 0.8.0'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['VISOR=VISOR.VISOR:main']}          
)
