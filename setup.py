from setuptools import setup, find_packages
def readme():
    with open('README.md') as f:
        return f.read()
long_description=''

setup(name='simphycompress',\
    version='1.0.1',\
    description='',\
    long_description=long_description,\
    url='https://github.com/merlyescalona/simphycompress',\
    download_url='https://github.com/merlyescalona/simphycompress/blob/master/dist/simphycompress.1.0.1.tar.gz',\
    author='Merly Escalona',\
    author_email='merlyescalona@uvigo.es',\
    license='GNU/GPL v3',\
    packages=['simphycompress'],\
    package_dir={'simphycompress': 'simphycompress'},\
    py_modules=[\
        'simphycompress.loggingformatter', \
        'simphycompress.msatools', \
        'simphycompress.simphycompress' \
    ],\
    install_requires=[\
        'argparse',\
        'ConfigParser',\
        'datetime',\
        'logging',\
        'numpy',\
        'setuptools',\
        'scipy'
    ],\
    classifiers=[\
        'Development Status :: 4 - Beta',\
        'Intended Audience :: Education',\
        'Intended Audience :: Science/Research',\
        'Intended Audience :: Developers',\
        'Topic :: Scientific/Engineering :: Bio-Informatics',\
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',\
        'Programming Language :: Python :: 2.7'\
    ],\
    keywords='biology phylogenomics next-generation-sequencing reference-sequences targeted-sequencing',\
    python_requires='~=2.7',\
    scripts=['bin/simphycompress'],\
    entry_points={
        'console_scripts':[\
            'simphycompress=simphycompress.__main__:main'\
        ]\
    },\
    zip_safe=False\
  )
