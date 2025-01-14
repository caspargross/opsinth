from setuptools import setup, find_packages

def readme():
    with open('README.md', 'r') as f:
        return f.read()

setup(
    name='opsinth',
    version='0.1',
    description='Opsin Analysis CLI',
    author='Caspar Gross',
    author_email='post@caspar.bio',
    url='https://github.com/caspargross/opsinth',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'opsinth': [
            'resources/*.bed',
            'resources/*.fa',
            'resources/*.fa.gz',
            'resources/*.fa.gz.fai',
            'resources/*.fa.gz.gzi',
            'resources/*.yaml',
            'resources/*.xml',
        ],
    },
    install_requires=[
        'pysam',
        'edlib',
        'matplotlib',
        'mappy',
        'pyyaml'
    ],
    entry_points={
        'console_scripts': [
            'opsinth=opsinth.__main__:main',
        ],
    },
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
) 