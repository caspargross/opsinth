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
    install_requires=[
        'pysam',
        'edlib'
        'matplotlib',
        'mappy'
    ],
    entry_points={
        'console_scripts': [
            'opsinth=opsinth.__main__:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
) 