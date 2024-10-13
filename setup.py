from setuptools import setup, find_packages

setup(
    name='cyberlabrat',
    version='0.1.10',
    description='Contains the classes used for generix data analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Jasmine Butler & Remi Corne',
    author_email='remi.corne@gmail.com',
    url='https://github.com/jjb-hub/phd',
    packages=find_packages(),
    package_data={"cyberlabrat": ["json/*.json"]},
    install_requires=[
        "pandas==2.0.3",
        "matplotlib==3.7.3",
        "outlier_utils",
        "pingouin",
        "scipy",
        "remi_statannotations_fork",
        "setuptools",
        "statsmodels",
        "openpyxl",
        "pydantic",
        "pytest",
        "tqdm",
        "ipykernel",
        "ipywidgets",
        "IPython",
        "networkx",
        ], 
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)

