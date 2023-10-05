import setuptools
setuptools.setup(
    name="table1",
    version="1.1.0",
    author="Hyung-Chul Lee",
    author_email="vital@snu.ac.kr",
    description="Python Libray for Table 1",
    long_description="Python Libray for Table 1",
    long_description_content_type="text/markdown",
    url="https://github.com/vitaldb/table1",
    install_requires=['numpy', 'pandas', 'rpy2', 'scipy'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)