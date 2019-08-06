import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirement.txt", 'r') as req:
    requirement = req.read()

setuptools.setup(
    name="edge_pydb",
    version="0.1.20",
    author="Tony Wong",
    author_email="tonywong94@gmail.com",
    description="Python based database for CARMA EDGE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tonywong94/edge_pydb",
    packages=['edge_pydb'],
    package_data={
        'edge_pydb': [
            'dat_*/*/*.csv',
            'img_*/*/*.csv',
            'img_*/*.hdf5',
            '_config.json'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license='MIT',
    include_package_data=True,
    install_requires=requirement,
)
