import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", 'r') as req:
    requirements = req.read()


#version_num = input("Input the new version number: ")
setuptools.setup(
    name="edge_pydb",
    version="1.5.0",  #version=version_num,
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
    exclude_package_data={"": ['*/build/*']},
    install_requires=requirements,
)
