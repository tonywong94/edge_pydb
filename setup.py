import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="edge_pydb",
    version="0.0.1",
    author="Tony Wong",
    author_email="tonywong94@gmail.com",
    description="Python based database for CARMA EDGE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tonywong94/edge_pydb",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
