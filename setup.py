import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='DEWPython',
    version='2.0.2',
    author='Mohit Melwani Daswani and contributors',
    author_email='38257523+mmelwani@users.noreply.github.com',
    description='Python-Implemented Deep Earth Water Model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/mmelwani/DEWPython',
    packages=setuptools.find_packages(),
    download_url='https://github.com/mmelwani/DEWPython/archive/2.0.2.tar.gz',
    python_requires='>=3.8',
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
      ],
    package_data={
        'DEWPython': ['resources/*'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent"
    ],
    include_package_data=True
)
