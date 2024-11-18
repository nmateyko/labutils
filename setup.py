from setuptools import setup

setup(
    name="labutils",
    version="1.1.0",
    packages=['sequence_analysis'],
    install_requires=['edlib'],
    entry_points={
        "console_scripts": [
            "get_plasmid_inserts=sequence_analysis.get_plasmid_inserts:main",
        ],
    },
    test_suite='tests',
    author="Nick Mateyko",
    author_email="nick.mateyko@gmail.com",
    description="Various scripts, mostly for sequence analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nmateyko/labutils",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
