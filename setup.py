from setuptools import setup

setup(
    name = 'labutils',
    version = '1.0.0',
    packages = ['sequence_analysis'],
    entry_points = {
            'console_scripts': [
            'get_plasmid_inserts = sequence_analysis.get_plasmid_inserts:main',
        ],
    },
)
