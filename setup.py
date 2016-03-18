from setuptools import setup

s_args = {
    'name': 'eqclustering',
    'version': '0.1.0',
    'description': 'Statistical earthquake clustering algorithm',
    'author': 'Mark Williams',
    'maintainer': 'Nevada Seismological Laboratory',
    'maintainer_email': 'nvseismolab@gmail.com',
    'url': 'https//github.com/NVSeismoLab/eqclustering',
    'py_modules': ['eqclustering'],
    'install_requires': [
        'numpy',
    ],
}

# Go
setup(**s_args)

