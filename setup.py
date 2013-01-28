try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Placevent solvent placement software',
    'author': 'Daniel J. Sindhikara',
    'url':'http://dansindhikara.com/Software/Entries/2012/6/22_Placevent_New.html',
    #'download_url':'Where to download it.',
    'author_email': 'sindhikara@gmail.com',
    'version': '1.2',
    'install_requires': ['numpy', 'argparse' ],
    'packages': ['placevent', 'placevent.tests'],
    'package_data': {'geometry': ['Eulersets/*'], 'pprism.tests': 
                                                ['data/AlaDP_small/*']},
    #'py_modules': ['rismgeometry'],
    'scripts': ['rismmap/rismmap.py', 'pprism/tests/pprism_tests.py', 'rismmap/tests/rismmap_test.sh'],
    'name': 'rismmap'
}
setup(**config)  
