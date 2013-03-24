try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Placevent solvent placement software',
    'author': 'Daniel J. Sindhikara',
    'url':'http://dansindhikara.com/Software/Entries/2012/6/22_Placevent_New.html',
    'download_url': 'http://dansindhikara.com/Software/Entries/2012/6/22_Placevent_New.html',
    'author_email': 'sindhikara@gmail.com',
    'version': '1.3',
    'install_requires': ['numpy'],
    'packages': ['placevent', 'placevent.tests', 'pgrid', 'ppdb'],
    'package_data': {'placevent.tests': ['1L2Y/diffs','1L2Y/*.*', '1L2Y/ref/*'],
                     'pgrid': ['shells.json']},
    'license': 'GPL',
    'scripts': ['placevent/placevent.py', 'placevent/tests/1L2Y/placevent_test.sh'],
    'name': 'placevent'
}
setup(**config)  
