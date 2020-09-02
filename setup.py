from setuptools import setup, find_packages

setup(
    name='franken_api',
    version='1.0.0',
    description='Json file to plot franken plots',
    url='',
    author='Vinay kaikala',

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        
    ],

    keywords='franken, Plots, Plotly, Json',

    packages=find_packages(),

    install_requires=['flask-restplus==0.12.1', 'Flask-SQLAlchemy==2.4.0', 'requests', 'click', "pytest==5.0.1", "PyMysql==0.9.3"],

    entry_points={
        'console_scripts': [
            'franken_api = franken_api.app:main'
        ]
    }
)
