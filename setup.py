from setuptools import setup, find_packages

setup(
    name='franken_api',
    version='1.0.0',
    description='Curator API\'s',
    url='',
    author='',

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.10',
        
    ],

    keywords='franken, Plots, Plotly, Json',

    packages=find_packages(),

    install_requires=['flask_restx', 'Flask-SQLAlchemy', 'pandas', 'click', 'requests'],

    entry_points={
        'console_scripts': [
            'franken_api = franken_api.app:main'
        ]
    }
)
