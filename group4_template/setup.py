#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['pandas',
                'requests',
                'click',
                'flask',
                'numpy',
                'werkzeug']

test_requirements = ['pytest>=3', ]

setup(
    author="Group 4",
    author_email='kritiamin6461@qmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Group Project Package",
    entry_points={
        'console_scripts': [
            'group4=group4.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='group4',
    name='group4',
    packages=find_packages(include=['group4', 'group4.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
