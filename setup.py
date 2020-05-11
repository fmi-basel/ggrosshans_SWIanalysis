from setuptools import setup, find_packages

contrib = [
    'Yannick Hauser', 'Markus Rempfler', 
]

# setup.
setup(
    name='Streamlit_GUI_SWI_final',
    version='0.1',
    description='SWI analyzer',
    author=', '.join(contrib),
    license='',
    packages=find_packages(exclude=[
        'tests',
    ]),
    install_requires=[
        'pandas','numpy', 'seaborn', 'matplotlib', 'scipy', 'uncertainties', 'streamlit', 
        'time', 'bokeh', 'PIL', 'procswi', 'plotting_SWI', 'pytest', 'io', 'os'
    ],
    # I believe zipping makes the import rather slow.
    zip_safe=False)