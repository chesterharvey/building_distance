# Building Distance Tools
This repo demonstrates Python tools for efficiently calculating *Building Distance*, the distance between street centerlines and the closest building. The measure is intended to provide a simple indicator of propensity for walkable streetscape design, as many walkable design qualities are theoretically promoted by proximity to buildings.

I presented the building distance concept and measurements across the 100 largest U.S. urban areas in a [paper](https://pheedloop.com/WSTLUR2021/site/sessions/?id=SESMAR7E5DQNYNCPH) at the [2021 World Symposium on Transport and Land Use Research (WSTLUR)](https://pheedloop.com/WSTLUR2021/site/home/).

## Key Contents
- **building_distance.py**: Python script containing core convenience functions for measuring building distance
- **building_distance_demo.ipynb**: Jupyter notebook providing a minimal example of measurement with demo data

## Installation
I recommend creating a [conda environment with OSMnx](https://osmnx.readthedocs.io/en/stable/), which includes almost all necessary dependencies:
```
conda config --prepend channels conda-forge
conda create -n my_env_name --strict-channel-priority osmnx
```
Activate that environment:
```
conda activate my_env_name
```
Then, install [StreetSpace](https://github.com/chesterharvey/StreetSpace), which is my own package of Python tools for measuring and analyzing streets, from GitHub using pip:
```
pip install git+https://github.com/chesterharvey/StreetSpace
```

## Data Sources
For convenience, the minimal example in the demo notebook uses data from OpenStreetMap (OSM) that are downloaded using [OSMnx](https://github.com/gboeing/osmnx).

In my paper, I used street centerline data from OSM and [building footprint data from Microsoft](https://github.com/microsoft/USBuildingFootprints). All data were stored in a PostGIS database for efficient spatial querying. Building distances were measured iteratively within 200x200 meter grid cells across each urban area. This grid size offered an efficient balance between the processing cost of calculating distances beteen all street-building node pairs, which increased exponentially with larger grid cells, and the cost of retrieving, transforming, and storing data for each iteratation, which would be economized with larger cells.
