import shapely as sh
from shapely.geometry import Point, LineString, Polygon
import pandas as pd
import geopandas as gpd
import numpy as np
import streetspace as sp

# Extract coordinates from geometries
def gdf_to_coordinates(gdf, key=None):
    geometries = gdf.geometry.tolist()
    if isinstance(geometries[0], sh.geometry.Point) or isinstance(geometries[0], sh.geometry.LineString):
        coordinate_sets = [geometries[i].coords.xy for i in range(len(geometries))]
        x = [x for xs,ys in coordinate_sets for x in xs]
        y = [y for xs,ys in coordinate_sets for y in ys]
        df = pd.DataFrame({'x':x, 'y':y})
        if key:
            keys = gdf[key].tolist()
            keys = [k for key, coordinate_set in zip(keys, coordinate_sets) for k in [key] * len(coordinate_set[0])]
            df[key] = keys
    # Polygon is unique because exterior method is required to access coordinates and last coordinte is redundant
    elif isinstance(geometries[0], sh.geometry.Polygon):
        coordinate_sets = [geometries[i].exterior.coords.xy for i in range(len(geometries))]
        # Exclude the last coordinates for polygons because they're redundant
        x = [x for xs,ys in coordinate_sets for x in xs[:-1]]
        y = [y for xs,ys in coordinate_sets for y in ys[:-1]]
        df = pd.DataFrame({'x':x, 'y':y})
        if key:
            keys = gdf[key].tolist()
            keys = [k for key, coordinate_set in zip(keys, coordinate_sets) for k in [key] * (len(coordinate_set[0])-1)]
            df[key] = keys
    return df
    #### TODO: add options for multipoint, multilinestring, and multipolygon
    
def construct_building_segments(building_coords, key):
    # Shift coordinates up one space to build pairs representing segment ends
    building_coords_grouped = building_coords.groupby(key)
    shifted = building_coords_grouped.shift(-1)
    # Bring coordiantes that fell off the bottom back to the top
    shifted = shifted.set_index(building_coords[key])
    shifted.x = shifted.x.fillna(building_coords_grouped.x.first())
    shifted.y = shifted.y.fillna(building_coords_grouped.y.first())
    shifted = shifted.set_index(building_coords.index).rename(columns={'x':'xb', 'y':'yb'})
    building_coords = building_coords.rename(columns={'x':'xa','y':'ya'}).copy()
    return pd.concat([building_coords, shifted], axis=1)
     
# Create Cartesian product of street and building points
def df_cartesian_product(df_a, df_b, label_a='', label_b='', label_shared_only=False):
    df_a = df_a.copy()
    df_b = df_b.copy()
    # Find shared column names and add labels
    if label_shared_only:
        shared_columns = [x for x in df_a.columns if x in df_b.columns]
        if len(shared_columns) > 0:
            df_a.columns = [f'{label_a}_{x}' if x in shared_columns else x for x in df_a.columns]
            df_b.columns = [f'{label_b}_{x}' if x in shared_columns else x for x in df_b.columns]
    else:
        df_a.columns = [f'{label_a}_{x}' for x in df_a.columns]
        df_b.columns = [f'{label_b}_{x}' for x in df_b.columns]
    # Construct Cartesian product of the dataframes
    df_a['_'] = 0
    df_b['_'] = 0
    output = df_a.merge(df_b, on='_')
    output = output.drop(columns=['_'])
    return output

def calculate_street_coordinates_and_angles(street_segments, key=None, points=False):
    # Get segment start and end coordinate and calculate angles between them
    u = gdf_to_coordinates(street_segments.reset_index(), 'index').groupby('index').agg('first')
    v = gdf_to_coordinates(street_segments.reset_index(), 'index').groupby('index').agg('last')
    angles = np.degrees(np.arctan2(v.y - u.y, v.x - u.x))

    # Compile street points dataframe with the start coordinates, angles, and key
    street_points = u    
    street_points['angle'] = angles
    if key:
        street_points[key] = street_segments[key]
    if points:
        street_points = gpd.GeoDataFrame(
            street_points, 
            geometry=gpd.points_from_xy(street_points.x, street_points.y),
            crs=street_segments.crs)
    return street_points 
    
def calculate_building_distances_and_angles(cartesian_product):
    # Calculate distances and angles between street and building points
    cartesian_product['building_distance'] = np.linalg.norm(
        (cartesian_product[['street_x','street_y']].values -
         cartesian_product[['building_x','building_y']].values), 
        axis=1)
    cartesian_product['building_angle'] = np.degrees(np.arctan2(
        cartesian_product.building_y - cartesian_product.street_y, 
        cartesian_product.building_x - cartesian_product.street_x))
    return cartesian_product

def calculate_building_sides(cartesian_product):
    # Determine which side of the street each building is on
    cartesian_product['angle_difference'] = cartesian_product.street_angle - cartesian_product.building_angle
    cartesian_product['angle_difference'] = sp.normalize_azimuth_array(cartesian_product['angle_difference'], zero_center=True)
    cartesian_product['building_side'] = np.where(cartesian_product['angle_difference'] < 0, 0, 1)
    return cartesian_product

def calculate_building_lat_lon_distances(cartesian_product):
    cartesian_product['building_lat_dist'] = np.absolute(
        np.sin(cartesian_product.angle_difference * np.pi / 180) * cartesian_product.building_distance)
    cartesian_product['building_lon_dist'] = np.absolute(
        np.cos(cartesian_product.angle_difference * np.pi / 180) * cartesian_product.building_distance)
    return cartesian_product

def gdf_points_along_lines_vectorized(gdf, spacing=None, random_density=None, 
    centered=False, include_start_point=False, include_end_point=True,
    return_advanced_fields=False, return_point_geometries=False):
    '''Constructs points along linestrings in a geodataframe using an efficient, vectorized approach.  
    
    spacing: constant interval at which points will be spaced along each line
    random_density: point density (points/linear unit) at which points will be uniformly randomly distributed along each line 
    centered: If True, a space between points is centered on the linestring (nearest points are 1/2 space away)
    include_start_point: If True, the start of each linestrings is included in output
    include_end_point: If True, the end of each linestrings is included in output
    return_advanced_fields: If True, data about point positions about linestring segments are included in output
    return_point_geometries: If True, point coordinates are converted to Shapely points are returned in a GeoDataFrame    
    ''' 
    ## Convert linestrings to coordinate-based segments and calculate lengths
    # Extract linestring coordinates from geodataframe
    linestring_coords = gdf_to_coordinates(gdf.reset_index(), key='index')
    # Construct pairs of coordinates to represent segments
    linestring_coords = pd.concat([
        linestring_coords.rename(columns={'x':'xa','y':'ya'}), 
        linestring_coords.groupby('index').shift(-1).rename(columns={'x':'xb','y':'yb'})], axis=1)
    # Calculate lengths of segments
    linestring_coords['seg_length'] = np.linalg.norm(
        (linestring_coords[['xa','ya']].values -
         linestring_coords[['xb','yb']].values), 
        axis=1)
    # Calculate angle of each segment
    linestring_coords['seg_angle'] = np.degrees(np.arctan2(
        linestring_coords.yb - linestring_coords.ya, 
        linestring_coords.xb - linestring_coords.xa))
    # Calculate start and end linear references for each segment
    linestring_coords_groups = linestring_coords.groupby('index')
    linestring_coords['end'] = linestring_coords_groups.seg_length.cumsum()
    linestring_coords['start'] = linestring_coords.end.shift(1).fillna(0)
    # Drop last records for each line, which have null values becauase they don't have an end point
    linestring_coords = linestring_coords.dropna()
    # Add up total line lengths as a framework for constructing linear references
    linestring_lengths = linestring_coords_groups.seg_length.sum()

    ## Construct linear references for points along linestrings
    # Function to construct equally spaced linear references as the sum of spacing and sequential row indices
    def create_spaced_linrefs():
        return (linestring_lengths.groupby('index').cumcount() * spacing).rename('lin_ref')
    # Function to center linear references on each linestring
    def center_spaced_linrefs():
        return lin_refs + (linestring_lengths % spacing / 2) - (spacing / 2)
    # Function to pull any linear references beyond each linestring's bounds back to its ends
    def constrain_linref_bounds():
        # Replace negative linear references with 0 and those greater than the linestring length with that length
        return lin_refs.mask(lin_refs < 0).fillna(0).mask(lin_refs > linestring_lengths).fillna(linestring_lengths)
    # Construct linear references based on even spacing or stratified random spacing
    if spacing:
        # Determine the final number of points (+ 2 for the ends)
        points_n = linestring_lengths // spacing + 2
        # Construct a seperate row for each point
        linestring_lengths = linestring_lengths.loc[linestring_lengths.index.repeat(points_n)]
        # Construct stratified random points
        if random_density:
            # Construct equally spaced points 
            lin_refs = create_spaced_linrefs()
            # Center if necessary
            if centered:
                lin_refs = center_spaced_linrefs()
                lin_refs = constrain_linref_bounds()
            # Calculate ends of intervals for drawing random numbers
            ends = lin_refs.groupby('index').shift(-1)
            # Combine intervals and drop those without an end
            lin_refs = pd.DataFrame({'a':lin_refs, 'b':ends}).dropna()
            # Draw random values between the intervals
            lin_refs = pd.Series(np.random.uniform(lin_refs.a, lin_refs.b), index=lin_refs.index, name='lin_ref')
        # Construct evenly spaced points
        else:
            lin_refs = create_spaced_linrefs()
            if centered:
                lin_refs = center_spaced_linrefs()
            lin_refs = constrain_linref_bounds()
    # Construct linear references that are random along entire linestring lengths
    elif random_density:
        # Determine total number of points
        points_n = linestring_lengths // (1 / random_density) + 2
        # Construct a seperate row for each point
        linestring_lengths = linestring_lengths.loc[linestring_lengths.index.repeat(points_n)]
        # Construct linear references as random distances within uniform distributions across the whole linestring lengths
        lin_refs = pd.Series(np.random.uniform(0, linestring_lengths, len(linestring_lengths)), index=linestring_lengths.index, name='lin_ref')
        # Set the first and last linear references to 0 and the linestring length to represent ends
        unit_vector = pd.Series([1 for x in range(len(lin_refs))], index=lin_refs.index)
        lin_refs = (lin_refs * unit_vector.groupby('index').shift(1)).fillna(0)
        lin_refs = (lin_refs * unit_vector.groupby('index').shift(-1)).fillna(linestring_lengths)
    # Remove start and end points
    if not include_start_point:
        filter_vector = pd.Series([True for x in range(len(lin_refs))], index=lin_refs.index)
        lin_refs = lin_refs[filter_vector.groupby('index').shift(1).fillna(False)]
    if not include_end_point:
        filter_vector = pd.Series([True for x in range(len(lin_refs))], index=lin_refs.index)
        lin_refs = lin_refs[filter_vector.groupby('index').shift(-1).fillna(False)]

    ## Attach linear references to segments and calculate cartesian positions     
    # Convert linear references to a dataframe
    lin_refs = pd.DataFrame({'lin_ref':lin_refs}).reset_index()
    # Join linear references with linestring segments based on linear position
    lin_refs = pd.merge_asof(
        lin_refs.sort_values('lin_ref'),
        linestring_coords.sort_values('start'),
        by='index',
        left_on='lin_ref',
        right_on='start',
        direction='backward').sort_values(['index','lin_ref']).reset_index(drop=True)
    # Calculate coordinates of linear reference along the segment
    lin_refs['seg_ref'] = lin_refs.lin_ref - lin_refs.start
    dx = lin_refs.xb - lin_refs.xa
    dy = lin_refs.yb - lin_refs.ya
    a = lin_refs.seg_ref/lin_refs.seg_length
    lin_refs['x'] = lin_refs.xa + (dx * a)
    lin_refs['y'] = lin_refs.ya + (dy * a)
    
    ## Organize up outputs
    basic = OrderedDict([
        ('index','linestring_index'),
        ('x','point_x'),
        ('y','point_y'),
        ('lin_ref','point_lin_ref'),
        ('seg_angle','linestring_azimuth_at_point'),      
    ])
    advanced = OrderedDict([
        ('xa','segment_start_x'),
        ('ya','segment_start_y'),
        ('xb','segment_end_x'),
        ('yb','segment_end_y'),
        ('seg_length','segment_lenth'),
        ('start','segment_lin_ref_start'),
        ('end','segment_lin_ref_end'),
        ('seg_ref','point_lin_ref_along_segment'),
    ])
    if return_advanced_fields:
        basic.update(advanced)
    lin_refs = lin_refs[[x for x in basic.keys()]].rename(columns=basic)
    # Add geometries if specified
    if return_point_geometries:
        lin_refs['geometry'] = [sh.geometry.Point(x,y) for x, y in zip(lin_refs.point_x, lin_refs.point_y)]
        lin_refs = gpd.GeoDataFrame(lin_refs, geometry='geometry', crs=gdf.crs)
    return lin_refs

def calculate_building_distances(streets, street_id, buildings, building_id, interval=10):
    # Set up id fields
    streets = streets[[street_id, 'geometry']]
    buildings = buildings[[building_id, 'geometry']]
    # Break streets into segments and convert to points
    street_segments = sp.gdf_split_lines(streets, interval)
    street_coords = calculate_street_coordinates_and_angles(street_segments, key=street_id)
    # Convert buildings to coordinates
    building_coords = gdf_to_coordinates(buildings, key=building_id)
    # Construct building segments
    building_segment_coords = construct_building_segments(building_coords, key=building_id)
    # Construct cartisian product of street points and building segments
    cp_coords = df_cartesian_product(street_coords, building_segment_coords, label_a='street', label_b='building')
    # Calculate closest points along building segments to each street segment
    cp_coords['building_x'], cp_coords['building_y'] = sp.closest_point_along_line_vectorized(
        (cp_coords.street_x, cp_coords.street_y),
        (cp_coords.building_xa, cp_coords.building_ya), 
        (cp_coords.building_xb, cp_coords.building_yb))
    # Calculate distances and angles to buildings
    cp_coords = calculate_building_distances_and_angles(cp_coords)
    # Calculate which side of each segment buildings are on
    cp_coords = calculate_building_sides(cp_coords)
    # Calculate latitudinal and longitudinal distances to buildings from each street point
    cp_coords = calculate_building_lat_lon_distances(cp_coords)
    # Get minimum building distance on each side
    cp_coords = cp_coords.loc[cp_coords.groupby(['street_x', 'street_y', 'building_side']).building_distance.idxmin()]
    return cp_coords


    