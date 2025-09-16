import ee
try:
    from .GEE_utils import initialize_ee, M2_TO_KM2, NLCD_YEARS, DRIVE_FOLDER, round4, get_counties
except ImportError:
    from GEE_utils import initialize_ee, M2_TO_KM2, NLCD_YEARS, DRIVE_FOLDER, round4, get_counties


initialize_ee('nlcd-469307')

print("\n⚙️  Starting per-year NLCD export...")

# County boundaries
counties = get_counties(add_area=False).select(['GEOID'])

for year in NLCD_YEARS:
    print(f"\n🏞️ Exporting NLCD for {year}...")
    nlcd_image = ee.Image(f'USGS/NLCD_RELEASES/2019_REL/NLCD/{year}').select('landcover')

    # Classes
    forest = nlcd_image.eq(41).Or(nlcd_image.eq(42)).Or(nlcd_image.eq(43))
    water = nlcd_image.eq(11)
    urban = nlcd_image.gte(21).And(nlcd_image.lte(24))
    agriculture = nlcd_image.eq(81).Or(nlcd_image.eq(82))
    cropland = nlcd_image.eq(82)
    pasture = nlcd_image.eq(81)
    wetlands_woody = nlcd_image.eq(90)
    wetlands_herb = nlcd_image.eq(95)
    wetlands = wetlands_woody.Or(wetlands_herb)
    shrub = nlcd_image.eq(52)
    grassland = nlcd_image.eq(71)
    barren = nlcd_image.eq(31)

    area_km2 = ee.Image.pixelArea().multiply(M2_TO_KM2)

    area_stack = forest.multiply(area_km2) \
        .addBands(water.multiply(area_km2)) \
        .addBands(urban.multiply(area_km2)) \
        .addBands(agriculture.multiply(area_km2)) \
        .addBands(cropland.multiply(area_km2)) \
        .addBands(pasture.multiply(area_km2)) \
        .addBands(wetlands.multiply(area_km2)) \
        .addBands(wetlands_woody.multiply(area_km2)) \
        .addBands(wetlands_herb.multiply(area_km2)) \
        .addBands(shrub.multiply(area_km2)) \
        .addBands(grassland.multiply(area_km2)) \
        .addBands(barren.multiply(area_km2)) \
        .rename([
            'nlcd_forest_km2', 'nlcd_water_km2', 'nlcd_urban_km2', 'nlcd_agriculture_km2', 'nlcd_cropland_km2', 'nlcd_pasture_km2',
            'nlcd_wetland_km2', 'nlcd_wetland_woody_km2', 'nlcd_wetland_herb_km2', 'nlcd_shrub_km2', 'nlcd_grassland_km2', 'nlcd_barren_km2'
        ])

    stats = area_stack.reduceRegions(
        collection=counties,
        reducer=ee.Reducer.sum(),
        scale=30
    )

    def fmt(feature):
        return feature.set({
            'Year': year,
            'nlcd_forest_km2': round4(feature.get('nlcd_forest_km2')),
            'nlcd_water_km2': round4(feature.get('nlcd_water_km2')),
            'nlcd_urban_km2': round4(feature.get('nlcd_urban_km2')),
            'nlcd_agriculture_km2': round4(feature.get('nlcd_agriculture_km2')),
            'nlcd_cropland_km2': round4(feature.get('nlcd_cropland_km2')),
            'nlcd_pasture_km2': round4(feature.get('nlcd_pasture_km2')),
            'nlcd_wetland_km2': round4(feature.get('nlcd_wetland_km2')),
            'nlcd_wetland_woody_km2': round4(feature.get('nlcd_wetland_woody_km2')),
            'nlcd_wetland_herb_km2': round4(feature.get('nlcd_wetland_herb_km2')),
            'nlcd_shrub_km2': round4(feature.get('nlcd_shrub_km2')),
            'nlcd_grassland_km2': round4(feature.get('nlcd_grassland_km2')),
            'nlcd_barren_km2': round4(feature.get('nlcd_barren_km2')),
        })

    stats = stats.map(fmt)

    task = ee.batch.Export.table.toDrive(
        collection=stats,
        description=f'WONDER_NLCD_Landuse_{year}',
        folder=DRIVE_FOLDER,
        fileNamePrefix=f'nlcd_landuse_{year}',
        fileFormat='CSV',
        selectors=['GEOID', 'Year', 'nlcd_forest_km2', 'nlcd_water_km2', 'nlcd_urban_km2', 'nlcd_agriculture_km2', 'nlcd_cropland_km2', 'nlcd_pasture_km2', 'nlcd_wetland_km2', 'nlcd_wetland_woody_km2', 'nlcd_wetland_herb_km2', 'nlcd_shrub_km2', 'nlcd_grassland_km2', 'nlcd_barren_km2']
    )
    task.start()
    print(f"✅ NLCD {year} export task started!")

print("\n🎉 All per-year NLCD export tasks have been submitted (to Drive folder 'WONDER').")


