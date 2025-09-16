import ee
try:
    from .GEE_utils import initialize_ee, M2_TO_KM2, YEARS, DRIVE_FOLDER, round4, get_counties
except ImportError:
    from GEE_utils import initialize_ee, M2_TO_KM2, YEARS, DRIVE_FOLDER, round4, get_counties


initialize_ee('nlcd-469307')

print("\n⚙️  Starting per-year JRC export (permanent & seasonal)...")

counties = get_counties(add_area=False).select(['GEOID'])
jrc_yearly = ee.ImageCollection('JRC/GSW1_4/YearlyHistory')

for year in YEARS:
    print(f"\n💧 Exporting JRC water for {year}...")
    img = jrc_yearly.filter(ee.Filter.eq('year', year)).first().select('waterClass')
    perm = img.eq(3).multiply(ee.Image.pixelArea()).multiply(M2_TO_KM2)
    seas = img.eq(2).multiply(ee.Image.pixelArea()).multiply(M2_TO_KM2)

    stats = perm.addBands(seas) \
        .rename(['jrc_permanent_water_km2', 'jrc_seasonal_water_km2']) \
        .reduceRegions(collection=counties, reducer=ee.Reducer.sum(), scale=30)

    def fmt(feature):
        return feature.set({
            'Year': year,
            'jrc_permanent_water_km2': round4(feature.get('jrc_permanent_water_km2')),
            'jrc_seasonal_water_km2': round4(feature.get('jrc_seasonal_water_km2')),
        })

    stats = stats.map(fmt)

    task = ee.batch.Export.table.toDrive(
        collection=stats,
        description=f'WONDER_JRC_Water_{year}',
        folder=DRIVE_FOLDER,
        fileNamePrefix=f'jrc_water_{year}',
        fileFormat='CSV',
        selectors=['GEOID', 'Year', 'jrc_permanent_water_km2', 'jrc_seasonal_water_km2']
    )
    task.start()
    print(f"✅ JRC {year} export task started!")

print("\n🎉 All per-year JRC export tasks have been submitted (to Drive folder 'WONDER').")


