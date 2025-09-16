import ee
try:
    from .GEE_utils import initialize_ee, DRIVE_FOLDER, get_counties
except ImportError:  # allow running as a script
    from GEE_utils import initialize_ee, DRIVE_FOLDER, get_counties


initialize_ee('nlcd-469307')

print("\n⚙️  Starting county base export...")

counties = get_counties(add_area=True)

task_county_base = ee.batch.Export.table.toDrive(
    collection=counties.select(['GEOID', 'total_area_km2']),
    description='WONDER_County_Base_Data',
    folder=DRIVE_FOLDER,
    fileNamePrefix='county_base',
    fileFormat='CSV'
)
task_county_base.start()
print("✅ County base data export task started!")


