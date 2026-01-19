---
name: marine_data_acquisition
description: Guides downloading oceanographic data from Copernicus, ERDDAP, and Bio-ORACLE. Useful when acquiring SST, chlorophyll, bathymetry, or other marine layers. Keywords: Copernicus, ERDDAP, Bio-ORACLE, DHW, SST.
---

# Marine Data Acquisition

This skill guides the autonomous acquisition and preprocessing of oceanographic data from major sources including the Copernicus Marine Service (CMEMS), NOAA CoastWatch (ERDDAP), and Bio-ORACLE. It ensures data integrity through standardized naming conventions, CRS verification, and efficient download management.

## When to use this skill

- Use this when acquiring environmental variables (e.g., SST, Chlorophyll-a, Bathymetry) for marine ecology analysis.
- Use this to automate data downloads from CMEMS or ERDDAP servers.
- Use this when you need standardized variable names across different data sources.
- Use this to calculate derived marine metrics like Degree Heating Weeks (DHW) for coral bleaching risk.

## How to use it

### Step 1: Pre-download Verification
Always check for existing local files before initiating an API call to avoid redundant downloads and respect server rate limits.

### Step 2: Data Acquisition
Use source-specific libraries (e.g., `copernicusmarine`, `erddapy`, `sdmpredictors`) with defined bounding boxes and date ranges.

### Step 3: Standardize Naming
Immediately rename acquired variables to consistent keys (e.g., `sst`, `chl_a`, `bathymetry`) to ensure downstream workflow compatibility.

### Step 4: CRS & Metadata Verification
Verify and explicitly set the Coordinate Reference System (typically WGS84 / EPSG:4326) and log download parameters for reproducibility.

### Step 5: Derived Metrics Calculation
Apply physical and ecological formulas (e.g., rolling windows for DHW) to create high-level analysis layers.

## Common pitfalls

1. **Hitting rate limits**: Implement `time.sleep()` between batch requests to avoid IP blocking from major servers.
2. **CRS mismatch**: ERDDAP and other sources may return varying projections; always reproject to a consistent CRS after download.
3. **Unit confusion**: Verify units (e.g., Kelvin vs. Celsius, mg/m³ vs. kg/m³) and convert to standard scientific units immediately.
4. **Time zone issues**: Most oceanographic sources use UTC; ensure local time offsets are handled if required for field data matching.
