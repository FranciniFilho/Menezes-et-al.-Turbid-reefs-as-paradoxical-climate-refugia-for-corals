
import requests

TOKEN = 'REDACTED_EARTHDATA_TOKEN'

url = 'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/AQUA_MODIS.20240101.L3m.DAY.CHL.chlor_a.4km.nc'
headers = {'Authorization': f'Bearer {TOKEN}'}

print(f"Testing NEW token against: {url}")
r = requests.get(url, headers=headers, timeout=30, allow_redirects=True)
print(f"Status Code: {r.status_code}")
print(f"Content-Type: {r.headers.get('Content-Type')}")

if r.status_code == 200:
    if 'text/html' in r.headers.get('Content-Type', ''):
        print(">>> TOKEN INVALID: Received HTML login page.")
    else:
        print(">>> TOKEN VALID: Received NetCDF data.")
elif r.status_code in [401, 403]:
    print(">>> TOKEN DENIED: Authentication error.")
else:
    print(f">>> UNEXPECTED STATUS: {r.status_code}")
