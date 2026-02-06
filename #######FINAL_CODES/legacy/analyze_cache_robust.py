
import json
import os
from pathlib import Path

cache_path = Path("audit_cache.json")
if not cache_path.exists():
    print("Cache file not found.")
    exit()

try:
    with open(cache_path, 'r') as f:
        cache = json.load(f)
except Exception as e:
    print(f"Error loading cache: {e}")
    exit()

total = len(cache)
ok_count = sum(1 for v in cache.values() if v.get('status') == 'OK')
error_count = total - ok_count

print(f"--- Cache Statistics ---")
print(f"Total Entries: {total}")
print(f"Valid (OK): {ok_count}")
print(f"Invalid/Errors: {error_count}")

# Check for path consistency
if total > 0:
    first_path = next(iter(cache.keys()))
    print(f"\nExample Path in Cache: {first_path}")
    print(f"Path exists on current system: {os.path.exists(first_path)}")
    
    # Check if paths are absolute or relative
    is_absolute = os.path.isabs(first_path)
    print(f"Paths are absolute: {is_absolute}")

# Check for mtime/size consistency for some files
samples = list(cache.keys())[:5]
print("\n--- Consistency Check (Sample 5) ---")
for p in samples:
    if os.path.exists(p):
        stat = os.stat(p)
        cached = cache[p]
        mtime_match = cached['mtime'] == stat.st_mtime
        size_match = cached['size'] == stat.st_size
        print(f"File: {Path(p).name}")
        print(f"  MTime Match: {mtime_match} (Cache: {cached['mtime']}, Current: {stat.st_mtime})")
        print(f"  Size Match: {size_match} (Cache: {cached['size']}, Current: {stat.st_size})")
    else:
        print(f"File not found: {p}")
