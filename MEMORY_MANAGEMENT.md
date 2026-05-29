# Memory Management in gwcatpy

## Overview

The gwcatpy library now includes **multi-layered memory management** to prevent out-of-memory (OOM) kills when processing large FITS files (gravitational wave sky maps). Memory checks occur at **multiple critical points**:

1. **Before each event** in `updateMaps()` loop - catches accumulated memory from previous events
2. **Before loading maps** in `calcAreas()` - last check before calling read_map()
3. **Inside `plotloc.read_map()`** during FITS conversion:
   - Before reading the file
   - After loading but before conversion
   - **Proactive estimation** of memory requirements before allocation
   - During the memory-intensive reordering operations
4. **After processing** each event with aggressive cleanup

This multi-layered approach ensures memory limits are enforced **before** operations that could cause memory explosions, not after the damage is done.

These features allow you to:

1. **Set memory limits** - Define maximum memory usage thresholds
2. **Abort operations automatically** - Stop processing when limits are exceeded (with early detection)
3. **Skip problematic files** - Continue processing other events even if one fails
4. **Monitor memory usage** - Track memory consumption during operations

## Quick Start

### Basic Usage with Memory Limits

```python
import gwcat_mod as gwcat

# Create catalog with memory limits
cat = gwcat.GWCat(
    fileIn='../data/events.json',
    max_memory_mb=2000,           # Process limit: 2GB
    max_memory_percent=80.0,      # System limit: 80%
    skip_on_memory_error=True,    # Skip events that exceed limits
    verbose=True
)

# Update maps - will skip events that cause memory issues
cat.updateMaps(verbose=True)
```

## Memory Limit Parameters

### GWCat Initialization Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_memory_mb` | float | None | Maximum process memory in MB before aborting operations |
| `max_memory_percent` | float | None | Maximum system memory percentage before aborting |
| `skip_on_memory_error` | bool | True | Skip events on memory error (True) or raise exception (False) |

### Example Configurations

#### 1. Process Memory Limit Only
```python
cat = gwcat.GWCat(
    fileIn='events.json',
    max_memory_mb=1500,  # 1.5 GB limit
    skip_on_memory_error=True
)
```

#### 2. System Memory Limit Only
```python
cat = gwcat.GWCat(
    fileIn='events.json',
    max_memory_percent=75.0,  # Stop if system RAM > 75%
    skip_on_memory_error=True
)
```

#### 3. Both Limits (Most Conservative)
```python
cat = gwcat.GWCat(
    fileIn='events.json',
    max_memory_mb=2000,
    max_memory_percent=80.0,
    skip_on_memory_error=True
)
```

## API Functions

### check_memory_limit()

Check if current memory usage exceeds specified limits.

```python
mem_info = gwcat.check_memory_limit(
    max_memory_mb=2000,
    max_memory_percent=80.0,
    verbose=True
)

print(f"Process memory: {mem_info['used_mb']:.2f} MB")
print(f"System memory: {mem_info['percent']:.1f}%")
print(f"Exceeds limit: {mem_info['exceeds_limit']}")
```

**Parameters:**
- `max_memory_mb` (float, optional): Maximum memory in MB
- `max_memory_percent` (float, optional): Maximum system memory percentage
- `verbose` (bool, optional): Print memory status

**Returns:**
- Dictionary with keys: `used_mb`, `percent`, `exceeds_limit`, `reason`

**Raises:**
- `MemoryLimitError`: If memory exceeds any specified limit

### memory_safe_operation()

Execute a function with memory limit checking.

```python
def my_function(arg1, arg2):
    # Your memory-intensive operation
    return result

result = gwcat.memory_safe_operation(
    my_function,
    arg1_value,
    arg2_value,
    max_memory_mb=2000,
    max_memory_percent=80.0,
    cleanup_func=my_cleanup,
    verbose=True
)
```

**Parameters:**
- `func` (callable): Function to execute
- `*args`: Positional arguments for func
- `max_memory_mb` (float, optional): Memory limit in MB
- `max_memory_percent` (float, optional): System memory limit percentage
- `cleanup_func` (callable, optional): Function to call for cleanup if limit exceeded
- `verbose` (bool, optional): Verbose output
- `**kwargs`: Keyword arguments for func

**Returns:**
- Result of func or None if memory limit exceeded

## Error Handling

### MemoryLimitError Exception

Raised when memory usage exceeds specified limits.

```python
try:
    cat = gwcat.GWCat(
        fileIn='events.json',
        max_memory_mb=1000,
        skip_on_memory_error=False  # Raise exception
    )
    cat.updateMaps()
    
except gwcat.MemoryLimitError as e:
    print(f"Memory limit exceeded: {e}")
    # Handle error - log it, notify user, etc.
```

### Skip vs. Crash Behavior

**Skip mode** (`skip_on_memory_error=True`):
- Prints warning message
- Skips the problematic event
- Continues processing other events
- Best for batch processing

**Crash mode** (`skip_on_memory_error=False`):
- Raises `MemoryLimitError`
- Stops all processing
- Allows custom error handling
- Best when you need to know about failures

## How It Works

### Memory Monitoring During calcAreas()

The `calcAreas()` function is the most memory-intensive operation because it:
1. Loads large FITS sky map files (can be >1GB)
2. Converts multi-order maps to healpix format
3. Calculates probability maps

The updated version:
1. **Checks memory before loading** - Ensures sufficient memory available
2. **Monitors during processing** - Catches memory spikes
3. **Cleans up after** - Explicitly deletes large objects and runs garbage collection
4. **Handles errors gracefully** - Skips problematic files or raises informative errors

```python
# Pseudo-code of what happens inside calcAreas()
def calcAreas(self, ev, verbose=False):
    try:
        # Check 1: Before loading map
        check_memory_limit(...)
        
        # Load potentially large FITS file
        map = plotloc.read_map(fitsFile)
        
        # Check 2: After loading map
        check_memory_limit(...)
        
        # Process the map
        totmap, a90 = plotloc.getProbMap(map, prob=0.9)
        
        # Clean up
        del map, totmap
        gc.collect()
        
    except MemoryLimitError:
        # Clean up and skip or raise
        ...
```

## Best Practices

### 1. Choose Appropriate Limits

```python
# For machines with 16GB RAM processing large files:
cat = gwcat.GWCat(
    max_memory_mb=4000,      # 4GB process limit
    max_memory_percent=80.0   # 80% system limit (~12.8GB)
)

# For machines with 8GB RAM:
cat = gwcat.GWCat(
    max_memory_mb=2000,      # 2GB process limit
    max_memory_percent=70.0   # 70% system limit (~5.6GB)
)
```

### 2. Enable Verbose Mode for Debugging

```python
cat = gwcat.GWCat(
    max_memory_mb=2000,
    verbose=True  # See which files cause issues
)
```

### 3. Process One Event at a Time for Large Files

```python
# Instead of processing all events:
# cat.updateMaps()

# Process individually:
for ev in ['GW230529_181500', 'GW150914', ...]:
    try:
        cat.updateMaps(event=ev, verbose=True)
    except gwcat.MemoryLimitError as e:
        print(f"Skipped {ev}: {e}")
```

### 4. Use Both Limits for Safety

```python
# Process limit catches runaway memory leaks
# System limit prevents system freeze
cat = gwcat.GWCat(
    max_memory_mb=3000,       # Process shouldn't exceed 3GB
    max_memory_percent=85.0   # System shouldn't exceed 85%
)
```

## Troubleshooting

### Problem: All events are being skipped

**Solution:** Limits may be too strict
```python
# Check current memory usage first
import psutil
import os

process = psutil.Process(os.getpid())
print(f"Current memory: {process.memory_info().rss / 1024**2:.2f} MB")
print(f"System memory: {psutil.virtual_memory().percent:.1f}%")

# Set limits above current usage
```

### Problem: Memory error for one specific file

**Solution:** Process it separately with higher limits or skip it
```python
# Skip the problematic file
cat = gwcat.GWCat(skip_on_memory_error=True)

# Or process with higher limit
cat_large = gwcat.GWCat(
    max_memory_mb=8000,
    skip_on_memory_error=False
)
cat_large.updateMaps(event='problematic_event')
```

### Problem: Still running out of memory

**Possible causes:**
1. Limits set too high for available RAM
2. Memory leak in processing code
3. File is genuinely too large

**Solutions:**
1. Lower the limits
2. Process files one at a time
3. Increase system RAM or use a machine with more memory
4. Consider downsampling the FITS files before processing

### Problem: Process killed during "Converting multiorder map"

This was a critical issue where the memory check happened **too late** - after the conversion had already started and memory had exploded.

**What was happening:**
```
Memory check: Process=309 MB, System=25%  ← Check passes
Loading map for GW240915_001357...
Converting multiorder map: ...            ← Memory explodes here
Killed                                    ← OOM killer strikes
```

**Why it happened:**
The multiorder→healpix conversion creates large intermediate arrays during `hp.ud_grade()` operations, causing memory to spike from 300MB to 7GB+ almost instantly.

**How it's fixed:**
Memory checks now happen **inside** `plotloc.read_map()` at multiple points:
1. Before reading the FITS file
2. After loading but before conversion
3. **Proactive estimation** of memory requirements
4. During the reordering loop (with cleanup between iterations)

**Now you'll see:**
```
Memory check: Process=309 MB, System=25%
Loading map for GW240915_001357...
Converting multiorder map: ...
Creating 12288 pixels at max nside=512 (total pixels: 3145728)
Estimated memory requirement: 240.0 MB
MEMORY LIMIT EXCEEDED: 309 MB + 240 MB > 500 MB
Skipping GW240915_001357 due to memory limit
```

The operation is **aborted before** the dangerous conversion starts, not after the OOM killer has already been triggered.

## Performance Impact

The memory checking adds minimal overhead:
- Memory checks take ~1-5ms
- Garbage collection is only called on cleanup
- No impact on processing speed when under limits

## See Also

- `example_memory_limits.py` - Working examples
- `testcat.py` - Test script using memory limits
- Python `psutil` documentation - Memory monitoring details
