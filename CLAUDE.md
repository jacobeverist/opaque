# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is **Opaque/DarkMapper**, a robotic navigation and mapping system for snake-like probe robots operating in unknown environments. The system uses a physics-based simulator (Bullet/Ogre) to test navigation algorithms that enable a flexible snake robot to explore, map, and navigate through complex environments using tactile sensing.

## Running the Simulator

### Main Entry Points

**Full Simulation (with visualization):**
```bash
cd src
python startSim.py
```

**Quick Mode (no physics visualization):**
```bash
cd src
python startQuick.py
```

### Common Command-Line Arguments

Both `startSim.py` and `startQuick.py` accept:
- `--mapFile <path>`: Load environment map from file (e.g., `--mapFile=testData/curveEnv/mapFile0004.txt`)
- `--maxNumPoses <int>`: Maximum number of poses to generate (default: 300 for startSim, varies for startQuick)
- `--numPoseParticles <int>`: Number of pose particles for particle filter (default: 40)
- `--bloomFeature` / `--no-bloomFeature`: Enable/disable bloom spatial features (default: enabled)
- `--bendFeature` / `--no-bendFeature`: Enable/disable bend spatial features (default: disabled)
- `--startOffset <float>`: Probe start displacement offset (default: 0.0)
- `--restoreConfig <bool>`: Restore Ogre config from file (startSim only)
- `--hideWindow <bool>`: Hide Ogre visualization window (startSim only)

### Batch Testing

Run multiple test scenarios:
```bash
cd src/tests
python batchTests.py --maxNumPoses 100 --startOffset 0.0 [--bendFeature]
```

This runs tests across predefined map sets and stores results in separate directories.

## Architecture

### Core Components

**Robot Models** (`src/opaque/robot/`):
- `BulletProbe.py`: Full physics simulation using Bullet engine
- `QuickProbe.py`: Fast kinematic simulation without full physics
- `SnakeProbe.py`: Base snake robot interface
- Each probe has ~40 segments with configurable length, width, and torque limits

**Control & Navigation** (`src/opaque/`):
- `TestNavigation.py`: Main navigation controller with exploration and path-following
- `TestProcess.py`: Alternative process-based control system
- `SnakeControl.py`: Base control interface

**Behaviors** (`src/opaque/behaviors/`):
Hierarchical motion primitives:
- `AdaptiveStep.py`: Adaptive stepping behavior for navigation
- `PathStep.py`: Follow computed paths
- `PokeWalls.py`: Tactile wall exploration
- `FrontExtend.py`, `HoldPosition.py`: Basic motion primitives
- `AnchorTransition.py`, `HoldTransition.py`: Behavior transitions
- Concertina gaits and curve-fitting behaviors for complex locomotion

**Mapping System** (`src/opaque/maps/`):

The mapping architecture uses a sophisticated probabilistic approach:

- `MapUser.py`: High-level mapping interface connecting robot to map algorithms
- `BayesMapper.py`: Core Bayesian mapping algorithm
- `MapState.py`: Global map state management with path graphs and junction detection
- `LocalNode.py`: Individual pose nodes with sensor data
- `Splices.py`: Path splicing and overlap detection
- `ParticleFilter.py`: Particle filter for pose estimation and localization
- `shoots.py`: Skeleton-based path representation and branching
- `gen_icp.py`: ICP (Iterative Closest Point) for pose alignment

**Key Mapping Concepts:**
- Poses are organized into paths that branch at junctions
- New branches are detected when the robot departs from known paths
- Paths are represented as spline-fitted skeletons with control points
- Multiprocessing pools handle computationally intensive operations (ICP, particle filtering, branch computation)

**Pose Estimation** (`src/opaque/pose/`):
- `AverageContacts.py`: Estimates robot pose from tactile contact sensors

**Visualization** (`src/opaque/`):
- `DrawThings.py`: Ogre-based 3D visualization
- `CmdDrawThings.py`: Command-line visualization fallback

### C++/Cython Extensions

**Physics Integration** (`src/modules/bulletprobe/`):
- Bullet physics wrapper for snake robot simulation
- Build with: `python setup.py build_ext --inplace`

**Math Extensions** (`src/modules/`):
- `transform.pyx`: Joint transformation computations
- `stability.pyx`: Stability analysis
- `reference.pyx`: Reference frame computations
- `servo.pyx`, `func.pyx`, `icp.pyx`: Various math utilities
- Build with: `cd src/modules && python setup.py build_ext --inplace`

**Geometric Libraries**:
- `alphashape/`: Alpha shape computation (CGAL-based)
- `medialaxis/`: Medial axis extraction
- `ogre/`: Ogre3D Python bindings

### Map File Format

Map files define wall geometries as lists of polylines:
```python
walls = [
    [[x1, y1], [x2, y2], ...],  # Wall 1
    [[x3, y3], [x4, y4], ...],  # Wall 2
]
```

Test maps are in `src/mapLibrary/` and `src/testData/curveEnv/`.

## Development Workflow

### Building Extensions

Before running simulations, build Cython extensions with Python 3:

**Main modules:**
```bash
cd src/modules
python3 setup.py build_ext --inplace
```

**Individual modules** (build each as needed):
```bash
cd src/modules/bulletprobe && python3 setup.py build_ext --inplace
cd src/modules/ogre && python3 setup.py build_ext --inplace
cd src/modules/alphashape && python3 setup.py build_ext --inplace
cd src/modules/medialaxis && python3 setup.py build_ext --inplace
cd src/modules/nelmin && python3 setup.py build_ext --inplace
cd src/modules/toro_module/trunk && python3 setup.py build_ext --inplace
```

**Note**: You'll need the appropriate C++ libraries installed (Bullet, Ogre3D, OpenCV, CGAL) for the modules to compile successfully.

### Python Environment

This codebase has been **fully migrated to Python 3.x**. Key dependencies:
- Bullet physics
- Ogre3D
- Cython
- NumPy, SciPy
- Matplotlib (pylab)
- Pillow (PIL fork for Python 3)
- CGAL (for alpha shapes)
- Multiprocessing for parallelization

**Migration Status:**
- ✅ **COMPLETE**: All Python files in `src/` converted using 2to3
- ✅ **COMPLETE**: All `src/modules/` files converted to Python 3
  - ✅ Updated all 7 setup.py files: `Cython.Distutils` → `Cython.Build`
  - ✅ Added `language_level="3"` directive to all Cython modules
  - ✅ Fixed all .pyx files (12 total):
    - transform.pyx: Fixed old Cython `for i from ... >= i > ...` syntax → `for i in range(...)`
    - ogre/ogreprobe.pyx: Fixed `print` statements → `print()` functions
    - All other .pyx files were already Python 3 compatible
  - C++ wrapper files unchanged (C++ is version-agnostic)
- When converting Python 2 to 3, key changes include:
  - `print` statements → `print()` function calls
  - `import cPickle as pickle` → `import pickle`
  - Relative imports → Absolute imports with dots (e.g., `from .module import`)
  - `dict.iteritems()` → `dict.items()`
  - `dict.iterkeys()` → `dict.keys()`
  - `dict.itervalues()` → `dict.values()`
  - `range()` returns iterator, use `list(range())` if list needed
  - Division: `/` is float division, use `//` for integer division

### System Path Setup

Both entry point scripts automatically configure Python paths:
```python
sys.path.insert(1, relPath + "/modules/")
sys.path.insert(1, relPath + "/modules/nelmin")
sys.path.insert(1, relPath + "/modules/bulletprobe")
sys.path.insert(1, relPath + "/modules/medialaxis")
sys.path.insert(1, relPath + "/modules/alphashape")
sys.path.insert(1, relPath + "/modules/ogre")
```

### Multiprocessing Cleanup

The system uses multiprocessing pools extensively. Both entry points register cleanup handlers:
- `gen_icp.overlapPool`: ICP computation workers
- `MapProcess.pool_*`: Various map processing pools
- `shoots.pool_branch`: Branch computation workers
- `ParticleFilter.pool_*`: Particle filter workers

Process pools are automatically terminated on exit via `atexit.register(cleanup)`.

## Key Algorithms

### Navigation Loop

1. **Pose Estimation**: Extract robot pose from joint angles and contact sensors
2. **Map Integration**: Add new pose to map, detect overlaps with existing paths
3. **Branch Detection**: Determine if robot departed from known paths (requires sufficient distance, unique angle, featured curve)
4. **Path Planning**: Compute navigation path through map graph
5. **Behavior Execution**: Execute motion primitives to follow path
6. **Repeat**: Continue until target reached or max poses exhausted

### Path Comparison & Merging

Sibling paths (sharing same parent) are compared via:
1. Stitch parent path with section between child junctions excised
2. ICP alignment between stitched and original parent
3. Merge if angular offset < 0.5 rad and junction angle difference < π/4

Parent-child path comparisons use splice-based ICP fitting.

### Recursion Limit

The code sets `sys.setrecursionlimit(10000)` due to deep recursion in mapping algorithms.

## Testing

Individual unit tests are in:
- `src/tests/`: Basic component tests (curves, maps, Dijkstra, etc.)
- `src/opaque/tests/`: Integration tests
- `src/opaque/behaviors/testAdaptiveCosine.py`: Behavior-specific tests

Run individual tests directly:
```bash
cd src
python tests/testMap.py
```

## Documentation

High-level algorithm description: `doc/algorithm.txt`
Research papers: `doc/*.pdf`
Dissertation materials: `doc/Dissertation/`

## Important Notes

- The camera follows segment 9 of the robot during visualization
- Default test environment uses junction maps from `testData/curveEnv/`
- MAX_NODES in TestNavigation is set to 60 (configurable: 60/100/300)
- Default recursion limit is 10000 for deep mapping computations
- Random seed is fixed to 0 for reproducibility in both main scripts
- All Python code (including `src/` and `src/modules`) has been fully migrated to Python 3.x
- This project is the experimental code used in my PhD dissertation.   The goal of this project is to map confined environments when all sensors have failed.  Instead it uses proprioception and movement of the robot's body to explore, map, and navigate the environment.  This code was written before I understood how structure a python project and manage dependencies.  I put many of the libraries within the repo as zip files, tarballs, or auto-installers.  Most of the development was done in Windows so a lot of the build instructions expect that Visual Studio C++ is available.  The docs/ directory contains the dissertation, research papers, and presentations as part of my PhD studies.  The ciba/ directory contains some of the downloadable libraries that are dependencies.  The src/ directory contains all code with the `src/opaque` directory being the core software for running simulations and producing maps as results.