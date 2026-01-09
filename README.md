# Opaque/DarkMapper

**Robotic Navigation and Mapping in Confined Environments Using Proprioception**

This repository contains the experimental code and simulation framework developed for my PhD dissertation research entitled, [Robot Mapping with Proprioceptive Spatial Awareness in Confined and Sensor-Challenged Environments](doc/Everist_Dissertation_2015.pdf).

The project investigates how snake-like robotic probes can explore, map, and navigate through unknown confined environments when all external sensors have failed, relying solely on proprioception (joint angles and body contact) and physical interaction with the environment.

## Research Overview

Traditional robotic mapping relies on exteroceptive sensors (cameras, LIDAR, sonar) that can fail in harsh environments—collapsed buildings, deep caves, narrow pipes, or disaster zones. This research demonstrates that a flexible snake robot can build accurate maps using only:

- **Proprioceptive sensing**: Joint angles and actuator states
- **Tactile feedback**: Contact detection from body-environment interaction
- **Motion planning**: Adaptive locomotion gaits and exploration behaviors

The system uses a **physics-based simulator** (Bullet Physics + Ogre3D) to validate navigation algorithms that enable exploration and mapping through probabilistic pose estimation, Bayesian mapping, and particle filtering.

### Key Contributions

- Novel mapping approach using only proprioception and body contact
- Particle filter-based localization without external landmarks
- Adaptive motion primitives for snake robot navigation
- Junction detection and branching path representation
- ICP-based path alignment and loop closure

## Documentation

- **Final Dissertation**: `doc/Everist_Dissertation_2015.pdf` - Submitted and Accepted Dissertation
- **Dissertation Defense**: `doc/defense/` - PowerPoint Slide Deck
- **Dissertation Materials**: `doc/Dissertation/` - Writing Materials, incomplete
- **Related Research Papers**: `doc/*.pdf` - Published conference/journal papers
- **Algorithm Overview**: `doc/algorithm.txt` - High-level algorithm description
- **Developer Guide**: `CLAUDE.md` - Detailed codebase documentation

## Quick Start

### Prerequisites

This project was developed primarily on **Windows** with **Visual Studio C++** and has been migrated to **Python 3.x**. You'll need:

- **Python 3.x** (tested with Python 3.6+)
- **Visual Studio C++** or equivalent C++ compiler
- **Cython** for building extension modules
- **NumPy, SciPy, Matplotlib**
- **Pillow** (PIL fork for Python 3)

**Physics & Graphics Libraries**:
- **Bullet Physics** - Physics simulation engine
- **Ogre3D** - 3D visualization (optional for quick mode)
- **OpenCV** - Computer vision utilities
- **CGAL** - Computational geometry (for alpha shapes)

> **Note**: Many dependencies are included in `ciba/` directory as zip files, tarballs, or installers. This project predates modern Python package management practices and includes vendored dependencies.

### Building Extensions

Before running simulations, compile the Cython extension modules:

```bash
# Build main modules
cd src/modules
python3 setup.py build_ext --inplace

# Build individual modules (as needed)
cd src/modules/bulletprobe && python3 setup.py build_ext --inplace
cd src/modules/ogre && python3 setup.py build_ext --inplace
cd src/modules/alphashape && python3 setup.py build_ext --inplace
cd src/modules/medialaxis && python3 setup.py build_ext --inplace
cd src/modules/nelmin && python3 setup.py build_ext --inplace
cd src/modules/toro_module/trunk && python3 setup.py build_ext --inplace
```

> **Warning**: Building these modules requires the corresponding C++ libraries (Bullet, Ogre, CGAL) to be installed and accessible to the compiler.

### Running Simulations

**Full simulation with 3D visualization:**
```bash
cd src
python3 startSim.py --mapFile=testData/curveEnv/mapFile0004.txt
```

**Quick mode (no graphics, faster):**
```bash
cd src
python3 startQuick.py --mapFile=testData/curveEnv/mapFile0004.txt
```

**Batch testing across multiple environments:**
```bash
cd src/tests
python3 batchTests.py --maxNumPoses 100 --startOffset 0.0
```

### Command-Line Options

- `--mapFile <path>`: Load environment map from file
- `--maxNumPoses <int>`: Maximum number of poses to generate (default: 300)
- `--numPoseParticles <int>`: Number of pose particles for particle filter (default: 40)
- `--bloomFeature` / `--no-bloomFeature`: Enable/disable bloom spatial features
- `--bendFeature` / `--no-bendFeature`: Enable/disable bend spatial features
- `--startOffset <float>`: Probe start displacement offset
- `--hideWindow <bool>`: Hide visualization window (startSim only)

## Project Structure

```
opaque/
├── src/                    # All source code
│   ├── opaque/            # Core navigation & mapping system
│   │   ├── robot/         # Snake robot models (BulletProbe, QuickProbe)
│   │   ├── behaviors/     # Motion primitives and gaits
│   │   ├── maps/          # Bayesian mapping, particle filter, path graphs
│   │   ├── pose/          # Pose estimation from proprioception
│   │   └── tests/         # Integration tests
│   ├── modules/           # C++/Cython extension modules
│   │   ├── bulletprobe/   # Bullet physics wrapper
│   │   ├── ogre/          # Ogre3D Python bindings
│   │   ├── alphashape/    # Alpha shape computation (CGAL)
│   │   ├── medialaxis/    # Medial axis extraction
│   │   └── *.pyx          # Math utilities (transforms, ICP, stability)
│   ├── mapLibrary/        # Test environment definitions
│   ├── testData/          # Test datasets and results
│   ├── tests/             # Unit tests
│   ├── startSim.py        # Main entry point (with visualization)
│   └── startQuick.py      # Quick mode entry point (no graphics)
├── doc/                   # Dissertation, papers, presentations
├── ciba/                  # Vendored dependencies (libraries, installers)
├── CLAUDE.md              # Detailed developer documentation
└── README.md              # This file
```

## Architecture Overview

### Robot Model

The snake robot consists of **~40 segments**, each with:
- Configurable length, width, and height
- Joint angle control with torque limits
- Contact sensors for tactile feedback
- Two simulation modes:
  - **BulletProbe**: Full physics simulation (slower, more accurate)
  - **QuickProbe**: Kinematic simulation (faster, no physics)

### Mapping System

The mapping architecture uses a sophisticated probabilistic approach:

1. **Pose Estimation** (`opaque/pose/`): Extract robot configuration from joint angles and contact points
2. **Bayesian Mapping** (`opaque/maps/BayesMapper.py`): Incremental map building with uncertainty
3. **Particle Filter** (`opaque/maps/ParticleFilter.py`): Localization within the map
4. **Path Graphs** (`opaque/maps/MapState.py`): Represent environment as branching paths with junctions
5. **ICP Alignment** (`opaque/maps/gen_icp.py`): Iterative Closest Point for path matching
6. **Branch Detection**: Identify when robot departs from known paths

**Key Concept**: Poses are organized into **paths** that branch at **junctions**. The system detects new branches when the robot explores previously unmapped areas, using spline-fitted skeletons to represent paths and ICP to align overlapping segments.

### Navigation Behaviors

Hierarchical motion primitives (`opaque/behaviors/`):
- **AdaptiveStep**: Adaptive stepping for general navigation
- **PathStep**: Follow computed paths through the map
- **PokeWalls**: Tactile wall exploration
- **Concertina gaits**: Snake-like locomotion patterns
- **FrontExtend**, **HoldPosition**: Basic movement primitives

### Parallelization

The system uses Python `multiprocessing` extensively for computationally intensive operations:
- ICP computation across pose pairs
- Particle filter updates
- Branch detection and path comparison

Process pools are automatically cleaned up on exit.

## Python 3 Migration

✅ This codebase has been **fully migrated from Python 2.x to Python 3.x**:

- All files in `src/` converted using `2to3`
- All Cython modules in `src/modules/` updated with Python 3 syntax
- Setup scripts migrated: `Cython.Distutils` → `Cython.Build`
- Language level directive added: `compiler_directives={'language_level': "3"}`
- Print statements → `print()` functions
- `cPickle` → `pickle`
- Dictionary iterators updated

See `CLAUDE.md` for detailed migration notes.

## Testing

**Unit tests**:
```bash
cd src
python3 tests/testMap.py
python3 tests/testDijkstra.py
```

**Behavior tests**:
```bash
cd src/opaque/behaviors
python3 testAdaptiveCosine.py
```

**Batch testing**:
```bash
cd src/tests
python3 batchTests.py --maxNumPoses 100
```

## Known Limitations

- **Platform**: Primarily developed and tested on Windows with Visual Studio C++
- **Dependencies**: Includes vendored libraries (non-standard for modern Python projects)
- **Build complexity**: Requires manual installation of C++ dependencies (Bullet, Ogre, CGAL)
- **Documentation**: Code predates modern documentation practices; see `CLAUDE.md` for guidance

## Citation

If you use this code in your research, please cite:

```bibtex
@phdthesis{everist_phd_2015,
  author = {Jacob Everist},
  title = {Robot Mapping with Proprioceptive Spatial Awareness in Confined and Sensor-Challenged Environments},
  school = {[University of Southern California]},
  year = {2015},
  note = {Available at: https://github.com/jacobeverist/opaque/doc/Everist_Dissertation_2015.pdf}
}
```

See `doc/*.pdf` for related papers.

## License

[Specify license here - e.g., MIT, GPL, Academic Use Only, etc.]

## Contact

For questions about this research or code:
- **Author**: Jacob Everist
- **Dissertation**: See `doc/Dissertation/`
- **Papers**: See `doc/*.pdf`

---

**Note**: This is research code developed during PhD studies. It demonstrates novel algorithms for proprioceptive mapping but was not designed for production use or as a reusable library. The unconventional project structure reflects the exploratory nature of dissertation research.
