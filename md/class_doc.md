# Classes Documentation

## SimEvent class
It generates a random number of tracks from a charged particle with a 
given angle and propagates them in the Z axis direction through NPLAN 
planes. Also digitizes the particle possition in each plane.


**Location:** *simulation/simulation.py*

### GenerateEvent.digitization() method
It simulates the digital answer in NPLAN planes of detectors, in which:
- the coordinates (nx,ny) of the crossed pad are determined.
- the flight time is determined integrating tint.

## TrackFinding class
It reconstructs the tracks using Kalman Filter. It Finds tracks for hits 
on any TRASGO-like detector

**Location:** *reconstruction/track_reconstruction.py*

## RootUnpacker class

**Location:** *real_data/root_unpacker.py*

## Represent3D class

**Location:** *represent/represent_3d.py*

