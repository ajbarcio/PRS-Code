# PRS Spring Designer

This is the repository for a monolithic planar rotary spring design tool. Currently, there is a manual parameter-tuning functionality that can be loaded with a couple different spring path (spiral, polynomial) and spring thickness (constant I, piecewise quadratic I) definitions, as well as form factor constraints to determine the geometry needed to achieve a desired stiffness and torque capacity. Under work is an in-loop thickness optimizer to produce uniformly loaded springs.

## Todo to improve:
 - Automated parameter tuning to achieve desired performance
 - Hysteresis predictor/minimizer (will require empirical test data)
 - Higher variety of built-in path/thickness definitions

## Versions

The main branch is currently based off of the modules branch, which is the most up-to-date working branch.

## Usage and Structure

Todo: Make this into an installable python package

### FOR NOW:

#### Important folders:

##### modules Folder

Contains any custom modules needed for this simulator. 

CRSCDEF, PATHDEF modules contain objects for defining the path and thickness of the spring profile. Additionally objects can be added to these, which, as long as they follow the abstract base classes included in each module, will work with the rest of the simulator.

spring module contains objects defining whole springs. Spring class requries defined path AND thickness, optimized spring simply requires a defined path. This module contains all methods and attributes required to deform a spring and process results (plotting, calculating stress, etc.)

materials is a database of materials to use. Currently it is very limited.

utils contains various mathematical and graphical utilities (plotting colorful lines, numerical methods, geometry math) used in the simulator

StatProfile is a module written by Dr. Gray Thomas used for timing iterative methods.

##### tests folder

the tests folder contains tests which can be run with pytest (if you have it installed). This should not require any plugins, let me know if it does. Currently, this tests the included crsc classes, path classes, and spring deformation methods.

##### debugAndLearning folder

this test contains some debugging/test scripts that don't run when pytest is used, as well as random scripts that I use to teach myself how computers work while I write this simulator.

#### TOP DIRECTORY SCRIPTS

I will add to these as time goes on

##### manualSpringDesigner

This script can be run to design any spring whos form factor restrictions are outlined in Spring_Constraints.ods. You have to manually adjust the parameters of whatever thickness, path definitions you choose to use.