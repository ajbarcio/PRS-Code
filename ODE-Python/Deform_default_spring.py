# import importlib.util
import sys
import argparse

# from Abaqus_Interface import Abaqus_Spring
import spring

# Create argument parser
parser = argparse.ArgumentParser(description='Example script with command line arguments')

# Add arguments
parser.add_argument('-n', '--n_option', type=str, help='name for spring')
# parser.add_argument('-b', '--b_option', type=str, help='Description of -b option')

# Parse arguments
args = parser.parse_args()

defaultSpring = spring.generate_default_spring(name="default_spring")
defaultSpringAbaqus = Abaqus_Spring("default_spring")

defaultSpringAbaqus.deform_with_abaqus()