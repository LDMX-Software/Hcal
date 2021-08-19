
# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

# Set the project name
project(Hcal VERSION 2.1.0
             DESCRIPTION "Module for all things HCal."
             LANGUAGES CXX
)


# Search and configure ROOT
find_package(ROOT CONFIG REQUIRED)

option(BUILD_EVENT_ONLY "Build the event library." ON)
if(BUILD_EVENT_ONLY)
  # Generate the files needed to build the event classes.
  setup_library(module Hcal name Event
                dependencies ROOT::Core 
                             Recon::Event)
  return()
endif()

setup_library(module Hcal
              dependencies ROOT::Physics 
              Framework::Framework Recon::Event Hcal::Event Tools::Tools DetDescr::DetDescr 
)

setup_python(package_name ${PYTHON_PACKAGE_NAME}/Hcal)

setup_test(dependencies Hcal::Hcal)
