
# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

# Set the project name
project(Hcal VERSION 2.1.0
             DESCRIPTION "Module for all things HCal."
             LANGUAGES CXX
)


if(BUILD_EVENT_ONLY)
  setup_event_library(module Hcal name Event
    linkdef include/Hcal/Event/LinkDef.h
    dependencies Recon::Event)
  return()
endif()

setup_library(module Hcal
              dependencies ROOT::Physics 
              Framework::Framework Hcal::Event Tools::Tools DetDescr::DetDescr 
)

setup_python(package_name ${PYTHON_PACKAGE_NAME}/Hcal)

setup_test(dependencies Hcal::Hcal)

setup_data(module Hcal)
