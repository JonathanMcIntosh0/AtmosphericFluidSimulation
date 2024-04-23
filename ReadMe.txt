Download the 3rd party zip file here: https://drive.google.com/file/d/1tLisdS6UnC0YngFZDfCtJ2hfo8Mnkvll/view?usp=sharing
Unzip into project directory then proceed with one of the 2 following options:
Option 1: Set a new environment variable for CMake named THIRDPARTY_DIR to the path of the 3rdparty directiory.
          I.e. when running CMake for the project you should have the environment variable: THIRDPARTY_DIR=[PATH_TO_DIR].
Option 2: Add the following line to the top of the CMakeLists.txt file:
            SET(ENV{THIRDPARTY_DIR} 3rdparty)
          Or use an absolute path like
            SET(ENV{THIRDPARTY_DIR} C:/...PATH_TO_PROJECT.../AtmosphericFluidSimulation/3rdparty)
