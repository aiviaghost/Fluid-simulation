add_library (interpolation STATIC)
target_sources (
       interpolation
       PUBLIC [[interpolation.hpp]]
       PRIVATE [[interpolation.cpp]]
)
target_link_libraries (interpolation PRIVATE CG_Labs_options glm)

add_library (parametric_shapes STATIC)
target_sources (
       parametric_shapes
       PUBLIC [[parametric_shapes.hpp]]
       PRIVATE [[parametric_shapes.cpp]]
)
target_link_libraries (parametric_shapes PRIVATE bonobo CG_Labs_options)



# Assignment 5
add_executable (EDAF80_Assignment5)
target_sources (
	EDAF80_Assignment5
	PRIVATE
		[[assignment5.hpp]]
		[[assignment5.cpp]]
)
target_link_libraries (EDAF80_Assignment5 PRIVATE assignment_setup parametric_shapes)
copy_dlls (EDAF80_Assignment5 "${CMAKE_CURRENT_BINARY_DIR}")



# Fluid
add_executable (EDAF80_Fluid)
target_sources (
	EDAF80_Fluid
	PRIVATE
		[[fluid.hpp]]
		[[fluid.cpp]]
)
target_link_libraries (EDAF80_Fluid PRIVATE assignment_setup parametric_shapes)
copy_dlls (EDAF80_Fluid "${CMAKE_CURRENT_BINARY_DIR}")


install (
	TARGETS
		EDAF80_Assignment5
	DESTINATION [[bin]]
)
