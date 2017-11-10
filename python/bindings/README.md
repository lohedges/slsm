# Python bindings

We use [pybind11](https://github.com/pybind/pybind11) to generate a full set of
Python bindings to LibLSM and provide an extension module, pyslsm. When built,
the module will be placed in the python directory, so it can be directly
imported by the included [example scripts](../README.md). Due to inconcistencies
with the naming of Python library paths, the module will not be installed
when invoking `make install`. As such, you will need to manually install it,
or simply copy it to the folder where any Python scripts are located.

The bindings are written in a way so as to provide maximum flexibility and
ease of use. A few important details are given below:

## Mutable types

Many common Python types, such as `int`, `float`, etc., are
[immutable](https://en.wikipedia.org/wiki/Immutable_object),
meaning that they cannot be changed. To allow consistency with the LibSLSM
C++ API, where methods are frequently passed `double` types by reference
we provide a simple `MutableFloat` that can be used in Python.

For example,

```python
# Import the module.
import pyslsm

# Create a mutable float object.
mut_float = pyslsm.MutableFloat(3)

# Print the value of the float.
print(mut_float.value)
```

While the member `value` is of type `float`, which is immutable,
by wrapping it in the user defined type `MutableFloat` allows
the value to be modified in a member function.

Methods expecting arguments of type `double&` in the C++ API
should be instead passed `MutableFloat`.

## STL containers

For simplicity, we have bound a variety of STL vector containers that
help map functionality with their Python counterpart, the list. For
example, to create a vector of doubles:

```python
# Import the module.
import pyslsm

# Create a STL vector of doubles from a list.
vec = pyslsm.VectorDouble([1, 2, 3, 4, 5])

# Print the length of the vector.
print(len(vec))

# Append another list to the vector.
vec.extend(pyslsm.VectorDouble([6, 7, 8, 9, 10]))

# Print the new length of the vector.
print(len(vec))
```

A `VectorDouble` can be passed to any method expecting a `std::vector<double>&`.

A full list of the container types is given below:

```cpp
pyslsm.VectorBool == std::vector<bool>
pyslsm.VectorDouble == std::vector<double>
pyslsm.VectorInt == std::vector<int>
pyslsm.VectorUnsignedInt == std::vector<unsigned int>

pyslsm.VectorBoundaryPoint == std::vector<slsm::BoundaryPoint>
pyslsm.VectorCoord == std::vector<slsm::Coord>
pyslsm.VectorElement == std::vector<slsm::Element>
pyslsm.VectorHole == std::vector<slsm::Hole>
pyslsm.VectorNode == std::vector<slsm::Node>
```

## Callback functions

Pybind11 provides fantastic support for `std::function` making it trivial to
bind a Python function to a C++ function wrapper. For example, if we have a
C++ callback defined as:

```cpp
typedef std::function<double (const BoundaryPoint&)> SensitivityCallback
```

(This is a wrapper to a function that takes a `BoundaryPoint` object
and returns a double.)

To expose the callback to Python we wrap it in a class:

```cpp
//! Wrapper structure to expose the callback function to Python.
class Callback
{
public:
    //! Constructor.
	Callback();

	/// The callback function.
	SensitivityCallback callback;
};
```

Then, in Python we could bind a function to this C++ callback as follows:

```python
# Import the module.
import pyslsm

# A python function declaration for the callback.
def callback(point):
	# Calculate a sensitivity based on the point's coordinates.

	# Return the sensitivity.
	return sens

# Initialise a sensitivity callback object.
cb = pyslsm.Callback()

# Bind the Python callback.
cb.callback = callback

# Can now compute sensitivites for the boundary points.
# This assumes bp is of type pylsm::BoundaryPoint.

# Initialise the sensitivity object.
sens = pyslsm.Sensitivity()

# Compute the sensitivity.
# Here we are setting the zeroth sensitivity, i.e. for the objective.
bp.sensitivity[0] = sens.computeSensitivity(bp, cb.callback)
```

Note that binding Python functions to a C++ callback has a significant
computational overhead. As such, using sensitivity callback functions
written in Python can be rather inefficent, hence we recommend using the
C++ API whenever performance is a concern.
