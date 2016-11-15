# Tests

A test suite is provided in the `tests` directory. To run all unit tests:

```cpp
make test
```

Note that this will compile the tests (and library) using the default compilation
flags (`release`). To build the library and tests in a specific mode, run, e.g.

```cpp
make devel test
```

The testing framework uses a modified version of [MinUnit.h](../src/MinUnit.h)
and [Debug.h](../src/Debug.h) adapted from Zed Shaw's
[Learn C The Hard Way](http://c.learncodethehardway.org/book).
