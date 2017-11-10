# Tests

A test suite is provided in the `tests` directory. Once LibSLSM has been built,
unit tests can be run as follows:

```cpp
./tests/runtests
```

The testing framework uses a modified version of [MinUnit.h](../src/MinUnit.h)
and [Debug.h](../src/Debug.h) adapted from Zed Shaw's
[Learn C The Hard Way](http://c.learncodethehardway.org/book).
