// Portability shim for Traveler, installed by vadr-install.sh as
// <traveler>/src/include/bits/stdc++.h
//
// One Traveler source file, src/utils/convex_hull.cpp, includes
// <bits/stdc++.h>. That header is a GNU libstdc++ implementation detail: it
// does not exist in libc++, which is the standard library used by the clang
// that MacOS/X ships as 'g++'. Traveler's own documentation says the program
// does not build outside a GCC container for this reason.
//
// Traveler's build passes -I<traveler>/src/include before any system include
// directory, so a file placed at <traveler>/src/include/bits/stdc++.h is found
// first on every platform. It simply includes the standard headers that
// libstdc++'s own <bits/stdc++.h> makes available, which is all that
// convex_hull.cpp needs from it. The same file is installed on Linux and on
// MacOS/X so that both platforms compile Traveler from identical sources.

#ifndef VADR_BITS_STDCXX_SHIM_H
#define VADR_BITS_STDCXX_SHIM_H

#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <chrono>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <exception>
#include <forward_list>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <random>
#include <ratio>
#include <regex>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <system_error>
#include <thread>
#include <tuple>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <valarray>
#include <vector>

#endif
