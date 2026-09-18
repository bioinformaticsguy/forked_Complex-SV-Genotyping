#ifndef READ_LENGTH_POLICY_HPP
#define READ_LENGTH_POLICY_HPP

#include <algorithm>

namespace ReadLengthPolicy
{
inline constexpr int tolerance = 1;

inline bool compatible(int first, int second)
{
    return std::max(first, second) - std::min(first, second) <= tolerance;
}
}

#endif
