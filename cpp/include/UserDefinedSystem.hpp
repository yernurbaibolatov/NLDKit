#pragma once

#include "DynamicalSystem.hpp"

class UserDefinedSystem : public DynamicalSystem {
public:
    UserDefinedSystem() = default;
    ~UserDefinedSystem() override = default;
};