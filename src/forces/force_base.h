/*
 * force_base.h
 */

#pragma once

#include "../rods/rod_dynamic.h"

namespace rodeo
{

namespace force
{

class ForceBase
{
public:

    ForceBase(RodState const& rod_state) : // TODO: use split of rod_state to pass only material_frames etc as mutable.
            m_rod(rod_state)
    {
    }

    virtual ~ForceBase()
    {
    }

    virtual double get_energy() const = 0;
    virtual VecXd get_force() const = 0;
    virtual BandLimitedMatXd get_jacobian() const = 0;

protected:

    RodState const& m_rod;

    friend class Composite;
    // Just for the rod consistency assertion in Composite();
};

}

}
