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

    ForceBase()
    {
    }

    virtual ~ForceBase()
    {
    }

    virtual void set_rod(RodState const* rod)
    {
        assert(rod != nullptr);
        m_rod = rod;
    }

    void accumulate_energy_force(RodState* rod)
    {
        assert(rod != nullptr);
        set_rod(rod);
        rod->m_potential_energy += get_potential_energy();
        rod->m_force_vector += get_force_vector();
    }

    virtual double get_potential_energy() const = 0;
    virtual VecXd get_force_vector() const = 0;
    virtual BandLimitedMatXd get_jacobian() const = 0;

protected:

    RodState const* m_rod = nullptr;

    friend class Composite;
    // Just for the rod consistency assertion in Composite();
};

}

}
