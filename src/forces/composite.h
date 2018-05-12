/*
 * composite.h
 */

#pragma once

#include "force_base.h"
#include <vector>

namespace rodeo
{

namespace force
{

class Composite: public ForceBase
{
public:

    Composite(RodState const& rod, std::vector<ForceBase *> const& components) :
            ForceBase(rod), m_components(components)
    {
        for (ForceBase const* component : m_components)
        {
            assert(&(component->m_rod) == &rod);
        }
    }

    double get_energy() const override
    {
        double energy = 0;
        for (ForceBase const* component : m_components)
        {
            energy += component->get_energy();
        }

        return energy;
    }

    VecXd get_force() const override
    {
        VecXd force(m_rod.get_rest_dofs()->size());

        for (ForceBase const* component : m_components)
        {
            force += component->get_force();
        }

        return force;
    }

    BandLimitedMatXd get_jacobian() const override
    {
        BandLimitedMatXd jacobian(m_rod.get_rest_dofs()->size(),
                m_rod.get_rest_dofs()->size());

        for (ForceBase const* component : m_components)
        {
            jacobian += component->get_jacobian();
        }

        return jacobian;
    }

private:
    std::vector<ForceBase *> m_components;
};

}

}
