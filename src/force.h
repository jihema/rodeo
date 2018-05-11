/*
 * force.h
 */

#pragma once

#include "rodmath.h"

namespace rodeo
{

class Force
{
public:

    Force(RodState& rod_state) : // TODO: use split of rod_state to pass only material_frames etc as mutable.
            m_rod_state(rod_state)
    {
    }

    virtual ~Force()
    {
    }

    virtual double get_energy() = 0;
    virtual VecXd get_force() = 0;
    virtual BandLimitedMatXd get_jacobian() = 0;

    virtual void accumulate_energy()
    {
        m_rod_state.get_energy() += get_energy();
    }
    virtual void accumulate_force(){}
    virtual void accumulate_jacobian(){}

protected:

    RodState& m_rod_state;
};

class Gravitation: public Force
{
public:

    Gravitation(RodState& rod_state, Vec3d const& g = Vec3d(0., 0., -9.81)) :
            Force(rod_state), m_g(g)
    {
    }

private:

    const Vec3d m_g;
};

class RodInternal: public Force
{
public:

};

}
