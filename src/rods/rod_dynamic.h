/*
 * rod_dynamic.h
 *
 *  Created on: 12 mai 2018
 *      Author: jma
 */

#pragma once

#include "rod_state.h"

namespace rodeo
{

class VelocityExtension
{
public:

    VelocityExtension(VecXd const& initial_velocity) :
            m_velocity(initial_velocity)
    {
    }

    virtual ~VelocityExtension()
    {
    }

    inline VecXd const& get_velocity() const
    {
        return m_velocity;
    }

    inline void set_velocity(VecXd const& velocity)
    {
        m_velocity = velocity;
    }

protected:

    VecXd m_velocity;
};

class ExplicitEulerRodState: public StaticRodState<1>, public VelocityExtension
{
public:
    using Solver = ExplicitEulerSolver;

    ExplicitEulerRodState(VecXd const& rest_dofs) :
            StaticRodState<1>(rest_dofs), //
            VelocityExtension(VecXd(rest_dofs.size()))
    {
    }

    ExplicitEulerRodState(VecXd const& rest_dofs, VecXd const& initial_velocity) :
            StaticRodState<1>(rest_dofs), VelocityExtension(initial_velocity)
    {
        assert(rest_dofs.size() == initial_velocity.size());
    }

};

class ImplicitEulerRodState: public StaticRodState<2>, public VelocityExtension
{
public:

    ImplicitEulerRodState(VecXd const& rest_dofs) :
            StaticRodState<2>(rest_dofs), VelocityExtension(
                    VecXd(rest_dofs.size()))
    {
    }

    ImplicitEulerRodState(VecXd const& rest_dofs, VecXd const& initial_velocity) :
            StaticRodState<2>(rest_dofs), VelocityExtension(initial_velocity)
    {
        assert(rest_dofs.size() == initial_velocity.size());
    }

};

using RodState = ExplicitEulerRodState;


}
