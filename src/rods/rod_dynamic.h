/*
 * rod_dynamic.h
 *
 *  Created on: 12 mai 2018
 *      Author: jma
 */

#pragma once

#include "rod_static.h"

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

    double get_time() const
    {
        return m_time;
    }

    void set_time(double time)
    {
        m_time = time;
    }

protected:
    void print(std::ostream& os) const
     {
         os << "Time = " << get_time();
     }

    VecXd m_velocity;
    double m_time = 0.;
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

    double get_kinetic_energy() const
    {
        return 0.5 * m_velocity.dot(m_mass * m_velocity);
    }

    void print(std::ostream& os) const
    {
        StaticRodState<1>::print(os);
        os << '\t';
        VelocityExtension::print(os);
        os << '\t';
        os << "kinetic energy = " << get_kinetic_energy();
    }
};

inline std::ostream& operator<<(std::ostream& os,
        ExplicitEulerRodState const& rod_state)
{
    rod_state.print(os);

    return os;
}

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

    double get_kinetic_energy()
    {
        return 0.5 * m_velocity.dot(m_mass * m_velocity);
    }

};

using RodState = ExplicitEulerRodState;

}
