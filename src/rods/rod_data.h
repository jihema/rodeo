/*
 * rod_data.h
 */

#pragma once

#include "../solvers/explicit_euler.h"
#include "rod_dynamic.h"

#include <boost/circular_buffer.hpp>

namespace rodeo
{

namespace force
{
class ForceBase;
}

class RodData
{
public:

    static constexpr size_t epochs = 3; // 0 = past; 1 = present; 2 = future.

    RodData(VecXd const& rest_dofs) :
            m_rest_dofs(rest_dofs), //
            m_rod_states(epochs)
    {
        for (size_t i = 0; i < epochs; ++i)
        {
            m_rod_states.push_back(RodState(m_rest_dofs));
            m_rod_states.back().set_dofs_to_rest();
            m_rod_states.back().pre_compute();
        }
    }

    void forward_step(double const h)
    {
        assert(m_rod_states.full());

        RodState& current_state = m_rod_states[1];
        RodState& future_state = m_rod_states[2];

        m_rod_states.push_back(m_rod_states.front()); // Rotating 2->1; 1->0; 0->2.

        current_state.pre_compute(); // Erases energy/force/jacobian, computes material frames.

        RodState::Solver::compute(m_force, current_state); // If m_force is 0 this will do nothing (force vector stays zero) but we still need to do the time step below for inertia.
        current_state.m_dirty = false;

        std::cout << current_state << '\n';

        RodState::Solver::time_step(future_state, current_state, h);
    }

//    void back_step(double const h)
//    {
//        assert(m_rod_states.full());
//
//        m_rod_states.push_front(m_rod_states.back()); // Reusing allocation.
//        RodState::Solver::time_step(m_rod_states[0], m_rod_states[1], -h);
////        m_rod_states[0].compute();
//    }

    void print(std::ostream& os) const
    {
//        os << "Past: " << m_rod_states[0].get_energy() << "; ";
//        os << "Present: " << m_rod_states[1].get_energy() << "; ";
//        os << "Future: " << m_rod_states[2].get_energy() << ".";
    }

    const VecXd& get_rest_dofs() const
    {
        return m_rest_dofs;
    }

    void set_rest_dofs(const VecXd& restDofs)
    {
        m_rest_dofs = restDofs;
    }

    void set_force(force::ForceBase* force)
    {
        m_force = force;
    }

private:

    VecXd m_rest_dofs;
    boost::circular_buffer<RodState> m_rod_states;
    force::ForceBase* m_force = nullptr;
};

inline std::ostream& operator<<(std::ostream& os, RodData const& rod_data)
{
    rod_data.print(os);
    return os;
}

}

