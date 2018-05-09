/*
 * rod_state.h
 *
 *  Created on: 9 May 2018
 *      Author: jau
 */

#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <eigen3/Eigen/Core>

#include <stddef.h>
#include <assert.h>
#include <stdexcept>

namespace rodeo
{

using VecXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using MatXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using BandLimitedMatXd = MatXd;
// FIXME

class DirtyAccess: public std::runtime_error
{
public:
    DirtyAccess(std::string const& what) :
            std::runtime_error(what)
    {
    }
};

class RodState
{
public:

    RodState(VecXd const& rest_dofs) :
            m_rest_dofs(&rest_dofs), //
            m_dirty(true), //
            m_dofs(rest_dofs.size()), //
            m_force(rest_dofs.size()), //
            m_jacobian(rest_dofs.size(), rest_dofs.size())
    {
    }

    void set_dofs(VecXd const& dofs)
    {
        m_dirty = true;
        m_dofs = dofs;
    }

    VecXd const& get_dofs() const
    {
        return m_dofs;
    }

    bool is_dirty() const
    {
        return m_dirty;
    }

    void compute();

    double get_energy() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_energy;
    }

    VecXd get_force() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_force;
    }

    BandLimitedMatXd get_jacobian() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_jacobian;
    }

private:

    VecXd const* m_rest_dofs = nullptr;
    bool m_dirty = true;
    VecXd m_dofs;
    double m_energy = 0.;
    VecXd m_force;
    BandLimitedMatXd m_jacobian;
};

}

