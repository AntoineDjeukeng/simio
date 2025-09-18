#include "../../include/utils/units.h"

/* length “unit size” in meters for each style */
static inline double meters_per_length_unit(t_units u)
{
    switch (u) {
        case UNITS_GROMACS:   return 1e-9;   /* nm */
        case UNITS_LMP_REAL:  return 1e-10;  /* Å */
        case UNITS_LMP_METAL: return 1e-10;  /* Å */
        case UNITS_LMP_SI:    return 1.0;    /* m  */
        default:              return 1.0;
    }
}

/* time “unit size” used in velocity for each style */
static inline double seconds_per_time_unit_for_velocity(t_units u)
{
    switch (u) {
        case UNITS_GROMACS:   return 1e-12;  /* ps */
        case UNITS_LMP_REAL:  return 1e-15;  /* fs */
        case UNITS_LMP_METAL: return 1e-12;  /* ps */
        case UNITS_LMP_SI:    return 1.0;    /* s  */
        default:              return 1.0;
    }
}

double ft_len_factor_to_fmt_units(t_units from, t_units to)
{
    const double m_from = meters_per_length_unit(from);
    const double m_to   = meters_per_length_unit(to);
    return m_from / m_to;
}

double ft_vel_factor_to_fmt_units(t_units from, t_units to)
{
    const double m_from = meters_per_length_unit(from);
    const double m_to   = meters_per_length_unit(to);
    const double s_from = seconds_per_time_unit_for_velocity(from);
    const double s_to   = seconds_per_time_unit_for_velocity(to);
    return (m_from / m_to) * (s_to / s_from);
}

