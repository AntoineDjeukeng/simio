#ifndef CORE_MODEL_H
#define CORE_MODEL_H
#include <stddef.h>
#include <stdint.h>

typedef enum e_fmt { FMT_NONE=0,FMT_GRO,FMT_XTC,FMT_TRR,FMT_LAMMPS_DATA } t_fmt;
typedef enum e_units { UNITS_UNKNOWN=0,UNITS_GROMACS,UNITS_LMP_REAL,UNITS_LMP_METAL,UNITS_LMP_SI } t_units;

typedef struct { double m[3][3]; } t_mat3x3;
typedef struct { double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz; int triclinic; } t_bounds;
typedef struct { t_mat3x3 h; int has_h; t_bounds b; int has_b; } t_box;

typedef struct { char s[9]; }  t_str8;
typedef struct { char s[17]; } t_str16;
typedef struct { char s[33]; } t_str32;

typedef struct {
    int id, mol_id, resid;
    t_str8 atom_name, res_name;
    t_str16 atype;
    double mass, charge;
    int lmp_type;
    int image[3];
} t_atom;

typedef struct { int i,j,type; }        t_bond;
typedef struct { int i,j,k,type; }      t_angle;
typedef struct { int i,j,k,l,type; }    t_dihedral;
typedef struct { int i,j,k,l,type; }    t_improper;

typedef struct {
    size_t natoms; t_atom *atoms;
    size_t nbonds; t_bond *bonds;
    size_t nangles; t_angle *angles;
    size_t ndihedrals; t_dihedral *dihedrals;
    size_t nimpropers; t_improper *impropers;
    size_t n_bond_types; t_str32 *bond_type_name;
    size_t n_angle_types; t_str32 *angle_type_name;
    size_t n_dihedral_types; t_str32 *dihedral_type_name;
    size_t n_improper_types; t_str32 *improper_type_name;
    char *title;
    t_units units;
} t_topology;

typedef enum { CH_NONE=0, CH_EXPLICIT_WALLS, CH_IMPLICIT_SLIT } t_channel_type;
typedef struct {
    t_channel_type type;
    int n_wall; int *wall_ids;
    double nrm[3], d_lo, d_hi; int defined;
} t_channel;

typedef struct {
    size_t natoms; double *x; double *v; double *f;
    int has_v, has_f; t_box box; double time_ps; int64_t step;
} t_frame;

typedef struct {
    int *full_to_dyn; size_t n_full;
    int *dyn_to_full; size_t n_dyn;
} t_index_map;

typedef struct {
    t_topology full; t_channel chan;
    t_frame static_full;
    t_index_map map; size_t n_dyn;
} t_system;

typedef struct {
    t_fmt fmt; int rw; size_t natoms; t_units units; void *impl;
} t_traj;

#endif