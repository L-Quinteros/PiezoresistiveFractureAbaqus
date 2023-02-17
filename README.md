# An electromechanical phase-field fracture modelling of piezoresistive materials.

The following repository contains an *user-subroutine* and the corresponding *.inp* files to simulate the piezoresistive fracture in *Abaqus*.

![](Fig_diagram.jpg?raw=true)


The system of equations of the solve is the following,

$$
\nabla \cdot \mathbf{\sigma} =\mathbf{0} \quad \text { on } \quad \Omega\\
$$

$$
\nabla \cdot \mathbf{\zeta} - \omega = 0 \quad \text { on } \quad \Omega\\
$$

$$
\nabla \cdot  \mathbf{J} =0 \quad \text{ on } \quad \Omega \\
$$

with the following boundary conditions,

$$
\mathbf{h} =\mathbf{\sigma}  \mathbf{n} \quad \text { on } \quad \partial \Omega_{h}
$$

$$
f_{\phi} = \mathbf{\zeta} \cdot \mathbf{n} \quad \text { on } \quad \partial \Omega_{f}
$$

$$
J_{n} =\mathbf{J} \cdot \mathbf{n} \quad \text{ on } \quad \partial \Omega_{ J_{n} } 
$$


A commonly used one-way coupling formulation for the piezoresistivity is adopted, which implies that the strain and phase-field variable affect the electrical field and not the opposite.
## Constitutive relationship
### Mechanical deformation
Small displacements  will be assumed throughout this work, so the strain field is expressed as:

$$
    \mathbf{\varepsilon}=\nabla^{sym} \mathbf{u}=\frac{1}{2} \left[\nabla \mathbf{u}^{\rm{T}}+\nabla \mathbf{u} \right],
$$


and a linear elastic relationship between the stress and strain is chosen as,

$$
\sigma=h_{1}(\phi)\mathbf{C} \colon \varepsilon,
$$

### Electrical conductivity

The relationship between the electric field and the electric potential is given by,

$$
E = -\nabla \varphi
$$

then, the linear relation between the conductivity and the electric field is.

$$
J = h_{2}(\phi,k,n) \sigma_{eff}(\varepsilon) E
$$

### Phase-field fracture

The Phase-field constitutive equation can be derived following thermodynamically consistent criteria, hence, the scalar microstress work  and the vectorial microstress work is obtained as,

$$
\omega =  \frac{\partial h_{1}}{\partial \phi} \psi_{0} + G_{c} \frac{1}{\ell} \phi
$$


$$
\zeta = G_{c} \ell \nabla \phi
$$

where the strain energy density for the undamaged solid can be stated as,


$$
\psi_{0} = \frac{1}{2} \varepsilon^T \colon C \colon \varepsilon
$$

## How to use the Fortran subroutine in Abaqus
