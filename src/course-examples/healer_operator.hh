#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

/** a local operator for solving the equation
 *
 *   - \Delta u + a*u = f   in \Omega
 *                  u = g   on \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 * \tparam BCType parameter class indicating the type of boundary condition
 */
template<class BCType>
class HealerLocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<HealerLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianVolume<HealerLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<HealerLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianBoundary<HealerLocalOperator<BCType> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };                                // assemble boundary

  HealerLocalOperator(const BCType& bctype_, // boundary cond.type
                          unsigned int intorder_=2) :
     bctype( bctype_ ), intorder( intorder_ )
   {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // assume Galerkin: lfsu == lfsv
    // This yields more efficient code since the local functionspace only
    // needs to be evaluated once, but would be incorrect for a finite volume
    // method

    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // select the two components (assume Galerkin scheme U=V)
    typedef typename LFSU::template Child<0>::Type LFSU0;       // extract components
    const LFSU0& lfsu0 = lfsu.template child<0>();           // with template magic
    typedef typename LFSU::template Child<1>::Type LFSU1;
    const LFSU1& lfsu1 = lfsu.template child<1>();

    // extract some types
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType Jacobian;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType Range;
    typedef Dune::FieldVector<RF,dimw> Gradient;
    typedef typename LFSU::Traits::SizeType size_type;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>&
      rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator
           it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions on reference element
        std::vector<Range> phi0(lfsu0.size());
        lfsu0.finiteElement().localBasis().evaluateFunction(it->position(),phi0);
        std::vector<Range> phi1(lfsu1.size());
        lfsu1.finiteElement().localBasis().evaluateFunction(it->position(),phi1);

        // compute u at integration point
        RF u=0.0;
        for (size_type i=0; i<lfsu0.size(); ++i)
          u += x(lfsu0,i)*phi0[i];
        for (size_type i=0; i<lfsu1.size(); ++i)
          u += x(lfsu1,i)*phi1[i];

        // evaluate gradient of basis functions on reference element
        std::vector<Jacobian> js0(lfsu.size());
        lfsu0.finiteElement().localBasis().evaluateJacobian(it->position(),js0);
        std::vector<Jacobian> js1(lfsu.size());
        lfsu1.finiteElement().localBasis().evaluateJacobian(it->position(),js1);

        // transform gradients from reference element to real element
        const typename EG::Geometry::JacobianInverseTransposed
  		  jac = eg.geometry().jacobianInverseTransposed(it->position());
	    std::vector<Dune::FieldVector<RF,dim> > gradphi0(lfsu0.size());
	    for (size_type i=0; i<lfsu0.size(); i++)
		  jac.mv(js0[i][0],gradphi0[i]);
	    std::vector<Dune::FieldVector<RF,dim> > gradphi1(lfsu1.size());
	    for (size_type i=0; i<lfsu1.size(); i++)
		  jac.mv(js1[i][0],gradphi1[i]);

        // compute gradient of u_0, u_1
        Dune::FieldVector<RF,dim> gradu0(0.0);
        for (size_type i=0; i<lfsu0.size(); i++)
          gradu0.axpy(x(lfsu0,i),gradphi0[i]);
        Dune::FieldVector<RF,dim> gradu1(0.0);
        for (size_type i=0; i<lfsu1.size(); i++)
          gradu1.axpy(x(lfsu1,i),gradphi1[i]);

        // define Elasticity coefficients;
        double elasticityCoeff= 1.0;

        // integrate both components
        RF factor = it->weight()*eg.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu0.size(); i++)
             r.accumulate(lfsu0,i,elasticityCoeff*(gradu0*gradphi0[i]+gradu1*gradphi1[i])*factor);
        for (size_type i=0; i<lfsu1.size(); i++)
             r.accumulate(lfsu1,i,elasticityCoeff*(gradu0*gradphi0[i]+gradu1*gradphi1[i]+9.81*gradphi1[i])*factor);
       }
  }
  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
	// select the two components (assume Galerkin scheme U=V)
	typedef typename LFSU::template Child<0>::Type LFSU0;       // extract components
	const LFSU0& lfsu0_s = lfsu_s.template child<0>();           // with template magic
	typedef typename LFSU::template Child<1>::Type LFSU1;
	const LFSU1& lfsu1_s = lfsu_s.template child<1>();

	// assume Galerkin: lfsu_s == lfsv_s
    // This yields more efficient code since the local functionspace only
    // needs to be evaluated once, but would be incorrect for a finite volume
    // method

    // some types
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU0::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType Range;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions
    const int dim = IG::dimension;

    // select quadrature rule for face
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>&
      rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin();
         it!=rule.end(); ++it)
      {
        // position of quadrature point in local coordinates of element
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());

        // evaluate basis functions at integration point
        std::vector<Range> phi0(lfsu0_s.size());
        lfsu0_s.finiteElement().localBasis().evaluateFunction(local,phi0);
        std::vector<Range> phi1(lfsu1_s.size());
        lfsu1_s.finiteElement().localBasis().evaluateFunction(local,phi1);

        // evaluate u (e.g. flux may depend on u)
        RF u=0.0;
        for (size_type i=0; i<lfsu1_s.size(); ++i)
          u += x_s(lfsu1_s,i)*phi1[i];

        // evaluate flux boundary condition
        Dune::FieldVector<RF,dim>
          globalpos = ig.geometry().global(it->position());
        RF j;
        if (globalpos[1]<0.5)
          j = 1.0;
        else
          j = -1.0; // some outflow

        // integrate j
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu1_s.size(); ++i)
          r_s.accumulate(lfsu1_s,i,j*phi1[i]*factor);
      }
  }

private:
  const BCType& bctype;
  unsigned int intorder;
};
