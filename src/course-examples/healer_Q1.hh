template<class GV> void healer_Q1 (const GV& gv)
{
  // <<<1>>> Choose domain and range field type
  static const unsigned int dim = GV::dimension;
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;


  Real time;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  CON con;
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,Real,2> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::no_blocking,3>
  VBE;
  typedef Dune::PDELab::VectorGridFunctionSpace<GV,FEM,dim,VBE,VBE,CON, Dune::PDELab::EntityBlockedOrderingTag> GFS;
  GFS gfs(gv,fem);
  gfs.name("displacement vector");
  gfs.update();

  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  cc.clear();

  typedef Healer_LinearElasticityParameters<GV,Real> LEPI;
  LEPI leparams;

  typedef Dune::PDELab::PowerConstraintsParameters<LEPI,dim>
      VectorConstraints;
  VectorConstraints vector_constraints(leparams);
  Dune::PDELab::constraints(vector_constraints,gfs,cc);
  std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;

  // <<<3>>> Make grid operator
  typedef Dune::PDELab::LinearElasticity<LEPI> LOP;
  LOP lop(leparams);
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(200);

  typedef Dune::PDELab::GridOperator<
    GFS,GFS,        /* ansatz and test space */
    LOP,            /* local operator */
    MBE,            /* matrix backend */
    Real,Real,Real, /* field types for domain, range and jacobian */
    CC,CC           /* constraints transformation  for ansatz and test space */
    > GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  typename GO::Traits::Jacobian jac(go);
  std::cout << jac.patternStatistics() << std::endl;

  // <<<4>>> Make FE function extending Dirichlet boundary conditions, and initial conditions
  typedef typename GO::Traits::Domain U;              // alternative way to extract vector type
  U uold(gfs,0.0);
  typedef Dune::PDELab::LinearElasticityDirichletExtensionAdapter<Healer_LinearElasticityParameters<GV,Real> > LEDEA;
  LEDEA ledea(gv,leparams);
  Dune::PDELab::interpolate(ledea,gfs,uold);

 // <<<5>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,true);

  // <<<6>>> define spatial problem solver
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  SLP slp(go,ls,uold,1e-10);
  Dune::PDELab::Alexander2Parameter<Real> method;
  Dune::PDELab::OneStepMethod<Real,GO,SLP,U,U> osm(method, go, slp);

  std::stringstream basename;
  basename << "healer";
  Dune::PDELab::FilenameHelper fn(basename.str());

  U unew(gfs, 0.0);
  unew = uold;
  while(time<0.1)
  {
	  osm.apply(time,0.01,uold,unew);
	  // <<<7>>> graphical output
	  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
	  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,unew);
	  vtkwriter.write(fn.getName(),Dune::VTK::appendedraw);
	  fn.increment();

	  uold=unew;
	  time+=0.01;
  }
}
