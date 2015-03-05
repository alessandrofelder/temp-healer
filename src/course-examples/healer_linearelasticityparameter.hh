/*
 * healer_linearelasticityparameter.hh
 *
 *  Created on: 27 Feb 2015
 *      Author: root
 */

#ifndef HEALER_LINEARELASTICITYPARAMETER_HH_
#define HEALER_LINEARELASTICITYPARAMETER_HH_

template <class GV, class Real>
class Healer_LinearElasticityParameters
	: public Dune::PDELab::LinearElasticityParameterInterface<Dune::PDELab::LinearElasticityParameterTraits<GV, Real>,Healer_LinearElasticityParameters<GV,Real>>
	  {
	  public:
			typedef Dune::PDELab::LinearElasticityParameterTraits<GV, Real> Traits;
			void f(const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeType & y, double t) const
	        {
				y[0] = 0.0;
				y[1] = -9.81*cos(t);
				y[2] = 0.0;
				return;
	        }

			template<typename I>
			bool isDirichlet(const I & ig,
			const typename Traits::IntersectionDomainType & coord
			) const
			{
				Dune::FieldVector<typename I::ctype, I::dimension> xg = ig.geometry().global(coord);
				if(xg[1]<=1.e-9)
					return true;
				else
					return false;
			}

			void g(const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeType & y) const
			{
				y[0]=0.0;
				y[1]=0.0;
				y[2]=0.0;
				return;
			}

			typename Traits::RangeFieldType
			lambda(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
			{
				return 6923077.0;
			}

			typename Traits::RangeFieldType
			mu(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
			{
				return 4615385.0;
			}



	  };

#endif /* HEALER_LINEARELASTICITYPARAMETER_HH_ */
