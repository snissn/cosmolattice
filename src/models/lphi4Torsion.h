#ifndef LPHI4TORSION_H
#define LPHI4TORSION_H

#include "CosmoInterface/cosmointerface.h"

namespace TempLat
{
    /////////
    // Model name and number of fields
    /////////
    struct ModelPars : public TempLat::DefaultModelPars {
        static constexpr size_t NScalars = 3;    // [phi0, phi1, Torsion]
        static constexpr size_t NPotTerms = 3;   // [phi4, interaction, torsion kinetic]
    };

    #define MODELNAME lphi4Torsion

    template<class R>
    using Model = MakeModel(R, ModelPars);

    class MODELNAME : public Model<MODELNAME>
    {
    private:
        double g, lambda, q;
        // mTorsion: Mass parameter for the torsion field.
        //   - Sets the "stiffness" or oscillation frequency of the torsion field T.
        //   - The quadratic term in the potential is (1/2) * mTorsion^2 * T^2.
        //   - Controls how fast T oscillates and how quickly it decays to zero in absence of sources.
        double mTorsion;

        // lambdaT: Quartic self-interaction coupling for the torsion field.
        //   - Sets the strength of the nonlinear self-interaction of T.
        //   - The quartic term in the potential is (1/4) * lambdaT * T^4.
        //   - When nonzero, introduces nonlinearity, enabling richer dynamics such as field trapping, solitons, or modified defect evolution.
        double lambdaT;


    public:

        // alpha_torsion: Strength of the dynamical coupling between the torsion field and the φ₀ field (main scalar).
        //   - Controls the magnitude of the additional term in the equation of motion for φ₀, proportional to Torsion * π₀ (where π₀ is the conjugate momentum of φ₀).
        //   - Physically, this parameter governs how strongly the torsion field can amplify, damp, or otherwise modify the dynamics of φ₀.
        //   - Set to zero to recover standard φ⁴ theory; increase to probe Einstein–Cartan-inspired effects or test amplification scenarios.
        double alpha_torsion; // Torsion coupling strength

        MODELNAME(ParameterParser& parser, RunParameters<double>& runPar, std::shared_ptr<MemoryToolBox> toolBox):
        Model<MODELNAME>(parser, runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL))
        {
            lambda = parser.get<double>("lambda");
            q = parser.get<double>("q");
            g = sqrt(q * lambda);

            // Torsion related parameters
            alpha_torsion = parser.get<double>("alpha_torsion", 0.0); // coupling strength
            mTorsion = parser.get<double>("mTorsion", 1.0);      // torsion mass parameter
            lambdaT  = parser.get<double>("lambdaT", 0.0);       // torsion quartic self-coupling

            // Initial conditions (homogeneous field values)
            fldS0 = parser.get<double, 3>("initial_amplitudes");     // [phi0, phi1, Torsion]
            piS0 = parser.get<double, 3>("initial_momenta", {0, 0, 0});

            alpha = 1;
            fStar = fldS0[0];
            omegaStar = sqrt(lambda) * fStar;

            setInitialPotentialAndMassesFromPotential();
        }

        /////////
        // Program potential (3 terms)
        /////////
        auto potentialTerms(Tag<0>) // phi0 quartic (phi4)
        {
            return 0.25 * pow<4>(fldS(0_c));
        }
        auto potentialTerms(Tag<1>) // Interaction between phi0 and phi1
        {
            return 0.5 * q * pow<2>(fldS(0_c) * fldS(1_c));
        }

        // potentialTerms(Tag<2>):
        //   Returns the torsion field's contribution to the potential energy density.
        //   Includes both a quadratic (mass) term and a quartic (self-interaction) term:
        //      V(T) = (1/2) * mTorsion^2 * T^2 + (1/4) * lambdaT * T^4
        //   - mTorsion: mass parameter for torsion field T
        //   - lambdaT: self-coupling strength for torsion field T
        auto potentialTerms(Tag<2>)
        {
            return 0.5 * mTorsion * mTorsion * pow<2>(fldS(2_c))
                 + 0.25 * lambdaT * pow<4>(fldS(2_c));
        }

        /////////
        // Derivatives of the program potential
        /////////
        auto potDeriv(Tag<0>) // dV/dphi0
        {
            // Standard terms plus coupling to torsion field’s “damping/amplification”
            // phi0: pow<3>(phi0) + q*phi0*phi1^2 + alpha_torsion * torsion * pi0 (see note below)
            return pow<3>(fldS(0_c)) + q * fldS(0_c) * pow<2>(fldS(1_c));
        }
        auto potDeriv(Tag<1>)
        {
            return q * fldS(1_c) * pow<2>(fldS(0_c));
        }
        // potDeriv(Tag<2>):
        //   Returns the derivative of the torsion potential with respect to T (the field value).
        //   Used in the equation of motion for the torsion field:
        //      dV/dT = mTorsion^2 * T + lambdaT * T^3
        auto potDeriv(Tag<2>)
        {
            return mTorsion * mTorsion * fldS(2_c)
                 + lambdaT * pow<3>(fldS(2_c));
        }


        /////////
        // Second derivatives (for mass initialization)
        /////////
        auto potDeriv2(Tag<0>)
        {
            return 3 * pow<2>(fldS(0_c)) + q * pow<2>(fldS(1_c));
        }
        auto potDeriv2(Tag<1>)
        {
            return q * pow<2>(fldS(0_c));
        }
        // potDeriv2(Tag<2>):
        //   Returns the second derivative of the torsion potential with respect to T.
        //   Used for mass and stability analysis and in initialization:
        //      d^2V/dT^2 = mTorsion^2 + 3 * lambdaT * T^2
        auto potDeriv2(Tag<2>)
        {
            return mTorsion * mTorsion
                 + 3.0 * lambdaT * pow<2>(fldS(2_c));
        }

        /////////
        // Custom: Torsion–φ₀ Coupling Implementation Note
        /////////
        // The Einstein–Cartan–inspired coupling term between the torsion field (T)
        // and the main scalar field's momentum (π₀) is **not implemented directly here**.
        // Instead, this term is injected by a custom modification to the evolution kernel:
        //    src/include/CosmoInterface/evolvers/velocityverlet.h
        // Specifically, in the kickScalar() function, for i == 0 (φ₀), the following term is added:
        //      + (w * dt / 2) * alpha_torsion * fldS(2_c) * piS(0_c)
        // This ensures the coupling enters d/dt(π₀) as required for the model.
        // See the corresponding patch or documentation for details.

    };
}

#endif // LPHI4TORSION_H

