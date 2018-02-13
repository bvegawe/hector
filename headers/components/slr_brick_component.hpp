/* Hector -- A Simple Climate Model
   Copyright (C) 2014-2015  Battelle Memorial Institute

   Please see the accompanying file LICENSE.md for additional licensing
   information.
*/
#ifndef SLR_BRICK_COMPONENT_H
#define SLR_BRICK_COMPONENT_H
/*
 *  slr_brick_component.hpp
 *  hector
 *
 *  Created by Ben on 31 January 2012.
 *
*/

#include "components/imodel_component.hpp"
#include "core/logger.hpp"
#include "data/tseries.hpp"
#include "data/unitval.hpp"

// Need to forward declare the components which depend on each other
#include "components/temperature_component.hpp"

namespace Hector {

//------------------------------------------------------------------------------
/*! \brief The sea level rise component.
 *
 * Computes sea level rise from mean global temperature, based on
 * Wong et al. (2017) in GMDD.
 */
class slrBRICKComponent : public IModelComponent {
    
public:
    slrBRICKComponent();
    ~slrBRICKComponent();
    
    //! IModelComponent methods
    virtual std::string getComponentName() const;

    virtual void init( Core* core );
    
    virtual unitval sendMessage( const std::string& message,
                                const std::string& datum,
                                const message_data info=message_data() ) throw ( h_exception );
    
    virtual void setData( const std::string& varName,
                          const message_data& ) throw ( h_exception );
    
    virtual void prepareToRun() throw ( h_exception );
    
    virtual void run( const double runToDate ) throw ( h_exception );
    
    virtual void shutDown();

    //! IVisitable methods
    virtual void accept( AVisitor* visitor );

    //! Reference period. No output occurs before we reach refperiod_high
	int refperiod_low, refperiod_high;
	int normalize_year;
	
private:
    virtual unitval getData( const std::string& varName,
                            const double valueIndex ) throw ( h_exception );
//
// TW, BVW -- check cm vs m -- BRICK routine outputs meters, might need to convert to cm for output and play nice with other components
//
    unitval	sl_rc;			// sea level rate of change, cm/yr
    unitval	slr;			// sea level rise, cm
    unitval	slr_gsic;		// sea level rise from glaciers and small ice caps, cm
    unitval	slr_te;			// sea level rise from thermal expansion, cm
    unitval	slr_gis;		// sea level rise from Greenland ice sheet, cm
    unitval	slr_ais;		// sea level rise from Antarctic ice sheet, cm
    unitval	sl_rc_no_ice;		// sea level rate of change, cm/yr, no ice
    unitval	slr_no_ice;		// sea level rise, cm, no ice

    //Time series arrays that are updated with each BRICK time-step
    //
    std::vector<double> tgav;
    std::vector<double> delta_ocheat;
    std::vector<double> vol_gis_out;
    std::vector<double> rad_ais_out;
    std::vector<double> vol_ais_out;
    std::vector<double> sl_out;
    std::vector<double> sl_gsic_out;
    std::vector<double> sl_te_out;
    std::vector<double> sl_gis_out;
    std::vector<double> sl_ais_out;

    double dt;                             // years per timestep (this is hardcoded into Hector)
    int ns;                             // number of timesteps
    double slope_Ta2Tg;                 // slope and intercept between global mean and Antarctic surface temp
    double intercept_Ta2Tg;             // from paleo reconstructions of Shaffer (2014), Ruckert et al (2017)
    double Tfrz;                        // freezing temperature of sea water, deg C
    double rho_w;                       // sea water density, kg/m3
    double rho_i;                       // ice density, kg/m3
    double rho_m;                       // rock density, kg/m3
    double Toc0;                        // reference high latitude ocean subsurface temperature, deg C
    double Rad0;                        // reference Antarctic ice sheet radius, m
    double Aoc;                         // area of ocean surface, m2
    double lf;                          // mean AIS sea-level rise fingerprint at AIS shore, -
    double includes_dSLais;		// whether AIS is part of SLR rate passed to AIS model. l244 in BRICK_coupledModel.R sets to 0 for full-BRICK runs, -
    double c_tee;			// heat capacity of conservative temperatures, from McDougle (2013)
    double rho_tee;			// ocean-average density kg/m3
    double sa_tee;			// same as Aoc above
    double secs_per_Year;		//
    double interior_area_frac;		// fraction of ocean surface area at top of interior layer (from DOECLIM)

    // BRICK model parameters (these should all be calibrated)

    // Glaciers and small ice caps model (from MAGICC; Wigley and Raper, 2005)
    unitval beta0_gsic;                 // initial mass balance temperature sensitivity, m/yr/deg C
    unitval V0_gsic;                    // initial GSIC volume ("initial" = year 1850), m
    unitval Gs0_gsic;                   // initial GSIC contribution to sea-level rise ("initial" = year 1850), m
    unitval n_gsic;                     // exponent for area-volume scaling, -
    unitval Teq_gsic;                   // equilibrium temperature (at which there is no GSIC SLR change), deg C
    double beta0_gsic_dbl;
    double V0_gsic_dbl;
    double Gs0_gsic_dbl;
    double n_gsic_dbl;
    double Teq_gsic_dbl;

    // Thermal expansion (Grinsted et al 2010, Mengel et al 2016)
    unitval a_te;                       // temperature sensitivity of equilibrium TE, m/deg C
    unitval b_te;                       // equilibrium TE for temperature Tg=0, m
    unitval invtau_te;                  // 1/timescale of thermal expansion, 1/yr
    unitval V0_te;                      // initial sea-level rise due to thermal expansion, m
    double a_te_dbl;
    double b_te_dbl;
    double invtau_te_dbl;
    double V0_te_dbl;

    // Explicit thermal expansion, using DOECLIM heat xfer and Roquet et al (2015) thermal expansion coeff.
    unitval a_tee;			// global ocean-avg coefficient of thermal expansion, kg/m3/deg C
    unitval luse_tee;			// whether to use the explicit calculation of thermal expansion
    double a_tee_dbl;
    int luse_tee_int;

    // Greenland ice sheet (SIMPLE, Bakker et al, 2016)
    unitval a_simple;                   // temperature sensitivity of equilibrium volume Veq, m SLE/deg C
    unitval b_simple;                   // equilibrium volume Veq (m SLE) for temperature Tg=0, m
    unitval alpha_simple;               // temperature sensitivity of exponential decay rate, 1/yr 1/deg C
    unitval beta_simple;                // exponential decay rate at Tg=0, 1/yr
    unitval V0_simple;                  // initial ice-sheet volume ("initial" refers to 1850), m
    double a_simple_dbl;
    double b_simple_dbl;
    double alpha_simple_dbl;
    double beta_simple_dbl;
    double V0_simple_dbl;

    // Antarctic ice sheet (DAIS, Shaffer 2014)
    unitval a_anto;                     // sensitivity of Antarctic ocean subsurf temp to global temp, deg C/deg C
    unitval b_anto;                     // Antarctic ocean subsurf temp when global mean temp = 0, deg C
    unitval b0_dais;                    // undisturbed bed height at the continent center, m
    unitval slope_dais;                 // slope of ice sheet bed before loading, -
    unitval mu_dais;                    // profile parameter for parabolic ice surface (related to ice stress), m0.5
    unitval h0_dais;                    // height of runoff line at AIS surface temperaure of 0 deg C, m
    unitval c_dais;                     // temperaure sensitivity of runoff line height, m/deg C
    unitval chr_dais;                   // uncertainty parameter for runoff line height, unitless
    unitval P0_dais;                    // annual precipitation for AIS surf temp Ta=0, m (ice equivalent)
    unitval kappa_dais;                 // coefficient for exponential dependency of precip on Ta, 1/degC
    unitval nu_dais;                    // propportionality constant relating runoff to precip, 1/m0.5 1/yr0.5
    unitval f0_dais;                    // proportionality constant for ice flow at groudning line, m/yr
    unitval gamma_dais;                 // power for relation of ice flow speed to water depth, -
    unitval alpha_dais;                 // partition parameter for efect of ocean subsurf temp on ice flux, -

    double a_anto_dbl;
    double b_anto_dbl;
    double b0_dais_dbl;
    double slope_dais_dbl;
    double mu_dais_dbl;
    double h0_dais_dbl;
    double c_dais_dbl;
    double chr_dais_dbl;
    double P0_dais_dbl;
    double kappa_dais_dbl;
    double nu_dais_dbl;
    double f0_dais_dbl;
    double gamma_dais_dbl;
    double alpha_dais_dbl;
    double daispar[21];	
        //! pointers to other components and stuff
    Core *core;

    //! logger
    Logger logger;
};

}

#endif // SLR_BRICK_COMPONENT_H
