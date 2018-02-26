/* Hector -- A Simple Climate Model
   Copyright (C) 2014-2015  Battelle Memorial Institute

   Please see the accompanying file LICENSE.md for additional licensing
   information.
*/
/*
 *  slr_brick_component.cpp
 *  hector
 *
 *  Created by Ben on 31 January 2012.
 *  Modified by Tony W on 17 February 2017 for BRICK SLR.
 *
 */

#include <boost/lexical_cast.hpp>
#include "components/slr_brick_component.hpp"
#include "components/temperature_component.hpp"
#include "core/core.hpp"
#include "core/dependency_finder.hpp"
#include "h_util.hpp"
#include "visitors/avisitor.hpp"
#include <iostream>

namespace Hector {
  
using namespace std;

extern "C" {
//void run_brick_ (int *ns, int *tstep, std::vector<double> *temperature_tseries,
//void run_brick_ (int *ns, int *tstep, double *temperature_tseries,
//                double *beta0_gsic , double *V0_gsic     , double *n_gsic      , 
//                double *Gs0_gsic   , double *Teq_gsic    , double *a_te        ,
//                double *b_te       , double *invtau_te   , double *V0_te       ,
//                double *a_simple   , double *b_simple    , double *alpha_simple,
//                double *beta_simple, double *V0_simple   , double *a_anto      ,
//                double *b_anto     , double *slope_Ta2tg , double *intercept_Ta2Tg,
//                double *b0_dais    , double *slope_dais  , double *mu_dais     ,
//                double *h0_dais    , double *c_dais      , double *chr_dais    ,
//                double *P0_dais    , double *kappa_dais  , double *nu_dais     ,
//                double *f0_dais    , double *gamma_dais  , double *alpha_dais  , 
//                double *Tfrz       , double *rho_w       , double *rho_i       ,
//                double *rho_m      , double *Toc0        , double *Rad0        ,
//                double *Aoc        , double *lf          ,
//                //std::vector<double> *sl_gsic_out, std::vector<double> *sl_te_out, std::vector<double> *sl_gis_out,
//                //std::vector<double> *sl_ais_out , std::vector<double> *sl_out);
//		double *sl_gsic_out, double *sl_te_out, double *sl_gis_out,
//               double *sl_ais_out , double *sl_out);
void run_brick_ (int *ns, double *tstep, double *temperature_tseries,
		double *delta_ocheat_tseries,
		double *beta0_gsic    , double *V0_gsic     , double *n_gsic      ,
		double *Gs0_gsic      , double *Teq_gsic    , double *sl_gsic_out ,
		double *a_te          , double *b_te        , double *invtau_te   , 
                double *V0_te         , double *sl_te_out   , 
		double *c_tee	      , double *a_tee	  , double *rho_tee     ,
	        double *sa_tee        , int *luse_tee       ,
		double *a_simple      , double *b_simple    , double *alpha_simple,
		double *beta_simple   , double *V0_simple   , double *sl_gis_out  ,
		double *vol_gis_out   ,
		double *a_anto        , double *b_anto      , double *slope_Ta2tg ,
		double *intercept_Ta2Tg, 
		int *luse_aisfastdyn  , double daispar[23]  ,
                //double *b0_dais    , double *slope_dais  , double *mu_dais     , 
		//double *h0_dais    , double *c_dais      , double *P0_dais     , 
		//double *kappa_dais  , double *nu_dais    , double *f0_dais     , 
		//double *gamma_dais  , double *alpha_dais , double *Tfrz        , 
		//double *rho_w       , double *rho_i      , double *rho_m       , 
		//double *Toc0        , double *Rad0       , double *Aoc         , 
		//double *lf	    , double *includes_dSLais , double *chr_dais  ,
		double *sl_ais_out    , double *rad_ais_out, double *vol_ais_out ,
		double *disint_ais_out, 
                double *sl_out	      );
}

//------------------------------------------------------------------------------
/*! \brief Constructor
 */
slrBRICKComponent::slrBRICKComponent() {
}

//------------------------------------------------------------------------------
/*! \brief Destructor
 */
slrBRICKComponent::~slrBRICKComponent() {
}

//------------------------------------------------------------------------------
// documentation is inherited
string slrBRICKComponent::getComponentName() const {
    const string name = SLR_BRICK_COMPONENT_NAME;
    return name;
}

//------------------------------------------------------------------------------
// documentation is inherited
void slrBRICKComponent::init( Core* coreptr ) {
    
    logger.open( getComponentName(), false, Logger::getGlobalLogger().getEchoToFile(), Logger::getGlobalLogger().getMinLogLevel() );
    H_LOG( logger, Logger::DEBUG ) << "hello " << getComponentName() << std::endl;
    
    core = coreptr;
    
    // Register the data we can provide
    core->registerCapability( D_SL_RC, getComponentName() );
    core->registerCapability( D_SLR, getComponentName() );
    core->registerCapability( D_SL_RC_NO_ICE, getComponentName() );
    core->registerCapability( D_SLR_NO_ICE, getComponentName() );
    core->registerCapability( D_SLR_GSIC, getComponentName() );
    core->registerCapability( D_SLR_TE, getComponentName() );
    core->registerCapability( D_SLR_GIS, getComponentName() );
    core->registerCapability( D_SLR_AIS, getComponentName() );

    // Register our dependencies - just need global mean surface temperature 
    core->registerDependency( D_GLOBAL_TEMP, getComponentName() );
    core->registerDependency( D_FLUX_MIXED, getComponentName() );
    core->registerDependency( D_FLUX_INTERIOR, getComponentName() );
 
    // Define the BRICK parameters

    // GSIC
    beta0_gsic.set( 0.000577, U_M_YR_C);    // initial mass balance temperature sensitivity
    V0_gsic.set( 0.4 , U_M);                // initial GSIC volume ("initial" = year 1850), m
    Gs0_gsic.set( 0.0 , U_M);               // initial GSIC contribution to sea-level rise ("initial" = year 1850), m
    n_gsic.set( 0.82 , U_UNITLESS);         // exponent for area-volume scaling, -
    Teq_gsic.set( -0.15 , U_DEGC);          // equilibrium temperature (at which there is no GSIC SLR change), deg C
    
    // TE
    a_te.set( 0.5 , U_M_C);                 // temperature sensitivity of equilibrium TE, m/deg C
    b_te.set( 0.0 , U_M);                   // equilibrium TE for temperature Tg=0, m
    invtau_te.set( 0.005 , U_1_YR);         // 1/timescale of thermal expansion, 1/yr
    V0_te.set( 0.0 , U_M);                  // initial sea-level rise due to thermal expansion, m

    // TEE
    a_tee.set( 0.16 , U_KG_M3_K);	    // Global ocean-avg coefficient of thermal expansion, kg/m3/deg C, from Roquet et al. (2015)
    luse_tee.set( 1 , U_UNITLESS);	    // Whether to use the explicit calculation of thermal expansion
 
    // GIS
    a_simple.set( -0.827 , U_M_C);          // temperature sensitivity of equilibrium volume Veq, m SLE/deg C
    b_simple.set( 7.242 , U_M);             // equilibrium volume Veq (m SLE) for temperature Tg=0, m
    alpha_simple.set( 0.000163 , U_1_YR_C); // temperature sensitivity of exponential decay rate, 1/yr 1/deg C
    beta_simple.set( 0.00002845 , U_1_YR);  // exponential decay rate at Tg=0, 1/yr
    V0_simple.set( 7.242 , U_M);            // initial ice-sheet volume ("initial" refers to 1850), m
    
    // AIS
    a_anto.set( 0.26 , U_DEGC_DEGC);        // sensitivity of Antarctic ocean subsurf temp to global temp, deg C/deg C
    b_anto.set( 0.62 , U_DEGC);             // Antarctic ocean subsurf temp when global mean temp = 0, deg C
    b0_dais.set( 775.0 , U_M);              // undisturbed bed height at the continent center, m
    slope_dais.set( 0.0006 , U_UNITLESS);   // slope of ice sheet bed before loading, -
    mu_dais.set( 8.7 , U_M05);              // profile parameter for parabolic ice surface (related to ice stress), m0.5
    h0_dais.set( 1471.0 , U_M);             // height of runoff line at AIS surface temperaure of 0 deg C, m
    c_dais.set( 95.0 , U_M_C);              // temperature sensitivity of runoff line height, m/deg C
    P0_dais.set( 0.35 , U_M);               // annual precipitation for AIS surf temp Ta=0, m (ice equivalent)
    kappa_dais.set( 0.04 , U_1_DEGC);       // coefficient for exponential dependency of precip on Ta, 1/degC
    nu_dais.set( 0.012 , U_1_M05_YR05);     // proportionality constant relating runoff to precip, 1/m0.5 1/yr0.5
    f0_dais.set( 1.2 , U_M_YR);             // proportionality constant for ice flow at grounding line, m/yr
    gamma_dais.set( 2.5 , U_UNITLESS);      // power for relation of ice flow speed to water depth, -
    alpha_dais.set( 0.5 , U_UNITLESS);      // partition parameter for efect of ocean subsurf temp on ice flux, -
    luse_aisfastdyn.set( 1, U_UNITLESS);    // whether to use the AIS fast dynamics emulator
    Tcrit_dais.set( -15.0, U_DEGC);         // trigger temperature, at which disintegration occurs, deg C
    lambda_dais.set( 0.01, U_M_YR);         // disintegration rate, m/yr
}

//------------------------------------------------------------------------------
// documentation is inherited
unitval slrBRICKComponent::sendMessage( const std::string& message,
                                  const std::string& datum,
                                  const message_data info ) throw ( h_exception )
{
    unitval returnval;
    
    if( message==M_GETDATA ) {          //! Caller is requesting data
        return getData( datum, info.date );
        
    } else if( message==M_SETDATA ) {   //! Caller is requesting to set data
        H_THROW("SLR BRICK sendMessage not yet implemented for message=M_SETDATA.");
        //TODO: call setData below
        //TODO: change core so that parsing is routed through sendMessage
        //TODO: make setData private
        
    } else {                        //! We don't handle any other messages
        H_THROW( "Caller sent unknown message: "+message );
    }
    
    return returnval;
}

//------------------------------------------------------------------------------
// documentation is inherited
void slrBRICKComponent::setData( const string& varName,
                            const message_data& data ) throw ( h_exception )
{
    using namespace boost;

    H_LOG( logger, Logger::DEBUG ) << "Setting " << varName << "[" << data.date << "]=" << data.value_str << std::endl;
    
    try {
        if( varName == D_BETA0_GSIC ) {
            H_ASSERT( data.date == Core::undefinedIndex() , "date not allowed" );
            beta0_gsic = unitval::parse_unitval( data.value_str, data.units_str, U_M_YR_C );
        } else if( varName == D_V0_GSIC ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            V0_gsic = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_GS0_GSIC ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            Gs0_gsic = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_N_GSIC ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            n_gsic = unitval::parse_unitval( data.value_str, data.units_str, U_UNITLESS );
        } else if( varName == D_TEQ_GSIC ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            Teq_gsic = unitval::parse_unitval( data.value_str, data.units_str, U_DEGC );
        } else if( varName == D_A_TE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            a_te = unitval::parse_unitval( data.value_str, data.units_str, U_M_C );
        } else if( varName == D_B_TE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            b_te = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_INVTAU_TE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            invtau_te = unitval::parse_unitval( data.value_str, data.units_str, U_1_YR );
        } else if( varName == D_V0_TE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            V0_te = unitval::parse_unitval( data.value_str, data.units_str, U_M );
	} else if( varName == D_A_TEE ) {
	    H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            a_tee = unitval::parse_unitval( data.value_str, data.units_str, U_KG_M3_K );
	} else if( varName == D_LUSE_TEE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            luse_tee = unitval::parse_unitval( data.value_str, data.units_str, U_UNITLESS );
        } else if( varName == D_A_SIMPLE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            a_simple = unitval::parse_unitval( data.value_str, data.units_str, U_M_C );
        } else if( varName == D_B_SIMPLE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            b_simple = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_ALPHA_SIMPLE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            alpha_simple = unitval::parse_unitval( data.value_str, data.units_str, U_1_YR_C );
        } else if( varName == D_BETA_SIMPLE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            beta_simple = unitval::parse_unitval( data.value_str, data.units_str, U_1_YR );
        } else if( varName == D_V0_SIMPLE ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            V0_simple = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_A_ANTO ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            a_anto = unitval::parse_unitval( data.value_str, data.units_str, U_DEGC_DEGC );
        } else if( varName == D_B_ANTO ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            b_anto = unitval::parse_unitval( data.value_str, data.units_str, U_DEGC );
        } else if( varName == D_B0_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            b0_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_SLOPE_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            slope_dais = unitval::parse_unitval( data.value_str, data.units_str, U_UNITLESS );
        } else if( varName == D_MU_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            mu_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M05 );
        } else if( varName == D_H0_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            h0_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_C_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            c_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M_C );
        } else if( varName == D_P0_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            P0_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M );
        } else if( varName == D_KAPPA_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            kappa_dais = unitval::parse_unitval( data.value_str, data.units_str, U_1_DEGC );
        } else if( varName == D_NU_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            nu_dais = unitval::parse_unitval( data.value_str, data.units_str, U_1_M05_YR05 );
        } else if( varName == D_F0_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            f0_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M_YR );
        } else if( varName == D_GAMMA_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            gamma_dais = unitval::parse_unitval( data.value_str, data.units_str, U_UNITLESS );
        } else if( varName == D_ALPHA_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            alpha_dais = unitval::parse_unitval( data.value_str, data.units_str, U_UNITLESS );
       } else if( varName == D_LUSE_AISFASTDYN ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            luse_aisfastdyn = unitval::parse_unitval( data.value_str, data.units_str, U_UNITLESS );            
       } else if( varName == D_TCRIT_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            Tcrit_dais = unitval::parse_unitval( data.value_str, data.units_str, U_DEGC );
        } else if( varName == D_LAMBDA_DAIS ) {
            H_ASSERT( data.date == Core::undefinedIndex(), "date not allowed" );
            lambda_dais = unitval::parse_unitval( data.value_str, data.units_str, U_M_YR );
        } else {
            H_THROW( "Unknown variable name while parsing " + getComponentName() + ": "
                    + varName );
        }
    } catch( bad_lexical_cast& castException ) {
        H_THROW( "Could not convert var: "+varName+", value: " + data.value_str + ", exception: "
                +castException.what() );
    } catch( h_exception& parseException ) {
        H_RETHROW( parseException, "Could not parse var: "+varName );
    }
}

//------------------------------------------------------------------------------
// documentation is inherited
void slrBRICKComponent::prepareToRun() throw ( h_exception ) {
    
    H_LOG( logger, Logger::DEBUG ) << "prepareToRun " << std::endl;

	//H_ASSERT( refperiod_high >= refperiod_low, "bad refperiod" );

// Tony TODO
// set in here the constants for BRICK
// define them in here: ./headers/components/temperature_component.hpp

    dt = 1.0;                // years per timestep (this is hard-coded into Hector)
    slope_Ta2Tg = 0.8364527; // slope and intercept between global mean and Antarctic surface temp
    intercept_Ta2Tg = 15.4235;   // from paleo reconstructions of Shaffer (2014), Ruckert et al (2017)
    Tfrz = -1.8;             // freezing temperature of sea water, deg C
    rho_w = 1030.0;          // sea water density, kg/m3
    rho_i = 917.0;           // ice density, kg/m3
    rho_m = 4000.0;          // rock density, kg/m3
    Toc0 = 0.72;             // Antarctic ocean subsurface reference temperature, deg C
    Rad0 = 1864000.0;        // Antarctic ice sheet reference radius, m
    Aoc = 3.619*pow(10.0,14.0); // ocean surface area, m2
    lf = -1.18;              // AIS shore mean AIS sea level fingerprint, -
    includes_dSLais = 0.0;   // whether AIS is part of SLR rate passed to AIS model. l244 in BRICK_coupledModel.R sets to 0 for full-BRICK runs, -
    Vmin = 19e15;            // minimum volume, below which there is no more ice to disintegrate (m^3), 18-20e15 

    c_tee   = 3991.86795711963; // heat capacity of conservative temperatures, from McDougle (2013)
    rho_tee = 1027.;	     // ocean-average density kg/m3
    sa_tee  = Aoc;            // ocean surface area, m2, from Eakins and Sharmin 2010
    secs_per_Year = 31556926.0;
    interior_area_frac = 0.95; // fraction of ocean surface area at top of interior layer (from DOECLIM)
    
    daispar[11] = Tfrz;
    daispar[12] = rho_w;
    daispar[13] = rho_i;
    daispar[14] = rho_m;
    daispar[15] = Toc0;
    daispar[16] = Rad0;
    daispar[17] = Aoc;
    daispar[18] = lf;
    daispar[19] = includes_dSLais;
    daispar[22] = Vmin;   
    // Convert unitval input parameters to doubles to pass to fortran function
    // GSIC
    beta0_gsic_dbl = double( beta0_gsic.value( U_M_YR_C ) );
    V0_gsic_dbl = double( V0_gsic.value( U_M ) );
    n_gsic_dbl = double( n_gsic.value( U_UNITLESS ) );
    Gs0_gsic_dbl = double( Gs0_gsic.value( U_M ) );
    Teq_gsic_dbl = double( Teq_gsic.value( U_DEGC ) );   
    // TE
    a_te_dbl = double( a_te.value( U_M_C ) );
    b_te_dbl = double( b_te.value( U_M ) );
    invtau_te_dbl = double( invtau_te.value( U_1_YR ) );
    V0_te_dbl = double( V0_te.value( U_M ) );
    // TEE
    a_tee_dbl = double( a_tee.value( U_KG_M3_K ) );
    luse_tee_int = int( luse_tee.value( U_UNITLESS ) );
    // GIS
    a_simple_dbl = double( a_simple.value( U_M_C ) );
    b_simple_dbl = double( b_simple.value( U_M ) );
    alpha_simple_dbl = double( alpha_simple.value( U_1_YR_C ) );
    beta_simple_dbl = double( beta_simple.value( U_1_YR ) );
    V0_simple_dbl = double( V0_simple.value( U_M ) );
    
    a_anto_dbl = double( a_anto.value( U_DEGC_DEGC ) );
    b_anto_dbl = double( b_anto.value( U_DEGC ) );
    // DAIS
    b0_dais_dbl = double( b0_dais.value( U_M ) );
    slope_dais_dbl = double( slope_dais.value( U_UNITLESS ) );
    mu_dais_dbl = double( mu_dais.value( U_M05 ) );
    h0_dais_dbl = double( h0_dais.value( U_M ) );
    c_dais_dbl = double( c_dais.value( U_M_C ) );
    P0_dais_dbl = double( P0_dais.value( U_M ) );
    kappa_dais_dbl = double( kappa_dais.value( U_1_DEGC ) );
    nu_dais_dbl = double( nu_dais.value( U_1_M05_YR05 ) );
    f0_dais_dbl = double( f0_dais.value( U_M_YR ) );
    gamma_dais_dbl = double( gamma_dais.value( U_UNITLESS ) );
    alpha_dais_dbl = double( alpha_dais.value( U_UNITLESS ) );
    luse_aisfastdyn_int = int( luse_aisfastdyn.value( U_UNITLESS ) );
    Tcrit_dais_dbl = double( Tcrit_dais.value( U_DEGC ) );
    lambda_dais_dbl = double( lambda_dais.value( U_M_YR ) );
    daispar[0] = b0_dais_dbl;
    daispar[1] = slope_dais_dbl;
    daispar[2] = mu_dais_dbl;
    daispar[3] = h0_dais_dbl;
    daispar[4] = c_dais_dbl;
    daispar[5] = P0_dais_dbl;
    daispar[6] = kappa_dais_dbl;
    daispar[7] = nu_dais_dbl;
    daispar[8] = f0_dais_dbl;
    daispar[9] = gamma_dais_dbl;
    daispar[10] = alpha_dais_dbl;
    daispar[20] = Tcrit_dais_dbl;
    daispar[21] = lambda_dais_dbl;

    // Initializing all model components that depend on the number of timesteps (ns)
    ns = core->getEndDate() - core->getStartDate();
    
    tgav.resize(ns);
    delta_ocheat.resize(ns);

    vol_gis_out.resize(ns);
    rad_ais_out.resize(ns);
    vol_ais_out.resize(ns);
    disint_ais_out.resize(ns);
    sl_gsic_out.resize(ns);
    sl_te_out.resize(ns);
    sl_gis_out.resize(ns);
    sl_ais_out.resize(ns);
    sl_out.resize(ns);
    
    for(int i=0; i<ns; i++) {
        tgav[i] = 0.0;
        delta_ocheat[i] = 0.0;
    }   
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// documentation is inherited
void slrBRICKComponent::run( const double runToDate ) throw ( h_exception ) {

    H_LOG( logger, Logger::DEBUG ) << "SLR run " << runToDate << std::endl;
 
    int tstep = runToDate - core->getStartDate() - 1;

    //tgav.set( runToDate, core->sendMessage( M_GETDATA, D_GLOBAL_TEMP ) );	// store global temperature
    tgav[tstep] = double(core->sendMessage( M_GETDATA, D_GLOBAL_TEMP ).value( U_DEGC ));   
    delta_ocheat[tstep] = double(core->sendMessage( M_GETDATA, D_HEAT_FLUX ).value( U_W_M2 )) * 
			  sa_tee * secs_per_Year;		// store delta ocheat
    
    if(tstep > 0){
	//H_LOG( logger, Logger::DEBUG ) << "executing run_brick_ " << runToDate << " " << ns << std::endl;
        //H_LOG( logger, Logger::DEBUG ) << "memory location vals, dais_params: " << "b0_dais: " << *(&b0_dais_dbl) << ". slope_dais: " << *(&slope_dais_dbl) << ". mu_dais: " << *(&mu_dais_dbl) << 
	//". h0_dais: " << *(&h0_dais_dbl) << ". c_dais: " << *(&c_dais_dbl) << 
	//". P0_dais: " << *(&P0_dais_dbl) << ". kappa_dais: " << *(&kappa_dais_dbl) << 
	//". nu_dais: " << *(&nu_dais_dbl) << ". f0_dais: " << *(&f0_dais_dbl) << 
	//". gamma_dais: " << *(&gamma_dais_dbl) << ". alpha_dais: " << *(&alpha_dais_dbl) << 
	//". Tfrz: " << *(&Tfrz) << ". rho_w: " << *(&rho_w) << ". rho_i: " << *(&rho_i) << 
	//". rho_m: " << *(&rho_m) << ". Toc0: " << *(&Toc0) << ". Rad0: " << *(&Rad0) << 
	//". Aoc: " << *(&Aoc) << ". lf: " << *(&lf) << ". includes_dSLais: " << *(&includes_dSLais) << 
	//". chr_dais: " << *(&chr_dais_dbl) << std::endl;
	//H_LOG( logger, Logger::DEBUG ) << daispar[0] << "| " << daispar[1] << "| " << daispar[2] 
	//	<< "| " << daispar[3] <<
	//	daispar[4] << "| " << daispar[5] << "| " << daispar[6] << "| " << daispar[7] << "| " << 
	//	daispar[8] << "| " << daispar[9] << "| " <<daispar[10] << "| " << daispar[11] << "| " << 
	//	daispar[12] << "| " << daispar[13] << "| " << daispar[14] << 
	//	"| " << daispar[15] << "| " << daispar[16] << "| " << daispar[17] << "| " << 
	//	daispar[18] << "| " << daispar[19] << "| " << daispar[20] << "| " << std::endl;
  //-------------
        run_brick_ (&ns         , &dt              , &tgav[0]          ,
		&delta_ocheat[0],
                &beta0_gsic_dbl , &V0_gsic_dbl     , &n_gsic_dbl       ,
                &Gs0_gsic_dbl   , &Teq_gsic_dbl    , &sl_gsic_out[0]   ,
		&a_te_dbl       , &b_te_dbl        , &invtau_te_dbl    , 
		&V0_te_dbl      , &sl_te_out[0]    , 
		&c_tee		, &a_tee_dbl	   , &rho_tee	       ,
		&sa_tee		, &luse_tee_int    , &a_simple_dbl     , 
		&b_simple_dbl   , &alpha_simple_dbl, &beta_simple_dbl  , 
		&V0_simple_dbl  , &sl_gis_out[0]   , &vol_gis_out[0]   ,
		&a_anto_dbl     , &b_anto_dbl      , &slope_Ta2Tg      , 
		&intercept_Ta2Tg,
                //&b0_dais_dbl    , &slope_dais_dbl  , &mu_dais_dbl      ,
                //&h0_dais_dbl    , &c_dais_dbl      , &P0_dais_dbl      , 
                //&kappa_dais_dbl , &nu_dais_dbl     , &f0_dais_dbl      , 
                //&gamma_dais_dbl , &alpha_dais_dbl  , &Tfrz             , 
                //&rho_w          , &rho_i           , &rho_m            , 
                //&Toc0           , &Rad0            , &Aoc              , 
                //&lf             , &includes_dSLais , &chr_dais_dbl     ,
                &luse_aisfastdyn_int, &daispar[0]      ,
                &sl_ais_out[0]  , &rad_ais_out[0]  , &vol_ais_out[0]   ,
                &disint_ais_out[0], &sl_out[0]       );
  //-------------

 
        sl_rc.set ( sl_out[tstep] - sl_out[tstep-1], U_M, 0.0 );
        sl_rc_no_ice.set ( sl_te_out[tstep] - sl_te_out[tstep-1], U_M, 0.0 );
	slr.set ( sl_out[tstep], U_M, 0.0 );
	slr_no_ice.set ( sl_te_out[tstep], U_M, 0.0 );
	slr_gsic.set ( sl_gsic_out[tstep], U_M, 0.0 );
	slr_te.set ( sl_te_out[tstep], U_M, 0.0 );
	slr_gis.set ( sl_gis_out[tstep], U_M, 0.0 );
	slr_ais.set ( sl_ais_out[tstep], U_M, 0.0 );	  
    } else {
	sl_rc.set ( 0.0, U_M, 0.0 );
        sl_rc_no_ice.set ( 0.0, U_M, 0.0 );
        slr.set ( 0.0, U_M, 0.0 );
        slr_no_ice.set ( 0.0, U_M, 0.0 );
        slr_gsic.set ( 0.0, U_M, 0.0 );
        slr_te.set ( 0.0, U_M, 0.0 );
        slr_gis.set ( 0.0, U_M, 0.0 );
        slr_ais.set ( 0.0, U_M, 0.0 );
    }
        H_LOG( logger, Logger::DEBUG ) << runToDate <<
                ". SLR_GIS: " << slr_gis << ". tstep: " <<
		tstep << ". tgav[tstep]: " << tgav[tstep] << std::endl;
}

//------------------------------------------------------------------------------
// documentation is inherited
unitval slrBRICKComponent::getData( const std::string& varName,
                              const double date ) throw ( h_exception ) {
    
    unitval returnval;
    
    //H_ASSERT( date != Core::undefinedIndex(), "Date required for all slr data" );
    H_ASSERT( date == Core::undefinedIndex(), "Only current temperatures provided" );
 
    if( varName == D_SL_RC ) {
        returnval = sl_rc;
    } else if( varName == D_SL_RC_NO_ICE ) {
        returnval = sl_rc_no_ice;
    } else if( varName == D_SLR ) {
        returnval = slr;
    } else if( varName == D_SLR_NO_ICE ) {
        returnval = slr_no_ice;
    } else if( varName == D_SLR_GSIC ) {
        returnval = slr_gsic;
    } else if( varName == D_SLR_TE ) {
        returnval = slr_te;
    } else if( varName == D_SLR_GIS ) {
        returnval = slr_gis;
    } else if( varName == D_SLR_AIS ) {
        returnval = slr_ais;
    } else {
        H_THROW( "Caller is requesting unknown variable: " + varName );
    }
    
    return returnval;
}

//------------------------------------------------------------------------------
// documentation is inherited
void slrBRICKComponent::shutDown() {
	H_LOG( logger, Logger::DEBUG ) << "goodbye " << getComponentName() << std::endl;
    logger.close();
}

//------------------------------------------------------------------------------
// documentation is inherited
void slrBRICKComponent::accept( AVisitor* visitor ) {
    visitor->visit( this );
}

}
