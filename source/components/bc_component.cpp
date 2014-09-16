/*
 *  bc_component.cpp
 *  hector
 *
 *  Created by Ben on 05/26/2011.
 *
 */

#include "components/bc_component.hpp"
#include "core/core.hpp"
#include "h_util.hpp"
#include "visitors/avisitor.hpp"

using namespace std;

//------------------------------------------------------------------------------
/*! \brief Constructor
 */
BlackCarbonComponent::BlackCarbonComponent() {
    BC_emissions.allowInterp( true );
    BC_emissions.name = BLACK_CARBON_COMPONENT_NAME;
}

//------------------------------------------------------------------------------
/*! \brief Destructor
 */
BlackCarbonComponent::~BlackCarbonComponent() {
}

//------------------------------------------------------------------------------
// documentation is inherited
string BlackCarbonComponent::getComponentName() const {
    const string name = BLACK_CARBON_COMPONENT_NAME;
    
    return name;
}

//------------------------------------------------------------------------------
// documentation is inherited
void BlackCarbonComponent::init( Core* coreptr ) {
    logger.open( getComponentName(), false, Logger::DEBUG );
    H_LOG( logger, Logger::DEBUG ) << "hello " << getComponentName() << std::endl;
    core = coreptr;
    
	// Inform core what data we can provide
    core->registerCapability( D_EMISSIONS_BC, getComponentName() );
}

//------------------------------------------------------------------------------
// documentation is inherited
unitval BlackCarbonComponent::sendMessage( const std::string& message,
                                          const std::string& datum,
                                          const message_data info ) throw ( h_exception )
{
    unitval returnval;
    
    if( message==M_GETDATA ) {          //! Caller is requesting data
        return getData( datum, info.date );
        
    } else if( message==M_SETDATA ) {   //! Caller is requesting to set data
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
void BlackCarbonComponent::setData( const string& varName,
                                    const message_data& data ) throw ( h_exception )
{
    H_LOG( logger, Logger::DEBUG ) << "Setting " << varName << "[" << data.date << "]=" << data.value_str << std::endl;
    
    try {
        if( varName ==  D_EMISSIONS_BC  ) {
            H_ASSERT( data.date != Core::undefinedIndex(), "date required" );
            BC_emissions.set( data.date, unitval::parse_unitval( data.value_str, data.units_str, U_KG ) );
        } else {
            H_THROW( "Unknown variable name while parsing " + getComponentName() + ": "
                    + varName );
        }
    } catch( h_exception& parseException ) {
        H_RETHROW( parseException, "Could not parse var: "+varName );
    }
}

//------------------------------------------------------------------------------
// documentation is inherited
void BlackCarbonComponent::prepareToRun() throw ( h_exception ) {
    
    H_LOG( logger, Logger::DEBUG ) << "prepareToRun " << std::endl;
    oldDate = core->getStartDate();
}

//------------------------------------------------------------------------------
// documentation is inherited
void BlackCarbonComponent::run( const double runToDate ) throw ( h_exception ) {
    H_ASSERT( !core->inSpinup() && runToDate-oldDate == 1, "timestep must equal 1" );
    oldDate = runToDate;
}

//------------------------------------------------------------------------------
// documentation is inherited
unitval BlackCarbonComponent::getData( const std::string& varName,
                                      const double date ) throw ( h_exception ) {
    
    unitval returnval;
    
    H_ASSERT( date != Core::undefinedIndex(), "Date required for bc_component" );
    
    if( varName == D_EMISSIONS_BC ) {
        returnval = BC_emissions.get( date );
    } else {
        H_THROW( "Caller is requesting unknown variable: " + varName );
    }
    
    return returnval;
}

//------------------------------------------------------------------------------
// documentation is inherited
void BlackCarbonComponent::shutDown() {
	H_LOG( logger, Logger::DEBUG ) << "goodbye " << getComponentName() << std::endl;
    logger.close();
}

//------------------------------------------------------------------------------
// documentation is inherited
void BlackCarbonComponent::accept( AVisitor* visitor ) {
    visitor->visit( this );
}