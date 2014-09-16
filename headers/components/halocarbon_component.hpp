#ifndef HALOCARBON_COMPONENT_HPP
#define HALOCARBON_COMPONENT_HPP
/*
 *  halocarbon_component.hpp
 *  hector
 *
 *  Created by Pralit Patel on 1/20/11.
 *
 */

#include "core/logger.hpp"
#include "data/tseries.hpp"
#include "data/unitval.hpp"
#include "imodel_component.hpp"

//------------------------------------------------------------------------------
/*! \brief Model component for a halocarbon.
 *
 *  A halocarbon model component that simply decays in the atmosphere.  Adapted
 *  from Bill Emanuel's python implementation.
 *
 */
class HalocarbonComponent : public IModelComponent {
    friend class INIRestartVisitor; // To allow creation of a restart file.
    friend class CSVOutputStreamVisitor;
    
public:
    HalocarbonComponent( std::string g );
    virtual ~HalocarbonComponent();
    
    
    // IModelComponent methods
    virtual std::string getComponentName() const;
    
    virtual void init( Core* core );
    
    virtual unitval sendMessage( const std::string& message,
                                const std::string& datum,
                                const message_data info=message_data() ) throw ( h_exception );
    
    virtual void setData( const std::string& varName,
                          const message_data& data ) throw ( h_exception );
    
    virtual void prepareToRun() throw ( h_exception );
    
    virtual void run( const double runToDate ) throw ( h_exception );
    
    virtual void shutDown();
    
    // IVisitable methods
    virtual void accept( AVisitor* visitor );
    
private:
    virtual unitval getData( const std::string& varName,
                            const double valueIndex ) throw ( h_exception );
    
    //! Who are we?
    std::string myGasName;
    
    //! Rate coefficient of loss [relates to year?]
    double tau;
    
    //! Radiative forcing efficiency [W/m^2/pptv]
    unitval rho;
    
    //! Forcing [W/m^2]
    tseries<unitval> hc_forcing;
    
    //    tseries<unitval> concentration; //! Time series of concentration, pptv
    tseries<unitval> emissions;     //! Time series of emissions, pptv
    unitval Ha;                     //! Current (ambient) concentration, pptv
    unitval H0;                     //! Preindustrial concentration, pptv
    
    double molarMass;
    
    //! The last calculated date value
    double calcDate;
    
    //! logger
    Logger logger;

	Core *core;
    double oldDate;
};

#endif // HALOCARBON_COMPONENT_HPP