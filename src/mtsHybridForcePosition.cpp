/*-*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-   */
/*ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab:*/

/*
  $Id: mtsHybridForcePosition.cpp 3181 2011-11-15 15:41:28Z sleonard $

  Author(s):  Simon Leonard
  Created on: 2013-07-22

  (C) Copyright 2012 Johns Hopkins University (JHU), All Rights Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/
/*
    osaJR3ForceSensor::Wrench zero( 0.0 );
    osaJR3ForceSensor 
    jr3.Open();
*/

#include "mtsHybridForcePosition.h"
#include <cisstParameterTypes/prmForceTorqueJointSet.h>
#include <cisstOSAbstraction/osaCPUAffinity.h>

osaJR3ForceSensor::Wrench 
convolve
( const std::list< osaJR3ForceSensor::Wrench >& stdft,
  const vctDynamicVector<double>& sg ){
    
    osaJR3ForceSensor::Wrench ft( 0.0 );
    std::list< osaJR3ForceSensor::Wrench >::const_iterator fti;
    size_t i=0;
    for( fti=stdft.begin(); fti!=stdft.end(); fti++, i++ )
        { ft = ft + sg[i] *(*fti); }

    return ft;
}

mtsHybridForcePosition::mtsHybridForcePosition
( const std::string& name,
  double period,
  const std::string& robotfilename, 
  const vctFrame4x4<double>& Rtw0,
  const vctFrame4x4<double>& Rtnt,
  const vctDynamicVector<double>& qinit,
  const vctDynamicVector<double>& qready,
  
  osaJR3ForceSensor* jr3,
  osaGravityCompensation* gc,
  osaHybridForcePosition* hfp ):
    
    mtsTaskPeriodic( name, period, true ),
    robot( robotfilename, Rtw0),
    tool(Rtnt),
    
    traj( NULL ),
    timer( 0.0 ),
    Rts( Rtnt.Rotation().Transpose() ),
    Rtnt(Rtnt),
    jr3( jr3 ),
    gc( gc ),
    hfp( hfp ),
    
    slave( NULL ),
    control( NULL ),

    qready( qready ),
    qsold( qinit ),
    tauold( 7, 0.0 ),
    
    state( DONOTHING ),
    enable( false ),
//-------------- RLS -------------------
    failstate(false),
    isMoving(false),
    prevTime(osaGetTime()),
    avgNum(100),
    haveFailed(0),
//----------------------------------
    fz( 0.0 ),
    sg( nmrSavitzkyGolay( 1, 0, 100, 0 ) ){

//-------------- For RLS -------------------------------

        // initial value of the estimated coeff. of friction and Fc
        vctFixedSizeVector<double,2> xinit(0.5,1);
        rls = new RLSestimator(xinit);

        ofsForceData.open("/home/lixiao/Desktop/Data1.txt");
        startTime = osaGetTime();
        
        timeStamps.push_back(startTime);
        // initialize jointPoses
        std::vector<double> temp(7, 0);
        for(int i = 0; i < avgNum; i++){jointPoses.push_back(temp);} 

        prevJointPos = qinit;
        
        xesti = xinit;
        rlsEstData.SetSize(6);
        rlsEstData.Assign((double)0,0,0,0,0,0);
        Festi = 0;
    //used to connect to mtsROS         
    rlsProvided = this->AddInterfaceProvided("Controller");
    if(rlsProvided){
        this->StateTable.AddData(rlsEstData,"rlsEstData");
        rlsProvided->AddCommandReadState(this->StateTable, rlsEstData, "GetRLSestimates");
    }
//-------------------------------------------------------
    robot.Attach( &tool );
    
    control = AddInterfaceRequired( "Control" );
    if( control ){
       
        control->AddEventHandlerVoid( &mtsHybridForcePosition::Move,
                                      this,
                                      "Move");
        control->AddEventHandlerVoid( &mtsHybridForcePosition::ToIdle,
                                      this,
                                      "ToIdle");
        control->AddEventHandlerVoid( &mtsHybridForcePosition::PrintTime, this, "PrintTime");
       }



    slave = AddInterfaceRequired( "Slave" );
    if( slave ){
        slave->AddFunction( "GetPositionMSR", mtsGetPosition );
        slave->AddFunction( "SetPositionCMD", mtsSetPosition );
    }

    
 }
 
    void mtsHybridForcePosition::Configure( const std::string& ){}
    void mtsHybridForcePosition::Startup(){
        osaCPUSetAffinity( OSA_CPU2 );
        Thread.SetPriority( 70 );
    }

    void mtsHybridForcePosition::Run(){ 
        
        // Timing stuff
        static double t1 = osaGetTime();
        double t2 = osaGetTime();
        dt.push_back( t2 - t1 );
        t1 = t2;
        
        if( 1.0/GetPeriodicity() < dt.size() ){
            
            std::list<double>::iterator it=dt.begin();
            double avg=0.0;
            double max=0.0;
            for( ; it!=dt.end(); it++ ){
                avg += *it;
                if( max < *it ) max = *it;
            }
           dt.clear();
        }
        
       ProcessQueuedCommands(); 
        ProcessQueuedEvents(); 

        if( state == MOVE )    { MoveTraj();      }
        else if(state == IDLE) { Idle();          }
        else                   { HybridControl(); }

    }

    void mtsHybridForcePosition::Idle(){

        timer = osaGetTime();

              // current joints
        prmPositionJointGet prmq; 
        mtsGetPosition( prmq );
        vctDynamicVector<double> q = prmq.Position();

        // current Cartesian pose
        vctFrame4x4<double> Rtwt = robot.ForwardKinematics( q );
       
            // extract the rotation of the tool
            vctMatrixRotation3<double> Rwt( Rtwt.Rotation() );
            
            // now get the orientation of the sensor
            vctMatrixRotation3<double> Rws( Rwt * Rts );

            osaJR3ForceSensor::Wrench w;
            bool valid = GetWrench( w );

  // ------------------------ For RLS ----------------------------------------


        double currtime = timer - startTime;

        bool wamNotMoving = WAMIsNotMoving(q, currtime);
    
     /** parameters: 
        *
        * Fn: normal force as measured by the force sensor (currently force in the z direction)
        * Fe: tangential force as measured by the force sensor (currently force in the x direction)
        * currtime: current time = osaGetTime() - startTime
        * currJointPos: current joint positions for all 7 joints
        * rlsEstData: is a vector containing [Fe, Fn, mu, Fc, Fest, haveFailed] that gets pushed to ROS (gets changed by the method)
        * failstate: used to indicate cutting failure (true if failed)  (gets changed by the method)
        * haveFailed: a double value used to draw a spike on rqt if failure happens (set to a large number)  (gets changed by the method)
        * wamNotMoving: is true when wam is not moving
        * motionDetection: a switch to toggle between with and without motion detection feature
        */
   rls->RLSestimate(w[2], 
                w[0], 
                currtime,
                q, 
                rlsEstData,
                failstate,
                haveFailed,
                wamNotMoving,
                "ON");

    rls->GetEstimates(xesti, Festi);

      ofsForceData<< timer - startTime << ", "<<w[0]<<", "<<w[2]<<", "<< xesti[0] << ", "<< xesti[1] <<", " << Festi <<std::endl;
             
}

//----------------------------------------------------------------------------
    

   
void mtsHybridForcePosition::MoveTraj(){

            prmPositionJointGet prmq; 
            mtsGetPosition( prmq );
            qready = prmq.Position();

            vctFrame4x4<double> Rtwt = robot.ForwardKinematics( qready );
        

            Rtwtsold = robot.ForwardKinematics( qready );
            Rtwtsoldcmd = robot.ForwardKinematics( qready );
            Rtwtsoldtrj = robot.ForwardKinematics( qready );

            vctFrame4x4<double> Rtwts(Rtwt);

            Rtwtsoldcmd = Rtwt;
            // move along the Y axis for 0.05m
            Rtwts[2-1][4-1] += 0.5; // CHANGE THIS to the correct direction and value

            // create a 10s trajectory from qready to Rtwts
            if( traj != NULL ) { delete traj; }
            traj = new robLinearSE3( Rtwtsoldcmd, Rtwts, 10.0 );
               
           // extract the rotation of the tool
            vctMatrixRotation3<double> Rwt( Rtwt.Rotation() );
            
            // now get the orientation of the sensor
            vctMatrixRotation3<double> Rws( Rwt * Rts );


             jr3->Zero( Rtwtsold );

            // Current force, compensate for sensor orientation
            osaJR3ForceSensor::Wrench ft;
            jr3->Read( ft, Rws, true, 3 );
            
            // filter the reading
            stdft.push_back( ft );
            if( sg.size() < stdft.size() ) { stdft.pop_front(); }
            ft = convolve( stdft, sg );
 
                    
           
            //std::cout<<ft[0]<<", "<<ft[1]<<", "<<ft[2]<<std::endl;
            // state = HYBRID;

        }




    void mtsHybridForcePosition::HybridControl(){
        
        
    }
    
bool mtsHybridForcePosition::GetPosition( vctVec& q ){

    // read the joint positions
    prmPositionJointGet prmq; 
    mtsGetPosition( prmq );
    q = prmq.Position();
    
    bool valid=false;
    prmq.GetValid( valid );

    return valid;
    
}

  bool mtsHybridForcePosition::GetPosition( vctFrm4x4& Rtwt ){

    // current joints
    vctVec q;
    bool valid = GetPosition( q );
    
    // current Cartesian pose
    Rtwt = robot.ForwardKinematics( q );
    
    return valid;

}

  bool mtsHybridForcePosition::GetWrench( osaJR3ForceSensor::Wrench& w ){
    
    vctFrm4x4 Rtwt;
    bool valid = GetPosition( Rtwt );
    
    // current orientation of the tool
    vctRot3 Rwt( Rtwt.Rotation() );
    // orientation of the sensor
    vctRot3 Rws( Rwt * Rts );
    jr3->Read( w, Rws, true, 3 );

    // filter the reading
    stdft.push_back( w );
    if( sg.size() < stdft.size() ) { stdft.pop_front(); }
    w = convolve( stdft, sg );
    /*
    if( 10 < fabs( w[2] ) )
        { std::cout << w[2] << std::endl; }
    */
    return valid;
    
}


