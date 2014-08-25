#ifndef _RLS_ESTIMATOR_H
#define _RLS_ESTIMATOR_H

#include<iostream>
#include<cmath>

class RLSestimator{

    public:

        // x[0] = mu is the cofficient of kinetic friction, x[1] = Fc is the estimated cutting force
        vctFixedSizeVector<double,2> x;  
        vctFixedSizeVector<double,2> xlast; // used for motionDetection feature
        // P is the covariance of x        
        vctFixedSizeMatrix<double,2,2> P;
        vctFixedSizeMatrix<double,2,2> Plast;

        // Fest is the estimated total tangential force 
        double Fest;
        double FestLast; // used for motionDetection feature

        // indicate a cutting failure mode
        bool fail;
        // indicate contate state
        bool contactState;
        // true if wam is not moving;
        bool wamMotionState;
        
        // SC: sliding and cutting
        // SWC: sliding without cutting
        // CWS: cutting without sliding
        enum scenerio {IDLE, SC, SWC, CWS};
        scenerio sc;

        // error threshold
        double cuttingFailureThreshold; // defines when cutting failure happends
        double contactThreshold;//

        //system matrix yk = Hk.transpose*x + v
        vctFixedSizeMatrix<double,2,1> Hk;
        // measurement error covariance
       // vctFixedSizeMatrix<double,2,2> Rk;
        double Rk;

        // yk is the measurement tangential force 
        double yk;
        // constructor
        RLSestimator(vctFixedSizeVector<double,2>& xinit): 
        
         x(xinit),
         xlast(xinit),
         P(vct2x2::Eye()),
         Plast(vct2x2::Eye()),   
         Rk(0.5),
         fail(false),
         contactState(true),
         wamMotionState(true), 
         sc(SWC),  
         cuttingFailureThreshold(5),
         contactThreshold(-2),
         Fest(0),
         FestLast(0),
         yk(0){
                Hk[0][0] = 1;
                Hk[1][0] = 1;
         }

        ~RLSestimator(){};
        
        //returns true if the cutter is in contact with the cutting surface
        bool inContact(const double& Fn){
            if(Fn < contactThreshold && contactState == true){
                std::cout<< "Cutter in Contact" << std::endl;
                contactState = false;
            }
            if(Fn > contactThreshold && contactState == false){
                std::cout<< "Cutter not in Contact!" << std::endl;
                contactState = true;
            }

            if(Fn < contactThreshold){
                return true;
            }
            if(Fn > contactThreshold){
                return false;
            }
        }

        // detects if WAM is in motion
        bool WAMnotMoving(const bool& wamNotMoving){
            if(wamNotMoving && wamMotionState == true){
                std::cout<< "WAM STOPPED!" << std::endl;
                wamMotionState = false;
                return true;
            }
            if(!wamNotMoving && wamMotionState == false){
                std::cout<< "WAM moving" << std::endl;
               wamMotionState = true;
                return false;
            }
        } 

        bool isSlidingWithoutCutting(const double& Fc, const bool& wamNotMoving, const double& currtime){
            /*if(wamNotMoving && wamMotionState == true){
                std::cout<< "WAM STOPPED!" << std::endl;
                wamMotionState = false;
                return false;
            }
            else{*/
                if(Fc > 1 && sc == SC){
                    std::cout<< "Cutting and sliding at Time" << currtime <<std::endl;
                    sc = SWC;
                    wamMotionState = true;
                    return false;
                }
                if( Fc <= 1 && sc == SWC){
                    std::cout<< "sliding without cutting at Time" << currtime << std::endl;
                    sc = SC;
                    wamMotionState = true;
                    return true;
                }
            //}    
        }      

        void GetEstimates (vctFixedSizeVector<double,2>& xesti, double& Festi) const{

            xesti[0] = x[0];
            xesti[1] = x[1];
            Festi = Fest; 
        }        
        
        // returns true if detects a cutting failure
        // parameters:
        //
        // Fn: measured normal force
        // Fe: measured tangential force
        // P: covariance matrix from last step 
        // xnew: newly estimated [mu, Fc]
        // Fnew: newly estimated tangential force
        // Cov: updated covariance matrix
        // isCuting: true if sliding and cutting, false if sliding but not cutting
       void Evaluate( const double &Fn, 
                      const double &Fe, 
                      const vctFixedSizeMatrix<double,2,2>& p, 
                      vctFixedSizeVector<double,2>& xnew, 
                      double& Fnew, 
                      vctFixedSizeMatrix<double,2,2>& Cov, 
                      bool isCutting){

               // xnew = xlast;
               // Fnew = Flast;
                
                // in the case of sliding without cutting
               /* if(!isCutting){
                    Hk[1][0] = 0;
                    xnew[1] = 0;
                }*/

                Hk[0][0] = Fn;
                yk = Fe;
                
                vctFixedSizeMatrix<double,2,1> K;
                vctFixedSizeMatrix<double,1,1> tempK;
                vctFixedSizeMatrix<double,1,1> invtempK;

                tempK = Hk.Transpose()*p*Hk + Rk;

                //Inverse(tempK,invtempK);
                //K = p*Hk*invtempK;
                invtempK[0][0] = 1/tempK[0][0];
                K = p*Hk*invtempK;

                xnew = xnew + K*(yk - Hk.Transpose()*xnew);

                Cov = (vct2x2::Eye() - K*Hk.Transpose())*p;
                
                double CovforbeniusNorm;
                
                //calculate frobenius norm of Cov
                CovforbeniusNorm = sqrt(Cov[0][0]*Cov[0][0]+ Cov[1][0]*Cov[1][0]+ Cov[0][1]*Cov[0][1]+ Cov[1][1]*Cov[1][1]);
          
                        
                    //std::cout<< Cov << std::endl;
                    //std::cout<< "froben: "<< CovforbeniusNorm <<std::endl;                
                   //to prevent P of becoming too small due to large estimation error (singularities, etc) and the estimator becomes unresponsive
                if(CovforbeniusNorm < 0.0001){Cov = vct2x2::Eye();}

                //P = (vct2x2::Eye() - K*Hk.Transpose())*P;
                
                Fnew = xnew[0]*Fn + xnew[1];
             
                //xlast = xnew;
                //Flast = Fnew;

               }


        // Evaluate the estimator comparing the outcome of three different cases
        // 1. sliding and cutting (simulataneously estimating mu and Fc)
        // 2. sliding without cutting (Fc = 0)
        bool EvaluateWithComparison( const double &Fn, const double &Fe ){

                vctFixedSizeVector<double,2> x1(xlast);
                vctFixedSizeVector<double,2> x2(xlast);

                double F1 = FestLast;
                double F2 = FestLast;

                vctFixedSizeMatrix<double,2,2> Cov1(Plast);
                vctFixedSizeMatrix<double,2,2> Cov2(Plast);

                // evaluate the first scenerio
                Evaluate(Fn, Fe, Plast, x1, F1, Cov1,true);
                // evaluate the second scenerio
                Evaluate(Fn, Fe, Plast, x2, F2, Cov2, false);
                
                double Fdiff[2] = {std::abs(F1-Fe), std::abs(F2-Fe)};

                // finding the case that best resembles the measured tangential force
                int minIndex = 0;
                double minTemp = 100;
                for(int j = 0; j < 2; ++j){

                    if(Fdiff[j] < minTemp){
                        minIndex = j;
                        minTemp = Fdiff[j];
                    }
                }
                
                // find the most suitable scenerio and update member variables
                if(minIndex == 0){
                    x = x1;
                    Fest = F1;
                    P = Cov1;
                }
                if(minIndex == 1){
                    x = x2;
                    Fest = F2;
                    P = Cov2;
                }

               
                xlast = x;
                FestLast = Fest;
                Plast = P;

              if(fabs(Fest - Fe) > cuttingFailureThreshold){

                return true;
            }
            else{

                return false;
            }


        }


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

        
   void RLSestimate(const double &Fn, 
                    const double &Fe, 
                    const double& currtime,
                    const vctDynamicVector<double>& currJointPos, 
                    vctDoubleVec& rlsEstData,
                    bool& failstate,
                    double& haveFailed,
                    const bool& wamNotMoving,
                    const std::string motionDetection = "OFF"){

        // detect is cutter is in contact
        //inContact(Fn);// Tested and works
        // detects if wam is motionDetection
        WAMnotMoving(wamNotMoving);
        //isSlidingWithoutCutting(x[1], wamNotMoving, currtime);

        if (motionDetection == "OFF"){

                   if((this->EvaluateWithComparison( Fn, Fe ) || Fe > 15) && !failstate){
                   //std::cout<<"Cutting Failure at time:"<< currtime <<std::endl;
                   failstate = true;
                   haveFailed = 10;
                   //to prevent of P becoming too small due to large estimation error (singularities, etc) and the estimator becomes unresponsive
                   P = vct2x2::Eye();
                }
                else{
                    failstate = false;
                    haveFailed = 0;
                }

                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
        }
        if (motionDetection == "ON"){

                
            if(wamNotMoving || !inContact(Fn)){

                x.Assign((double)0,0);
                Fest = 0;
                haveFailed = 0;
                xlast[1] = 0;
                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
      
               }
              else{
                if((this->EvaluateWithComparison( Fn, Fe ) || Fe > 15) && !failstate){
                   //std::cout<<"Cutting Failure at time:"<< currtime <<std::endl;
                   failstate = true;
                   haveFailed = 10;
                   P = vct2x2::Eye();
                }
                else{
                    failstate = false;
                    haveFailed = 0;
                }

                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
            }
    }
 }

      void Inverse(vctFixedSizeMatrix<double,2,2>& M, vctFixedSizeMatrix<double,2,2>& Minv){
            double detM;
            detM = M[0][0]*M[1][1] - M[0][1]*M[1][0];
            
            if(abs(detM) < 0.0001){
                detM = 0.0001;
            }

            Minv[0][0] = M[1][1];
            Minv[0][1] = -M[0][1];
            Minv[1][0] = -M[1][0];
            Minv[1][1] = M[0][0];

            Minv = (1/detM)*Minv;


        }

};




#endif
