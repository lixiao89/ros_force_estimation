#ifndef _RLS_ESTIMATOR_H
#define _RLS_ESTIMATOR_H

#include<iostream>

class RLSestimator{

    public:

        // x[0] = mu is the cofficient of kinetic friction, x[1] = Fc is the estimated cutting force
        vctFixedSizeVector<double,2> x;  
        vctFixedSizeVector<double,2> xlast; // used for motionDetection feature
        // P is the covariance of x        
        vctFixedSizeMatrix<double,2,2> P;
        // Fest is the estimated total tangential force 
        double Fest;
        double FestLast; // used for motionDetection feature

        // indicate a cutting failure mode
        bool fail;
        // error threshold
        double threshold;

        //system matrix yk = Hk.transpose*x + v
        vctFixedSizeMatrix<double,2,2> Hk;
        // measurement error covariance
        vctFixedSizeMatrix<double,2,2> Rk;
        
        // yk is the measurement vector with first element measured tangential force and second element 0
        vctFixedSizeVector<double,2> yk;
        // constructor
        RLSestimator(vctFixedSizeVector<double,2>& xinit): 
        
         x(xinit),
         xlast(xinit),
         P(vct2x2::Eye()),
         Rk(vct2x2::Eye()),
         fail(false),
         threshold(5),
         Hk(vct2x2::Eye()),
         Fest(0),
         FestLast(0),
         yk(0,0){
                
               Hk[0][0] = 0;
               Hk[0][1] = 0;
               Hk[1][0] = 1;
               Hk[1][1] = 0;
         }

        void GetEstimates(vctFixedSizeVector<double,2>& xesti, double& Festi){

            xesti[0] = x[0];
            xesti[1] = x[1];
            Festi = Fest; 
        }        
        
        // returns true if detects a cutting failure
       bool Evaluate( const double &Fn, const double &Fe ){

                x = xlast;
                Fest = FestLast;

                Hk[0][0] = Fn;
                yk[0] = Fe;
                

                vctFixedSizeMatrix<double,2,2> K;
                vctFixedSizeMatrix<double,2,2> tempK;
                vctFixedSizeMatrix<double,2,2> invtempK;

                tempK = Hk.Transpose()*P*Hk + Rk;

                Inverse(tempK,invtempK);
                K = P*Hk*invtempK;

                x = x + K*(yk - Hk.Transpose()*x);

                P = (vct2x2::Eye() - K*Hk.Transpose())*P;
                
                Fest = x[0]*Fn + x[1];
                
                xlast = x;
                FestLast = Fest;

                if(fabs(Fest - Fe) > threshold){
                
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

        if (motionDetection == "OFF"){
            
            if(this->Evaluate( Fn, Fe ) && !failstate){
                   std::cout<<"Cutting Failure at time:"<< currtime <<std::endl;
                   failstate = true;
                   haveFailed = 70;
                }
                else{
                    haveFailed = 0;
                }

                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
        }
        if (motionDetection == "ON"){

                
            if(wamNotMoving){

                x.Assign((double)0,0);
                Fest = 0;
                haveFailed = 0;
                rlsEstData.Assign((double)Fe, Fn, x[0], x[1], Fest, haveFailed);
      
               }
              else{
                if(this->Evaluate( Fn, Fe ) && !failstate){
                   std::cout<<"Cutting Failure at time:"<< currtime <<std::endl;
                   failstate = true;
                   haveFailed = 70;
                }
                else{
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
