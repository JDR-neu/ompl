    #include <ompl/control/SpaceInformation.h>
    #include <ompl/base/goals/GoalState.h>
    #include <ompl/base/spaces/SE2StateSpace.h>
    #include <ompl/control/StatePropagator.h>
    #include <ompl/control/spaces/RealVectorControlSpace.h>


    namespace ob = ompl::base;
    namespace oc = ompl::control;
    
    const double b= 0.54;
    const double dt=0.1;
    const int max_size=10001;
    const double end_controls=2000;



 struct UnicycleState {
        
        // pose
        double x;
        double y;
        double yaw;


        UnicycleState() { x = 0.0; y = 0.0; yaw = 0.0;}
        UnicycleState(double xx, double yy, double zz) : x(xx), y(yy),yaw(zz) {}

        UnicycleState(double* p) { x = p[0]; y = p[1]; yaw = p[2];}

        UnicycleState& operator=(UnicycleState const& copy)
        {
            x = copy.x;
            y = copy.y;
            yaw = copy.yaw;


            return *this;
        }

        UnicycleState(const UnicycleState& copy)
        {
            x = copy.x;
            y = copy.y;
            yaw = copy.yaw;
        }

        ~UnicycleState() {

        } 


    };

    struct UnicycleControl {
        
        // v, translational velocity
        double v; 
        // w, rotational velocity
        double w;



        UnicycleControl() { v = 0.0; w = 0.0;}
        UnicycleControl(double vv, double ww) : v(vv), w(ww) {}

        UnicycleControl(double* p) { v = p[0]; w = p[1];}

        UnicycleControl& operator=(UnicycleControl const& copy)
        {
            v = copy.v;
            w = copy.w;

           return *this;
        }

        UnicycleControl(const UnicycleControl& copy)
        {
            v = copy.v;
            w = copy.w;
        }

        ~UnicycleControl() {


        } 




    };



    class POSQ : public oc::StatePropagator {

    public:

      double *result;
      double *intRes;
      
      oc::Control *control_res;
      oc::RealVectorControlSpace *controlSpace;

      POSQ(const oc::SpaceInformationPtr &si) : oc::StatePropagator(si)
      {
        
        ob::StateSpacePtr space(new ob::SE2StateSpace());
        controlSpace=new oc::RealVectorControlSpace(space,max_size);
        control_res= (*controlSpace).allocControl();
        
        std::cout <<"POSQ Steer function activated"<<std::endl;
      	result= (double*)malloc(sizeof(double)*5);
        intRes= (double*)malloc(sizeof(double)*5);

        std::cout <<"POSQ: result allocated"<<std::endl;

      }
      

      void setRes(double *r) const {

        this->intRes[0]=r[0];
        
        this->intRes[1]=r[1];

        this->intRes[2]=r[2];

        this->intRes[3]=r[3];

        this->intRes[4]=r[4];


      }




      void propagate(const ob::State *state, const oc::Control *control, const double duration, ob::State *result) const
      {
        
         // std::cout<<"Call propagate"<<std::endl;


        const ob::SE2StateSpace::StateType *se2state = state->as<ob::SE2StateSpace::StateType>();
        const double* pos = se2state->as<ob::RealVectorStateSpace::StateType>(0)->values;
        const double rot = se2state->as<ob::SO2StateSpace::StateType>(1)->value;



        oc::Control *control_int;
        control_int= (*controlSpace).allocControl();
        controlSpace->copyControl(control_int, control);
        // controlSpace->printControl(control_res, std::cout);


        oc::RealVectorControlSpace::ControlType *rcontrol = static_cast< oc::RealVectorControlSpace::ControlType*>(control_int);

        // // controlSpace->printControl(rcontrol, std::cout);
        const double* ctrl = rcontrol->as<oc::RealVectorControlSpace::ControlType>()->values;
        double x=pos[0];
        double y=pos[1];
        double yaw=rot;
        
        int i=0;


        while(true){

           if(ctrl[2*i]==end_controls || ctrl[2*i+1]==end_controls){
                     std::cout<<"value Ctrls = "<<ctrl[2*i]<<" "<< ctrl[2*i+1]<<std::endl;
                     break;
                 }

            x =  x + ctrl[2*i] * duration * cos(yaw),
            y =  y + ctrl[2*i] * duration * sin(yaw);
            yaw   =  yaw    + ctrl[2*i+1] * duration;

            i++;
            // std::cout<<"Ctrl "<<ctrl[2*i]<<" "<<ctrl[2*i+1]<<std::endl;
             std::cout<<"Intermediate State added: "<<x<<" "<<y<<" "<<yaw<<" "<<std::endl;
        }


         result->as<ob::SE2StateSpace::StateType>()->setXY(x,y);
         result->as<ob::SE2StateSpace::StateType>()->setYaw(yaw);

         // std::cout<<"Final State added: "<<x<<" "<<y<<" "<<yaw<<" "<<std::endl;






      }

      
      double normangle(double a, double mina) const {

        double ap,minap;
        ap=a;
        minap=mina;
       
       while (ap>= (minap+M_PI*2)){
            ap=ap-M_PI*2;
        }
        while(ap<minap){
            ap=ap+M_PI*2;

        }

        return ap;
      }

      

      double* posctrlstep (double x_c, double y_c, double t_c, 
                   double x_end, double y_end, double t_end, double ct, double b, int dir) const {
        

        /** This function will generate a vector of double as output:

         *  [0] Vl velocity of the left wheel;
         *  [1] Vr velocity of the right wheel;
         *  [2] V translational Velocity;
         *  [3] W Angular Velocity.
         *  [4] EOT End Of Trajectory

        **/


        static double oldBeta;

        double Krho,Kalpha,Kbeta,Kphi,Kv,Vmax,RhoEndCondition;


        Kalpha  = 6.91;
        Kbeta   = -1;
        Kphi    = -1;

        Krho    = 0.2;
        Kv	    = 5.9;


        Vmax    = Krho;


        RhoEndCondition = 0.15;




        if(ct==0){
            oldBeta=0;

        }


        double dx,dy,rho,fRho,alpha,phi,beta,v,w,vl,vr,eot;

        // rho
        eot=1;
        dx=x_end-x_c;
        dy=y_end -y_c;
        rho=sqrt(dx*dx+dy*dy);
        fRho=rho;

        if(fRho>(Vmax/Krho)){
            fRho=Vmax/Krho;
        }


        //alpha

        alpha=atan2(dy,dx)-t_c;
        alpha=normangle(alpha,-M_PI);

        //direction

        if (dir==0){


            if(alpha>(M_PI/2)) {
                fRho=-fRho;
                alpha=alpha-M_PI;
            }else if(alpha<=-M_PI/2){
                fRho=-fRho;
                alpha=alpha+M_PI;
            }
        }
        else if(dir==-1){
            fRho=-fRho;
            alpha=alpha+M_PI;
            if(alpha>M_PI){
                alpha=alpha-2*M_PI;
            }
        }


        //phi

        phi=t_end-t_c;
        phi=normangle(phi, -M_PI);

        beta=normangle(phi-alpha, -M_PI);

        if ((abs(oldBeta-beta)>M_PI)){
            beta=oldBeta;
        }
        oldBeta=beta;

        //set speed



        v=Krho*tanh(Kv*fRho);
        w=Kalpha*alpha+Kbeta*beta;


        if (rho<RhoEndCondition){

            eot=1;
            //    cout<<"eot evaluated... "<<eot<<endl;
        }
        else {
            eot=0;
            //   cout<<"eot evaluated... "<<eot<<endl;

        }

        if(eot){
            w=0.;
        }


        //Convert speed to wheel speed

        vl=v-w*b/2;

        if(abs(vl)>Vmax){

            if(vl<0){
                vl=Vmax*-1;}
            else{
                vl=Vmax;}
        }

        vr=v+w*b/2;

        if(abs(vr)>Vmax){
            if(vr<0){
                vr=Vmax*-1;}
            else{
                vr=Vmax;}
        }





        result[0]=vl;
        result[1]=vr;
        result[2]=(vl+vr)/2;
        result[3]=(vr-vl)/b;
        result[4]=eot;


        return result;

    }



       
      bool steer(const ob::State* from, const ob::State* to, oc::Control* c_result, double& duration) const
      {

       //std::cout <<"POSQ Steer function activated"<<std::endl;
        
        double sl,sr,oldSl,oldSr,t,eot,dSl,dSr,dSm,dSd,vl,vr,enc_l,enc_r,dir;
    	

        dir=1;
        enc_l=0;
        enc_r=0;
        sl=0;
        sr=0;
        oldSl=0;
        oldSr=0;
        eot=0;
        t=0;
        vl=0;
        vr=0;

        double x,y,th;

        double x_fin,y_fin,th_fin;
     
        x= from->as<ob::SE2StateSpace::StateType>()->getX();
        y= from->as<ob::SE2StateSpace::StateType>()->getY();
        th=from->as<ob::SE2StateSpace::StateType>()->getYaw();

        x_fin= to->as<ob::SE2StateSpace::StateType>()->getX();
        y_fin= to->as<ob::SE2StateSpace::StateType>()->getY();
        th_fin= to->as<ob::SE2StateSpace::StateType>()->getYaw();

        // std::cout<<"From state "<<x<<" "<<y<<" "<<" "<<th<<" to "<<" "<<x_fin<<" "<<" "<<y_fin<<" "<<th_fin<<" "<<std::endl;

        double vv,ww;
        vv=0;
        ww=0;


        std::vector<UnicycleControl> controls;
        std::vector<UnicycleState> states;
        controls.clear();
        states.clear();
        
        int n_steps;
        n_steps=0;


        while(eot==0){
                // calculate distance for both wheels
                dSl=sl-oldSl;
                dSr=sr-oldSr;
                dSm=(dSl+dSr)/2;


                dSd=(dSr-dSl)/b;





                x=x+dSm*cos(th+dSd/2);
                y=y+dSm*sin(th+dSd/2);
                th=normangle(th+dSd, -M_PI);


                // intRes= posctrlstep (x,y,th,x_fin,y_fin,th_fin, t,b,dir);
                setRes(posctrlstep (x,y,th,x_fin,y_fin,th_fin, t,b,dir));
                //Save the velocity commands,eot
                vv=intRes[2];
                ww=intRes[3];

                eot=intRes[4];
                vl=intRes[0];
                vr=intRes[1];



                //Increase the timer
                t=t+dt;
    	
        	    // Count the number of steps
        	    n_steps++;

                // keep track of previous wheel position
                oldSl=sl;
                oldSr=sr;


                // increase encoder values
                enc_l=enc_l+dt*vl;
                enc_r=enc_r+dt*vr;

                sl=enc_l;
                sr=enc_r;

       


        if(eot==1){


                // Add current values to the Trajectory
    	        states.push_back(UnicycleState(x,y,th));
                controls.push_back(UnicycleControl(vv,ww));
            
        	    /// save the last state!!!
        	    double xf,yf,yawf,vf,wf;
                vf=intRes[2];
                wf=intRes[3];
                dSl=sl-oldSl;
                dSr=sr-oldSr;
                dSm=(dSl+dSr)/2;
                dSd=(dSr-dSl)/b;
                xf=x+dSm*cos(th+dSd/2);
                yf=y+dSm*sin(th+dSd/2);
                yawf=normangle(th+dSd, -M_PI);

                // std::cout<<"adding UnicycleState EOF"<<std::endl;
    	        states.push_back(UnicycleState(xf,yf,yawf));
                // std::cout<<"adding UnicycleControl EOF  "<<vf<<" "<<wf<<std::endl;
                controls.push_back(UnicycleControl(vf,wf));

                }
                else

                {
                // Add current values to the Trajectory
                // std::cout<<"adding UnicycleState"<<std::endl;
    	        UnicycleState s;
                s.x=x;
                s.y=y;
                s.yaw=th;
                // states.push_back(UnicycleState(x,y,th));
                states.push_back(s);

                // std::cout<<"adding UnicycleControl "<<vv<<" "<<ww<<std::endl;
                UnicycleControl c;
                c.v=vv;
                c.w=ww;
                // controls.push_back(UnicycleControl(vv,ww));
                controls.push_back(c);

                }

        }

       
        oc::RealVectorControlSpace::ControlType *rcontrol = static_cast< oc::RealVectorControlSpace::ControlType*>(control_res);
        
        // std::cout<<"Size control for Steering: "<<(int)controls.size()<<std::endl;

        for(int i=0; i<(int)controls.size();i++){
        //std::cout<<" controls for Steering: " << controls[i].v<<" " << controls[i].w<<std::endl;
    	
    	(*rcontrol).values[2*i]   =  controls[i].v;
        (*rcontrol).values[2*i+1] =  controls[i].w;


        }

        (*rcontrol).values[2*(int)controls.size()]   =  end_controls;
        
        c_result=(*controlSpace).allocControl();
        controlSpace->copyControl(c_result, control_res);
        // controlSpace->printControl(control_res, std::cout);

        duration=(2*(int)controls.size()+1)*dt;

        return true;



    }

      
      bool canSteer() const
    {
        return true;
       
    }



    };


