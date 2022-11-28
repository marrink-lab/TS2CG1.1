


#include "Argument.h"
#include "Job.h"
#include "Nfunction.h"
#include "BackMap.h"
#include "Shape_1DSinBuilder.h"
#include "Membrane_Builder.h"

Job::Job(std::vector <std::string> argument)
{



Argument a(argument);
    std::string exacutable=argument.at(0);
    int  n=exacutable.size();
    char L1 = exacutable.at(n-1);
    char L2 = exacutable.at(n-2);
    char L3 = exacutable.at(n-3);
    
    std::string function = a.GetFunction();
    bool condition = true;
if(a.GetHealth()==true )
{
    if (n>3)
    {
        char L4 = exacutable.at(n-4);
        if (L4!='/')
        {
            condition = false;
            std::cout<<argument.at(0)<<" 1. executable is unknown \n";
        }
        
    }
    if(condition == true)
    {
        if( L3 == 'P' && L2 == 'C' && L1 == 'G')
        {
            if(function=="backmap")
            BackMap B(&a);
            else if(function=="1dsin")
            Shape_1DSinBuilder B(&a);
            else if(function=="analytical_shape")
                Membrane_Builder B(&a);
            else
            std::cout<<function<<" function is not regognized \n";

            
        }
        else if(n<3)
        {
            condition  = false;
            std::cout<<argument.at(0)<<"2. executable is unknownt \n";
        }
        
    }

}


}
Job::~Job()
{
    
}

