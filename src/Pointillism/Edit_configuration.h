#if !defined(AFX_Edit_configuration_H_8F4B21B8_D13C_9321_TT23_124095086224__INCLUDED_)
#define AFX_Edit_configuration_H_8F4B21B8_D13C_9321_TT23_124095086224__INCLUDED_
/*
 *  Copyright Weria Pezeshkian (weria.pezeshkian@gmail.com), 2020.
 * This is a class to call different functions: at the moment, only one exist
 */
#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "exclusion.h"

class Edit_configuration
{
public:
    
	Edit_configuration( std::vector <std::string> arg);
	 ~Edit_configuration();

private:

Vec3D *m_pBox;
    std::vector<vertex*>      m_pAllV;
    std::vector<vertex>       m_AllV;
    std::vector<triangle*>    m_pAllT;
    std::vector<links*>       m_pAllLinks;
    std::vector<links*>       m_pHalfLinks1;
    std::vector<links*>       m_pHalfLinks2;
    std::vector<inclusion*>   m_pInc;
    std::vector<exclusion*>   m_pExc;

    int  m_Iteration;  // how many time we should increase the number of trinagles each time, *4T, therefore N_T = N_T*4^m_Iteration
    Vec3D m_Zoom ;
    std::string m_Folder ;  // Folder to write output in
    double m_AP ;
    std::string m_Shape;
    bool m_FindnewBox;
    double m_minRoughness;
    bool m_smooth;
    int m_monolayer;
    void Rescaling(Vec3D zoom );   // rescale the position and the box
    void UpdateGeometry( );  // updates curvature, area etc of each triangle, vertex etc

    std::string m_MosAlType;
    void BackMapOneLayer(int layer , std::string file, double);
    void MakeFlatMonolayer(int layer , std::string file, double);
    bool check(std::string file);     // a function to check how the ts file looklike and do nothing
    void Minimize(std::string file);   // may not work well, needs optimization
    
    double  PPBCM_Cluster(double , std::vector <double>);
    bool FileExist (const std::string& name); // check if a file exist



};


#endif
