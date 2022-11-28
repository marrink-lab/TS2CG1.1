#if !defined(AFX_Domain_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_Domain_H_9D4B21B8_C13C_5648_BF23_124095086234__INCLUDED_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include "LMatrix.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "point.h"

struct DomainLipid {
    std::string Name;           // name of the lipid
    double Ap;                  // AP of this lipid structure
    double Ratio;               // how many ratio in the domain
    int no_created;          // how many has been created
    int MaxNo;                  // how many of this type we should create
 //   int DynamicMaxNo;
    
} ;

class Domain
{
public:
    
	Domain(int domaintypeid, std::vector<point*>  point);
	virtual ~Domain();
    
            inline  std::vector<point*>                 GetDomainPoint()                const  {return m_point;}
            inline  std::vector<DomainLipid>            GetDomainLipids()               const  {return m_AllDomainLipids;}
            inline  std::vector<DomainLipid*>           GetpDomainLipids()               const  {return m_pAllDomainLipids;}
            inline  bool                                GetHealth()                     const  {return m_Health;}
            inline  int                                 GetDomainID()                   const  {return m_DomainTypeID;}
            inline  int                                 GetDomainTotalLipid()                   const  {return m_DomainTotalLipid;}


public:
    void AddADomainLipid(std::string name, double Ap, double Ratio);
    void Configure(bool);
private:
    
    
    

private:
    int m_DomainTypeID;
    std::vector<point*>  m_point;
    std::vector<DomainLipid>  m_AllDomainLipids;
    std::vector<DomainLipid*>  m_pAllDomainLipids;
    bool m_Health;
    int m_DomainTotalLipid;



};


#endif
