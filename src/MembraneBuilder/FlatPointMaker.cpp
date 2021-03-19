  #if !defined(AFX_FlatPointMaker_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_FlatPointMaker_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "FlatPointMaker.h"
#include "GroFile.h"
#include "ReadDTSFolder.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"
FlatPointMaker::FlatPointMaker(Argument *pArgu)
{
    m_WallBox.push_back(0);
    m_WallBox.push_back(1);
    m_WallBox.push_back(0);
    m_WallBox.push_back(1);

    m_Box(0) = 10;m_Box(1) = 10;m_Box(2) = 10;
    m_pBox = &m_Box;
    m_Thickness = 4;
    m_Density = 0.4;
    m_WallDensity = 0.1;
    std::string ifilename = pArgu->GetStructureFileName();
    Initialize(ifilename);


    //=============== We need the H wall distance, and the only info needed about the wall here
    Wall CWall = pArgu->GetWall();
    double Hwall = CWall.GetH();

    //*********

    m_Point1 = CalculateArea_MakePoints(1, m_Density,m_Thickness/2,false);
    m_Point2 = CalculateArea_MakePoints(-1, m_Density,m_Thickness/2,false);

    
    m_WallPoint1 = CalculateArea_MakePoints(1, m_WallDensity,m_Thickness/2+Hwall,true);
    m_WallPoint2 = CalculateArea_MakePoints(-1, m_WallDensity,m_Thickness/2+Hwall,true);

}
FlatPointMaker::~FlatPointMaker()
{
    
}
std::vector<point> FlatPointMaker::CalculateArea_MakePoints(int layer, double APL,double H, bool wall)
{
    std::vector<point> Cpoints;
    double dt = 0.0001/m_Box(1);
    int Nup = int((*m_pBox)(0)/dt+H/(2*dt))+1;
    int Ndown = int(H/(2*dt))+1;
    std::vector <double> c;  // we set curvature to zero as it is not important for this function
    c.push_back(0);
    c.push_back(0);
    std::vector<Vec3D> Pos,NormalV,T1vec,T2vec;
    std::vector<Vec3D> V1, V2,V3,V4;
    Vec3D OLDV1;
    Vec3D t2(0,1,0);

    double xmin = 0;
    double ymin = 0;
    double xmax = (*m_pBox)(0);
    double ymax = (*m_pBox)(1);
    
    if(wall==true)
    {
        xmin = (m_WallBox.at(0))*((*m_pBox)(0));
        xmax = (m_WallBox.at(1))*((*m_pBox)(0));
        ymin = (m_WallBox.at(2))*((*m_pBox)(1));
        ymax = (m_WallBox.at(3))*((*m_pBox)(1));


    }
    
    for (int i=-Ndown;i<Nup;i++)
    {
        double t = i*dt;
        Vec3D X = F(t,layer,H);
        
        
        if(X(0)>=0 && X(0)<(*m_pBox)(0))
        {
            V1.push_back(X);
            Vec3D T = X-OLDV1;
            T = T*(1/T.norm());
            V3.push_back(T);
            V4.push_back(t2);
            OLDV1 = X;
            Vec3D n = T*t2;
            n=n*(1/n.norm());
            if(layer==-1)
                n=n*(-1);
            V2.push_back(n);
            
        }
        
    }
    double area = 0;
    for (int i=1;i<V1.size();i++)
    {
        double DX = (V1.at(i-1))(0)-(V1.at(i))(0);
        double DZ = (V1.at(i-1))(2)-(V1.at(i))(2);
        area+=sqrt(DX*DX+DZ*DZ);
    }
    double Lenght = (*m_pBox)(0);
    area = Lenght*(*m_pBox)(1);
    double DL = sqrt(APL);
    Vec3D DXOLD = V1.at(0);
    Pos.push_back(V1.at(0));
    NormalV.push_back(V2.at(0));
    T1vec.push_back(V3.at(0));
    T2vec.push_back(V4.at(0));
    for(int i=1;i<V1.size();i++)
    {
        Vec3D DX = DXOLD-V1.at(i);
        if(DX.norm()>=DL)
        {
            DXOLD = V1.at(i);
            
            Pos.push_back(V1.at(i));
            NormalV.push_back(V2.at(i));
            T1vec.push_back(V3.at(i));
            T2vec.push_back(V4.at(i));
        }
    }
    int beadid = 0;
    double dy = sqrt(APL);
    int NY = m_Box(1)/dy;
    APL = area/double(Pos.size()*NY);
    for (int i=0;i<Pos.size();i++)
    {
        for (int j=0;j<NY;j++)
        {
                double y=(double(j))*dy+0.5*double(i%2)*dy;
                (Pos.at(i))(1) = y;
               if((Pos.at(i))(0)>=xmin && (Pos.at(i))(0)<=xmax && y>=ymin && y<=ymax)
               {
                point p(beadid, APL, Pos.at(i), NormalV.at(i), T1vec.at(i), T2vec.at(i) , c );
                Cpoints.push_back(p);
                beadid++;
               }
        }
    }
    for (std::vector<point >::iterator it = Cpoints.begin() ; it != Cpoints.end(); ++it)
    {
            it->UpdateArea(area/double(Cpoints.size()));
    }

    
    return Cpoints;

}
Vec3D FlatPointMaker::F(double t, int layer,double H)
{
    double pi = acos(-1);
    Vec3D F;
    F(0) =t;
    F(2) =(*m_pBox)(2)/2;
    F = F + Normal(t)*(H*double(layer)/(Normal(t)).norm());
    return F;
}
Vec3D FlatPointMaker::Normal(double t)/// normal function to the mid surface at point denotaed by t
{
    Vec3D N(0,0,1);

    return N;
}
void FlatPointMaker::Initialize(std::string filename)
{
    
    bool OK=true;
    Nfunction f;
    std::ifstream file;
    file.open(filename.c_str());
    bool flag = false;
    std::string str;

    while (true)
    {
        std::getline (file,str);
        if(file.eof())
            break;
        
        std::vector<std::string> Line = f.split(str);
        if(Line.size()!=0 && (Line.at(0)).at(0)!=';')
        {
            if((Line.at(0)).at(0)=='[' && flag==false)
            {
                str = f.trim(str);
                str.erase(std::remove(str.begin(), str.end(), '['), str.end());
                str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
                str.erase(std::remove(str.begin(), str.end(), ' '), str.end());

                if(str=="ShapeData")
                    flag = true;
            }
            else if((Line.at(0))=="End" && flag==true)
            {
                flag=false;
            }
            else if(flag==true)
            {
                if(Line.at(0)=="Box")
                {
                    if(Line.size()<4)
                        std::cout<<" Error: Box information in the str file is not correct \n";
                    else
                    {
                        m_Box(0) = f.String_to_Double(Line.at(1));
                        m_Box(1) = f.String_to_Double(Line.at(2));
                        m_Box(2) = f.String_to_Double(Line.at(3));

                    }
                }
                else if(Line.at(0)=="ShapeType")
                {
                }
                else if(Line.at(0)=="WallRange")
                {
                    if(Line.size()<5)
                        std::cout<<" Error: Wall range information in the str file is not correct \n";
                    else
                    {
                        m_WallBox.at(0) = f.String_to_Double(Line.at(1));
                        m_WallBox.at(1) = f.String_to_Double(Line.at(2));
                        m_WallBox.at(2) = f.String_to_Double(Line.at(3));
                        m_WallBox.at(3) = f.String_to_Double(Line.at(4));

                    }
                }
                else if(Line.at(0)=="Density")
                {
                    if(Line.size()<3)
                        std::cout<<" Error: Density information in the str file is not correct \n";
                    else
                    {
                        m_Density = f.String_to_Double(Line.at(1));
                        m_WallDensity = f.String_to_Double(Line.at(2));
                        
                    }
                }
                else if(Line.at(0)=="Thickness")
                {
                    if(Line.size()<2)
                        std::cout<<" Error: Thickness information in the str file is not correct \n";
                    else
                    {
                        m_Thickness = f.String_to_Double(Line.at(1));
                        
                    }
                }
            }


        }
            
            
        }
        
    
    file.close();

    
}






#endif



