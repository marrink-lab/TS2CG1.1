  #if !defined(AFX_BackMap_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_BackMap_CPP_9T4A21B7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include <math.h>
#include "BackMap.h"
#include "GroFile.h"
#include "ReadDTSFolder.h"
#include "GenerateUnitCells.h"
#include "Def.h"
#include "PDBFile.h"
#include "GenDomains.h"


/*
 
 further extentions
 
 1). in random inc, we should prevent them to be close
 2.) angle and theta is not used yet
 3.) pattern based protein insertion
 
 
 
 */




BackMap::BackMap(Argument *pArgu)
{
    m_monolayer = false;

    srand (pArgu->GetSeed());
    std::cout<<"\n";
    std::cout<<"=========================================================================================================="<<"\n";
    std::cout<<"****************** "<< SoftWareName <<" ******************   "<<"\n";
    std::cout<<"******************  Version:  "<<SoftWareVersion<<" ****************** \n";
    std::cout<<"=========================================================================================================="<<"\n";

    Nfunction f;      // In this class there are some useful function and we can use it.


    
    //===================== OutPut file name declaration and finding input file names ========================================
    std::string gname = pArgu->GetGeneralOutputFilename();
    m_FinalOutputGroFileName =gname+".gro";
    std::string m_FinalTopologyFileName=gname+".top";
    std::string dtsfolder = pArgu->GetDTSFolder();
    std::string strfilename = pArgu->GetStructureFileName();
    m_ResID = 1;

    m_Renormalizedlipidratio = pArgu->GetRenorm();
    m_Iter = pArgu->GetIter();
    double RCutOff = pArgu->GetRCutOff();     /// That will be counted as a cutoff for the protein_lipid distance
    m_InclusionDirectionType = pArgu->GetInclusionDirectionType(); //Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually
    //std::cout<<"The inclusion direction type is: "<<m_InclusionDirectionType<<"\n";
    //==========================================================================================================

    GenerateMolType  MOLTYPE(pArgu);
    m_MoleculesType = MOLTYPE.GetMolType();

    if(MOLTYPE.GetHealth()==false)
    {
        std::cout<<"-> aborted! You are allowed to try one more time. Kidding, please do not :) \n";
        std::exit(0);
        
    }

    //================== Reading DTS folder to get the points ====================================================
    ReadDTSFolder ReadDTSFile(dtsfolder);

    std::vector<inclusion*>  m_pInc = ReadDTSFile.GetInclusion();
    std::vector<exclusion*>  m_pExc = ReadDTSFile.GetExclusion();

    m_point1 = ReadDTSFile.GetUpperPoints();
    m_point2 = ReadDTSFile.GetInnerPoints();
    

    
    if(m_point2.size()==0)
    m_monolayer = true;

    Vec3D *pBox= ReadDTSFile.GetBox();
    m_Box (0)=(*pBox)(0);     m_Box (1)=(*pBox)(1);     m_Box (2)=(*pBox)(2);
    m_pBox = pBox;

    
    
    //=============== make wall; Wall info and data
    Wall CWall = pArgu->GetWall();
    CWall.UpdateBox(pBox);
    CWall.CreateWall(m_point1,m_point2);
    std::vector<bead> WB = CWall.GetWallBead();
    std::vector<point*> WPoint = CWall.GetWallPoint();
    //==========  We write the wall at the end

    
    //********************** Finding the total area of the layers
    m_TotalAreaUp = 0.0;
    m_TotalAreaDown = 0.0;
    for ( std::vector<point*>::iterator it = m_point1.begin(); it != m_point1.end(); it++ )
        m_TotalAreaUp+= (*it)->GetArea();
    for ( std::vector<point*>::iterator it = m_point2.begin(); it != m_point2.end(); it++ )
        m_TotalAreaDown+= (*it)->GetArea();

    if(m_monolayer == false)
    std::cout<<"--> Note: the total upper monolayer area is "<<m_TotalAreaUp<<" and the total lower monolayer area is "<<m_TotalAreaDown<<"\n";
    else
    std::cout<<"--> Note: the total monolayer area is "<<m_TotalAreaUp<<" \n";

    

    
    //********* Lets exclude the points based on exclusion
    
    if(m_pExc.size()!=0)
    {
        std::cout<<" Note: we are excluding points based on exclusion, If it is slow, contact the developer \n";
        // m_pExc
        //m_point2
        
        for ( std::vector<exclusion*>::iterator it = m_pExc.begin(); it != m_pExc.end(); it++ )
        {
            int pointid=(*it)->GetPointID();
            if(pointid<0 || pointid>m_point1.size())
            std::cout<<"error 23456 \n";
            
            point *Up_p1=m_point1.at(pointid);
            Vec3D Pos = Up_p1->GetPos();
            Vec3D N = Up_p1->GetNormal();
            double R = (*it)->GetRadius();
            if(R!=0)
            Up_p1->UpdateArea(0);

            for ( std::vector<point*>::iterator it1 = m_point1.begin(); it1 != m_point1.end(); it1++ )
            {
                Vec3D Pos1 = (*it1)->GetPos();
                Vec3D DP = Pos1-Pos;
                double dist = DP.norm();
                Vec3D UnitDP =DP*(1/dist);
                double sinT = fabs((UnitDP*N).norm());
                double cosT = fabs(N.dot(UnitDP,N));

                if(dist*sinT<=R && dist*cosT<6)
                {
                    (*it1)->UpdateArea(0);

                }

                
            }
            for ( std::vector<point*>::iterator it1 = m_point2.begin(); it1 != m_point2.end(); it1++ )
            {
                Vec3D Pos1 = (*it1)->GetPos();
                Vec3D DP = Pos1-Pos;
                double dist = DP.norm();
                Vec3D UnitDP =DP*(1/dist);
                double sinT = fabs((UnitDP*N).norm());
                double cosT = fabs(N.dot(UnitDP,N));

                if(dist*sinT<=R && dist*cosT< 6)
                {
                    (*it1)->UpdateArea(0);
                    
                }
                
            }
            

        }

        
    }
    

    
    //**********************
    
    

   
    //*****************************
    //*********************************
    // read strc file to find protein

    if( FindProteinList(strfilename)==false)
        std::exit(0);
    


    if(m_pInc.size()==0)
    {
        CreateRandomInclusion();
        for (std::vector<inclusion>::iterator it2 = m_RandomInc.begin() ; it2 != m_RandomInc.end(); ++it2)
            m_pInc.push_back(&(*it2));
    }
    else 
    {

    	for ( std::map<int,ProteinList>::iterator it = m_TotalProteinList.begin(); it != m_TotalProteinList.end(); it++ )
    	{
        	int number = 0;
    		int id = (it->second).ID;
    		for ( std::vector<inclusion*>::iterator it1 = m_pInc.begin(); it1 != m_pInc.end(); it1++ )
    		{

			if(id==(*it1)->GetTypeID())
			{
				number++;

			}
    		}
        	(it->second).created = number;
	}

    }

    //== Place the inclusions

    for ( std::map<int,ProteinList>::iterator it1 = m_TotalProteinList.begin(); it1 != m_TotalProteinList.end(); it1++ )
    {

        int plistid=it1->first;

        for ( std::vector<inclusion*>::iterator it = m_pInc.begin(); it != m_pInc.end(); it++ )
        {

            int id = (*it)->GetTypeID();
            if(plistid==id)
            {

            std::string ptype=(m_TotalProteinList.at(id)).ProteinName;

            int pointid = (*it)->GetPointID();
            Vec3D  Dir =  (*it)->GetDirection();
            point *Up_p1=m_point1.at(pointid);
            Vec3D N = Up_p1->GetNormal();
            Vec3D Pos = Up_p1->GetPos();
            Vec3D T1 =   Up_p1->GetP1();
            Vec3D T2 =   Up_p1->GetP2();
                if (m_MoleculesType.count(ptype) == 0)
                    std::cout << "Error:-----> molecule name " <<ptype<<" does not exist in the attached gro files \n";
                
                
            GenProtein(m_MoleculesType.at(ptype), id, Pos, N, Dir, T1,T2);
                

            }
        }
        plistid++;
    }
    
    std::cout<<" We have placed the proteins, now is time to add the lipids \n";
    std::vector<bead*> tempropbeads;
    std::vector<bead> temprobeads;
    
//=== m_FinalBeads contain all the beads, and will recieve more. So, to avoid change in the pointer reference we make a copy for Unit Cell check
    for ( std::vector<bead>::iterator it = m_FinalBeads.begin(); it != m_FinalBeads.end(); it++ )
        temprobeads.push_back((*it));
    for ( std::vector<bead>::iterator it = temprobeads.begin(); it != temprobeads.end(); it++ )
        tempropbeads.push_back(&(*it));
    
        GenerateUnitCells GCNT(tempropbeads, pBox,RCutOff,1.0);
    
        GCNT.Generate();
    
    
    
    std::vector<point*>  p1;       // points that do not touch proteins
    std::vector<point*>  p2;      // should be deleted
    
    


    // Here, we try to remove the points that are covered by the proteins. Since the lipid will be placed with Pcc=A/Ap; setting A=0 make it removed
    for ( std::vector<point*>::iterator it = m_point1.begin(); it != m_point1.end(); it++ )
    {
        bool rem = false;
        Vec3D Pos1 = (*it)->GetPos();
        Vec3D N =   (*it)->GetNormal();
        Vec3D Pos2 = Pos1 - N*1.5;
        rem = GCNT.anythingaround(Pos1);  //121945
        if(rem==false)
        rem = GCNT.anythingaround(Pos2);
        if(rem==true)
        (*it)->UpdateArea(0);
        else
        p1.push_back(*it);
    }
    if(m_monolayer == false)
    {
    for ( std::vector<point*>::iterator it = m_point2.begin(); it != m_point2.end(); it++ )
    {
        bool rem = false;
        Vec3D Pos1 = (*it)->GetPos();
        Vec3D N =   (*it)->GetNormal();
        Vec3D Pos2 = Pos1 - N*1.5;
        rem = GCNT.anythingaround(Pos1);  //121945
        if(rem==false)
        rem = GCNT.anythingaround(Pos2);
        if(rem==true)
        (*it)->UpdateArea(0);
        else
        p2.push_back(*it);
    }
    }
    
    // Make all the domain containing different lipids
    GenDomains GENDOMAIN(strfilename,p1,p2,m_Renormalizedlipidratio);
    std::vector<Domain*> pAllDomain = GENDOMAIN.GetDomains();
    std::cout<<" Number of the domains defined in the input file  "<< pAllDomain.size()/2 <<"\n";
    int layer = 0;
        std::cout<<"***************************  we aim to generate  ********************** \n";
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {

        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
        std::cout <<"*     For domain with ID "<<(*it)->GetDomainID() <<" \n";
        for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
        {
            std::cout <<"*     "<<(*it2).MaxNo<<"  "<<(*it2).Name<<"     "<<std::endl ;

        }
        }
        if(layer%2!=0 && m_monolayer == false)
        std::cout <<"   in the upper monolayer \n";
        else if(layer%2==0 && m_monolayer == false)
        std::cout <<"   in the lower monolayer \n";
        if(layer%2!=0 && m_monolayer == true)
        std::cout <<"   in the  monolayer \n";
    }
    

//    for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
    {
        
    }
    
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
            std::vector<DomainLipid*> DL = (*it)->GetpDomainLipids();
       /*     std::vector<DomainLipid> DLL = (*it)->GetDomainLipids();

            int s=0;
            for ( std::vector<DomainLipid*>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                std::cout <<s<<" *   begining   "<<(*it2)->MaxNo<<"  "<<(*it2)->Name<<"     "<<std::endl ;
                s++;
            }
         s=0;
        for ( std::vector<DomainLipid>::iterator it2 = DLL.begin(); it2 != DLL.end(); it2++ )
        {
            std::cout <<s<<" *   beginingQ   "<<(*it2).MaxNo<<"  "<<(*it2).Name<<"     "<<std::endl ;
            s++;
        }*/
        int iteration = 0;
        int madetotallipids = 0;
        int lipidlistID = 0;
        int NoMadeLipid = 0 ;  // temperory because our strcture cannot increase it
        int CreateTotalLipids = (*it)->GetDomainTotalLipid();

        std::vector<point*>  dpoint = (*it)->GetDomainPoint();
        while(true)
        {
            iteration++;

                if( lipidlistID>DL.size()-1)
                    break;
                if(madetotallipids==CreateTotalLipids)
                    break;
                if(iteration>m_Iter*(dpoint.size()))
                {
                std::cout<<" Warning: With "<< m_Iter <<" iterations, we could not place the expected number of the lipids \n";
                std::cout<<" if you are unhappy, increase the number of the iteration with option -iter, or regenerate the points \n";
                break;
                }



            
            int pointid = rand()%(dpoint.size());
            point* tempoint = dpoint.at(pointid);

            DomainLipid *LL=DL.at(lipidlistID);
            int RNG=(rand()%CreateTotalLipids)+1;
            Vec3D  Dir(0,0,0);
            Vec3D N = tempoint->GetNormal();
            Vec3D T1 =   tempoint->GetP1();
            Vec3D T2 =   tempoint->GetP2();

            std::string ltype = LL->Name;
            Vec3D Pos = tempoint->GetPos();
            double area = tempoint->GetArea();
            double rn = double(rand()%(1000000))/1000000.0;
            double prob=area/(LL->Ap);
           // std::cout<<DL.size()<<"  "<<LL->Ap<<"We 2222get here \n";

            if(prob>rn && LL->MaxNo>NoMadeLipid)
            {
                madetotallipids++;
                NoMadeLipid++;
                int t = LL->no_created;
                ((*it)->GetpDomainLipids()).at(lipidlistID)->no_created=t+1;
                //std::cout<< ((*it)->GetpDomainLipids()).at(lipidlistID)->no_created<<" "<<t+1<<"no beads \n";
                if (m_MoleculesType.count(ltype) == 0)
                    std::cout << "Error:-----> molecule name " <<ltype<<" does not exist in the lib files \n";
                
                GenLipid(m_MoleculesType.at(ltype), 0, Pos, N, Dir, T1, T2);
                (dpoint.at(pointid))->UpdateArea(0);
            }
            if(LL->MaxNo==NoMadeLipid )
            {
           //  std::cout<<"domain id "<<(*it)->GetDomainID()<<"  name "<<LL->Name<<" created for the domain "<<madetotallipids<<"  domain "<<CreateTotalLipids<<"  "<<lipidlistID<<"  "<<NoMadeLipid<<"  "<<LL->MaxNo<<"\n";
                lipidlistID++;
                NoMadeLipid = 0;


            }
            
        }

        
    }
    
    std::cout<<"*************************** Number of created Lipids,   ********************** \n";
    layer = 0;
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
            std::cout <<"*     For domain with ID "<<(*it)->GetDomainID() <<" \n";
            for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                std::cout <<"*     "<<(*it2).no_created<<"  "<<(*it2).Name<<"     "<<std::endl ;
                
            }
        }
        if(layer%2!=0 && m_monolayer == false)
        std::cout <<"   in the upper monolayer \n";
        else if(layer%2==0 && m_monolayer == false)
        std::cout <<"   in the lower monolayer \n";
        if(layer%2!=0 && m_monolayer == true)
        std::cout <<"   in the  monolayer \n";
    }

    //=============== write the wall info
    if(WPoint.size()>0 && CWall.GetState()==true)
    {
        PDBFile pdb;
        std::string pdbfile = "wall.pdb";
        pdb.WritePDBFile(pdbfile, WPoint);
    }
    else if(CWall.GetState()==true)
    {
        std::cout<<"Note ----> No wall.pdb file will be generated since the total created wall beads are zero \n";
    }
    for (std::vector<bead>::iterator it = WB.begin() ; it != WB.end(); ++it)
    {
        m_FinalBeads.push_back((*it));
    }
    //============== End Wall info and data

    
    WriteFinalGroFile();
    
    
    //==========================================================================================================
    //=============== Open Topology files and make mols
    //==========================================================================================================
    std::ofstream Topgro;
    Topgro.open(m_FinalTopologyFileName.c_str());
    Topgro<<" ;This file was generated by TS Back Mapping \n";
    Topgro<<" [ system ] \n";
    Topgro<<" Expect a large membrane \n";
    Topgro<<" [ molecules ] \n";
    
    for ( std::map<int,ProteinList>::iterator it = m_TotalProteinList.begin(); it != m_TotalProteinList.end(); it++ )
        Topgro<<(it->second).ProteinName<<"   "<<(it->second).created<<"\n";

    layer = 0;
    
    
    for ( std::vector<Domain*>::iterator it = pAllDomain.begin(); it != pAllDomain.end(); it++ )
    {
        layer++;
        std::vector<DomainLipid> DL = (*it)->GetDomainLipids();
        
        if(layer%2!=0 || m_monolayer == false)
        {
            Topgro <<"; domain "<<(*it)->GetDomainID() <<" \n";
            if(layer%2!=0 && m_monolayer == false)
            Topgro  <<" ;  in the upper monolayer \n";
            else if(layer%2==0 && m_monolayer == false)
            Topgro <<" ;  in the lower monolayer \n";
            for ( std::vector<DomainLipid>::iterator it2 = DL.begin(); it2 != DL.end(); it2++ )
            {
                Topgro <<"     "<<(*it2).Name<<"  "<<(*it2).no_created<<"     "<<std::endl ;
                
            }
        }

    }
    if(WB.size()!=0)
    Topgro<<"Wall    "<<WB.size()<<"\n";



}
BackMap::~BackMap()
{
    
}
void BackMap::GenLipid(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
{
    //

    
    
     Tensor2 LG = TransferMatLG(Normal, t1, t2);
    
    
    
      std::vector<bead> vbeads = moltype.Beads;
            for ( std::vector<bead>::iterator it = vbeads.begin(); it != vbeads.end(); it++ )
            {
                Vec3D BPos((*it).GetXPos(),(*it).GetYPos(),(*it).GetZPos());
                Vec3D vX = LG*BPos+ Pos;
                int beadid = (m_FinalBeads.size()+1);
                bead TemB(beadid, (*it).GetBeadName(), (*it).GetBeadType(), (*it).GetResName(), m_ResID, vX(0), vX(1),vX(2));
                m_FinalBeads.push_back(TemB);
            }
    
    m_ResID++;
    
    



    
}
void BackMap::GenProtein(MolType moltype, int listid, Vec3D Pos, Vec3D Normal, Vec3D Dir,Vec3D t1,Vec3D t2)
{
    //
    

        

        Tensor2 LG = TransferMatLG(Normal, t1, t2);
        Tensor2 GL = LG.Transpose(LG);
    

        //===== to fit to the protein diretion
        Vec3D LocalDir;
    
        if(m_InclusionDirectionType=="Global") //Note: in the normal condition PLM always write global so only applicable if you want to change the point folder manually
        LocalDir = GL*Dir;
        else if(m_InclusionDirectionType=="Local")
        LocalDir = Dir;

    std::cout<<LocalDir(0)<<"  GOLB "<<LocalDir(1)<<"   "<<LocalDir(2)<<"  \n ";

        double C=LocalDir(0);
        double S= LocalDir(1);
        Tensor2 Rot=Rz(C,S);
        Vec3D DH= Normal*((m_TotalProteinList.at(listid)).Z0);
        double phi = (m_TotalProteinList.at(listid)).Phi;
        double theta = (m_TotalProteinList.at(listid)).Theta;

        //
        std::vector<bead> vbeads = moltype.Beads;
        for ( std::vector<bead>::iterator it = vbeads.begin(); it != vbeads.end(); it++ )
        {
            Vec3D BPos((*it).GetXPos(),(*it).GetYPos(),(*it).GetZPos());
            Vec3D vX = LG*(Rot*BPos)+ Pos+DH;
            int beadid = (m_FinalBeads.size()+1);
            bead TemB(beadid, (*it).GetBeadName(), (*it).GetBeadType(), (*it).GetResName(), m_ResID, vX(0), vX(1),vX(2));
            m_FinalBeads.push_back(TemB);
        }

        m_ResID++;
  
    
}
double BackMap::dist2between2Points(Vec3D X1,Vec3D X2)
{
    
    double dist2=0;
    
    double x1=X1(0);
    double y1=X1(1);
    double z1=X1(2);
    
    double x2=X2(0);
    double y2=X2(1);
    double z2=X2(2);
    
    
    double dx=x2-x1;
    double dy=y2-y1;
    double dz=z2-z1;
    
    if(fabs(dx)>(*m_pBox)(0)/2.0)
    {
        if(dx<0)
            dx=(*m_pBox)(0)+dx;
        else if(dx>0)
            dx=dx-(*m_pBox)(0);
    }
    if(fabs(dy)>(*m_pBox)(1)/2.0)
    {
        if(dy<0)
            dy=(*m_pBox)(1)+dy;
        else if(dy>0)
            dy=dy-(*m_pBox)(1);
    }
    if(fabs(dz)>(*m_pBox)(2)/2.0)
    {
        if(dz<0)
            dz=(*m_pBox)(2)+dz;
        else if(dz>0)
            dz=dz-(*m_pBox)(2);
    }

    dist2=dx*dx+dy*dy+dz*dz;
    return dist2;
}
void BackMap::WriteFinalGroFile()
{
    
   

    FILE *fgro;
    fgro = fopen(m_FinalOutputGroFileName.c_str(), "w");
    
    
    /// resid  res name   noatom   x   y   z
    const char* Title=" System ";
    int Size=m_FinalBeads.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<bead>::iterator it = m_FinalBeads.begin() ; it != m_FinalBeads.end(); ++it)
    {
        
        i++;
        double x=(*it).GetXPos();
        double y=(*it).GetYPos();
        double z=(*it).GetZPos();
        
        
        const char* A1=((*it).GetResName()).c_str();
        const char* A2=((*it).GetBeadName()).c_str();
        int resid=(*it).GetResid();
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%100000,A1,A2,i%100000,x,y,z );
        
    }
    
    
    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
    
    
}

Tensor2 BackMap::Rz(double cos, double sin)
{
    Tensor2  R;
    
    R(0,0) = cos;
    R(0,1) = -sin;
    R(1,0) = sin;
    R(1,1) = cos;
    R(2,2) = 1;

    
    
    return R;
}
void BackMap::CreateRandomInclusion()
{
    int totinc = 0;
    int totcreated = 0;
    
    if(m_TotalProteinList.size()!=0)
    {
        std::cout<<m_TotalProteinList.size()<<" According to the data and area we generate  proteins \n";
    }
    for ( std::map<int,ProteinList>::iterator it = m_TotalProteinList.begin(); it != m_TotalProteinList.end(); it++ )
    {
        (it->second).created = 0;
        double ratio = (it->second).Ratio;
        std::string type = (it->second).ProteinName;
        double area = (m_MoleculesType.at(type)).molarea;
        int neededno = m_TotalAreaUp/(area)*ratio;
        (it->second).Maxno = neededno;
        totinc+=neededno;
        
        std::cout<<" We will try to generate  "<<neededno<<"  "<<type<<" protein \n" ;

    }
    

    int id=0;
    int s=0;
    while (totcreated<totinc && s<(m_point1.size()))
    {
        s++;
        bool accept = true;
        int RNG1=(rand()%totinc)+1;
        int RNG2=(rand()%m_point1.size());
        int pointid = (m_point1.at(RNG2))->GetID();
        int tid=0;
        int l=0;
        bool chosen =false;
        int max = 0;
        double R = 0;
        
        
        ///======== here we only try to find one of the protein type randomly based on RNG1
        for ( std::map<int,ProteinList>::iterator it = m_TotalProteinList.begin(); it != m_TotalProteinList.end(); it++ )
        {
             max =  (it->second).Maxno;
            int nocreated = (it->second).created;
            //if(l<max && RNG1<=max)
            if(RNG1>l && RNG1<=max+l && nocreated<max)
            {
                tid = (it->second).ID;
                chosen = true;
                
                std::string type = (it->second).ProteinName;
                double area = (m_MoleculesType.at(type)).molarea;
                
                R=sqrt(area/acos(-1));
                
            }
            l=l+max;

        }

        
        
        //===================
        
        for ( std::vector<ExcludedVolumeBeads>::iterator it = m_ExcludeBeads.begin(); it != m_ExcludeBeads.end(); it++ )
        {
            
            Vec3D XP1 = it->X;
            double R1 = it->R;
            
            Vec3D XP2 =(m_point1.at(RNG2))->GetPos();
            
            if(XP2.dot((XP2-XP1),(XP2-XP1))<(R1+R)*(R1+R))
                accept = false;

            
        }
        
        //======================
        
        if(accept==true)
        {

            double d1= double(rand()%1000)/1000;
            double d2= double(rand()%1000)/1000;
            
            
            //
            Vec3D D(d1,d2,0);
            D=D*(1/(D.norm()));

            Vec3D N = (m_point1.at(RNG2))->GetNormal();
            Vec3D T1 =   (m_point1.at(RNG2))->GetP1();
            Vec3D T2 =   (m_point1.at(RNG2))->GetP2();
            Tensor2 LG = TransferMatLG(N,T1,T2);
            D=LG*D;
            //
            id++;
        inclusion inc(id, tid,pointid,D);
        m_RandomInc.push_back(inc);
            totcreated++;
            (m_TotalProteinList.at(tid)).created = (m_TotalProteinList.at(tid)).created +1;
            
            
            ExcludedVolumeBeads Ex;
            Ex.X =(m_point1.at(RNG2))->GetPos();
            Ex.R = R;
            m_ExcludeBeads.push_back(Ex);
            
            
        }
        

    }
    
}
/*void BackMap::CreateWallBead(std::vector<point*>  p1, std::vector<point*>  p2)
{
    std::string file="Wall.gro";
    FILE *fgro;
    fgro = fopen(file.c_str(), "w");
    
    
    /// resid  res name   noatom   x   y   z
    const char* Title=" Wall of Water Beads ";
    int Size=p1.size()+p2.size();
    
    fprintf(fgro,  "%s\n",Title);
    fprintf(fgro, "%5d\n",Size);
    int i=0;
    for (std::vector<point*>::iterator it = p1.begin() ; it != p1.end(); ++it)
    {
        
        i++;
        Vec3D X=(*it)->GetPos();
        Vec3D N=(*it)->GetNormal();
        
        X=X+N*(0.7);
        const char* A1="W2";
        const char* A2="PO4";
        int resid=1;
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%10000,A1,A2,i%10000,X(0),X(1),X(2) );
        
    }
    for (std::vector<point*>::iterator it = p2.begin() ; it != p2.end(); ++it)
    {
        
        i++;
        Vec3D X=(*it)->GetPos();
        Vec3D N=(*it)->GetNormal();
        
        X=X+N*(0.7);
        const char* A1="W2";
        const char* A2="PO4";
        int resid=1;
        fprintf(fgro, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",resid%10000,A1,A2,i%10000,X(0),X(1),X(2) );
        
    }
    
    fprintf(fgro,  "%10.5f%10.5f%10.5f\n",m_Box(0),m_Box(1),m_Box(2) );
    fclose(fgro);
    
}
 */
bool BackMap::FindProteinList(std::string filename)
{
    bool OK=true;
    Nfunction f;
    std::ifstream strfile;
    strfile.open(filename.c_str());
    std::string str,pname;
    double ap;
    int pid;
    
    
    
    int l=0;
    bool lipidflag=false;
    bool proteinflag=false;
    bool flag = false;
    while (true)
    {
       
        l++;
        std::getline (strfile,str);
        
        if(strfile.eof())
            break;
        
        std::vector<std::string> Line = f.split(str);
        //****************
        if(Line.size()!=0 && (Line.at(0)).at(0)!=';')
        {
            
            if((Line.at(0)).at(0)=='[' && flag==false)
            {
                str = f.trim(str);
                str.erase(std::remove(str.begin(), str.end(), '['), str.end());
                str.erase(std::remove(str.begin(), str.end(), ']'), str.end());
                str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
                flag=true;
                

                 if(str=="ProteinList")
                {
                    lipidflag = false;
                    proteinflag =true;
                }
                else if(str=="LipidsList")
                {
                    lipidflag = true;
                    proteinflag =false;
                }
                else if(str=="ShapeData")
                {
                    lipidflag = false;
                    proteinflag =false;
                }
                else
                {

                    std::cout<<" Error7787: line "<<l<<" <"<<str<<"> of the "<<filename<< " file \n";
                    OK=false;
                }
                
            }
            else if((Line.at(0))=="End" && flag==true)
            {
                flag=false;
            }
            else if(flag==true && proteinflag==true)
            {
                if(Line.size()>=6)
                {
                    ProteinList CP;
                    CP.ID = f.String_to_Int(Line.at(1));
                    CP.ProteinName = Line.at(0);
                    CP.Ratio = f.String_to_Double(Line.at(2));
                    CP.Phi=f.String_to_Double(Line.at(3));
                    CP.Theta=f.String_to_Double(Line.at(4));
                    CP.Z0=f.String_to_Double(Line.at(5));
                    m_TotalProteinList.insert(std::pair<int,ProteinList>(CP.ID, CP));
                }
                else
                {

                    std::cout<<" Error9837: line "<<l<<" <"<<str<<"> of the "<<filename<< " file \n";
                    OK=false;
                    
                }
                
            }
            else if(flag==true && lipidflag==true)
            {


            }

            
        }
        
    }
    strfile.close();

    
    
    
    
  if(m_TotalProteinList.size()!=0)
  {
    std::cout<<"*************************** Protein List ID ********************** \n";
    for ( std::map<int,ProteinList>::iterator it = m_TotalProteinList.begin(); it != m_TotalProteinList.end(); it++ )
    {
        std::cout <<"*    "<< it->first  <<" ---> "<< (it->second).ProteinName<< std::endl ;
    }
    std::cout<<"************************************************************** \n";
  }
    
    return OK;
}
Tensor2  BackMap::TransferMatLG(Vec3D Normal, Vec3D t1, Vec3D t2)
{
    
    
   // if(Normal.dot((t1*t2),Normal)<0)
    //t2=t2*(-1);

    
    Tensor2  GL(t1,t2,Normal);
    Tensor2 LG=GL.Transpose(GL);
    
    
    
    
    //======== test, ignore it
    /* Vec3D V1(1,2,4);
    Vec3D V2(-1,3,5);
    
    if(fabs((LG*(V1*V2))(0)-((LG*V1)*(LG*V2))(0))>0.01)
    {
    std::cout<<" Before "<<(LG*(V1*V2))(0)<<"   "<<(LG*(V1*V2))(1)<<"   "<<(LG*(V1*V2))(2)<<"   \n";
    std::cout<<" Before "<<((LG*V1)*(LG*V2))(0)<<"   "<<((LG*V1)*(LG*V2))(1)<<"   "<<((LG*V1)*(LG*V2))(2)<<"   \n";

    std::cout<<" After "<<((GL*t1))(0)<<"   "<<((GL*t1))(1)<<"   "<<((GL*t1))(2)<<"   \n";
    
    std::cout<<" Before "<<V1.dot(V1,V2)<<"  after dot "<<V1.dot(LG*V1,LG*V2)<<"   "<<"   \n";
    
    
    std::cout<<"  "<<(t1*t2)(0)<<"   "<<(t1*t2)(1)<<"   "<<(t1*t2)(2)<<"   \n";
    std::cout<<"  "<<(Normal)(0)<<"   "<<(Normal)(1)<<"   "<<(Normal)(2)<<"   \n";

    std::cout<<" ==================  \n";
    
    }
   */
    
    
    return LG;
    
}







#endif



