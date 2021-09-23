/*
 <TOVTK4>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 
 <TOVTK5>  Modified by Dr Xin Lai. (see https://github.com/lxidea/ToVTK/)

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

#include "Functions.h"
#include "JException.h"
#include "JCfgRun.h"
#include "TypesDef.h"
#include "JPartDataBi4.h"
#include "JVtkLib.h"
#include "JDataArrays.h"
#include "JOutputCsv.h"
#include "JXml.h"
#include "JSpaceCtes.h"
#include "JSpaceEParms.h"
#include "JSpaceParts.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <filesystem>
#include "omp.h"

//using namespace std;
using std::string;
using std::exception;

typedef struct
{
    unsigned i,j,k,l,m,n,o,p;
}conn;
const unsigned maxlist = 270;
const char *APP_NAME="ToVtk v5.0.028 (03-09-2021)";

//==============================================================================
// Invoca una excepcion referente a gestion de ficheros
//==============================================================================
void ExceptionFile(string msg,string file){
  msg=msg+"\nFile:["+file+"]\n";
  throw msg;
}

//==============================================================================
// Invoca una excepcion generica
//==============================================================================
void ExceptionText(string msg){
  throw msg;
}

void RunFilesGrid(const JCfgRun *cfg){
  const bool outtype=true,outmk=true,outpres=true;
  const bool outdp=false;

  //-Reads XML file.
  JSpaceCtes* xmlctes=NULL;
  JSpaceEParms* xmleparms=NULL;
  JSpaceParts* xmlparts=NULL;
  if(outmk){
    if(!fun::FileExists(cfg->FileXml))ExceptionFile("XML file of case was not found.",(cfg->FileXml.empty()? string("???"): cfg->FileXml));
    JXml xml; xml.LoadFile(cfg->FileXml);
    xmlctes=new JSpaceCtes();     xmlctes->LoadXmlRun(&xml,"case.execution.constants");
    xmleparms=new JSpaceEParms(); xmleparms->LoadXml(&xml,"case.execution.parameters");
    xmlparts=new JSpaceParts();   xmlparts->LoadXml(&xml,"case.execution.particles");
  }

  //-Prepares input files.
  const bool onefile=!cfg->FileIn.empty();
  const string casein=cfg->FileIn;
  const string dirin=fun::GetDirWithSlash(cfg->DirIn);
  const int last=cfg->Last;
  int part=(onefile||cfg->First<0? 0: cfg->First);
  bool cont = cfg->Continue;
  string file=fun::FileNameSec(cfg->SaveGrid,part);
  file=fun::AddExtension(file,"vtk");
  while (cont && fun::FileExists(file)) {
    part++;
    file=fun::FileNameSec(cfg->SaveGrid,part);
    file=fun::AddExtension(file,"vtk");
  }

  unsigned npiece=0;
  unsigned casenp=0,casenbound=0,casenfixed=0,casenmoving=0,casenfloat=0,casenfluid=0;
  unsigned casenpmax=0;
  double cteb=0,rhop0=0,gamma=0;

  // for grid rebuild
  unsigned gnp=0,gconn=0;
  tfloat3 *gpos=NULL,*gvel=NULL;
  float *gvol=NULL,*gmass=NULL;
  unsigned *gidp=NULL;
  float *grhop=NULL,*gpres=NULL;
  byte *gtype=NULL,*gmk=NULL;
  float *gdropr=NULL,*ggrav=NULL;
  unsigned **neiblist=NULL;
  unsigned short *neibNum=NULL;
  tuint8* connectivity=NULL;

  double timestep=0;
  unsigned *idp=NULL;
  tfloat3 *pos=NULL,*vel=NULL;
  tdouble3 *posd=NULL;
  float *rhop=NULL,*vol=NULL;
  unsigned *ridp=NULL;
  byte *type=NULL,*mk=NULL;
  float *mass=NULL;
  float *pres=NULL;
  float *dropr=NULL,*grav=NULL;

  //-Processes input files.
  byte npie=0;
  file=JPartDataBi4::GetFileData(casein,dirin,part,npie);
  if(file.empty())ExceptionText("Error: Data files not found.");
  bool firstdata=true;
  while((last<0||part<=last) && fun::FileExists(file)){
    if(!npiece){//-Reads initial data.
      JPartDataBi4 pd;
      if(onefile)pd.LoadFileCase("",casein,0,npie);
      else pd.LoadFilePart(dirin,part,0,npie);
      npiece=pd.GetNpiece();
      //maxx=ToTFloat3(pd.Get_DomainMax());
      //minx=ToTFloat3(pd.Get_DomainMin());
      if(npiece>1)ExceptionText("Error: The number of pieces is higher than 1.");
      casenp=std::max((unsigned)pd.Get_CaseNp(),(unsigned)pd.Get_Npok());
      casenpmax=(unsigned)pd.Get_NpTotal();
      //if(pd.Get_CaseNp()!=casenp)ExceptionText("Error: The number of particles is too big.");
      if(!pd.Get_IdpSimple())ExceptionText("Error: Only Idp (32 bits) is valid at the moment.");
      casenfluid=(unsigned)pd.Get_CaseNfluid();
      casenfixed=(unsigned)pd.Get_CaseNfixed();
      casenmoving=(unsigned)pd.Get_CaseNmoving();
      casenfloat=(unsigned)pd.Get_CaseNfloat();
      casenbound=casenp-casenfluid;
      cteb=pd.Get_B();
      rhop0=pd.Get_Rhop0();
      gamma=pd.Get_Gamma();
      pos=new tfloat3[casenpmax];
      posd=new tdouble3[casenpmax];
      vel=new tfloat3[casenpmax];
      vol=new float[casenpmax];
      rhop=new float[casenpmax];
      idp=new unsigned[casenpmax];
      mass=new float[casenpmax];
      //-Allocates memory to other variables.
      if(outmk)mk=new byte[casenpmax];
      if(outtype)type=new byte[casenpmax];
      if(outpres)pres=new float[casenpmax];
      if(outmk)ridp=new unsigned[pd.Get_NpTotal()];
      if(outdp) {dropr=new float[casenpmax];
                 grav=new float[casenpmax];}
      //printf("CaseNp:%u,Npok:%u,NpTotal:%u\n",pd.Get_CaseNp(),pd.Get_Npok(),pd.Get_NpTotal());
    }

    //-Reads particle data.
    unsigned np=0;
    double dp;
    tfloat3 minx={FLT_MAX,FLT_MAX,FLT_MAX};
    tfloat3 maxx={FLT_MIN,FLT_MIN,FLT_MIN};
    tfloat3 posmin,posmax;
    for(unsigned cp=0;cp<npiece && cp<1;cp++){
      file=(onefile? JPartDataBi4::GetFileNameCase(casein,cp,npiece): dirin+JPartDataBi4::GetFileNamePart(part,cp,npiece));
      if(!cp)printf("load> %s\n",file.c_str());
      JPartDataBi4 pd;
      if(onefile)pd.LoadFileCase("",casein,0,npiece);
      else pd.LoadFilePart(dirin,part,0,npiece);
      if(!cp)timestep=pd.Get_TimeStep();
      posmin=ToTFloat3(pd.Get_MapPosMin());
      posmax=ToTFloat3(pd.Get_MapPosMax());
      const bool possimple=pd.Get_PosSimple();
      const unsigned npok=pd.Get_Npok();
      const bool data2d=pd.Get_Data2d();
      dp=pd.Get_Dp();
      if(npok){
        // test arrays
        //printf("Arrays Count:%d\n",pd.ArraysCount());
        //for(unsigned i=0;i<pd.ArraysCount();i++)
        //  printf("%d:%s\n",pd.ArrayName(i));
        //-Loads data from PART.
        pd.Get_Idp(npok,idp);
        pd.Get_Vel(npok,vel);
        //pd.Get_Vol(npok,vol);
        pd.Get_Rhop(npok,rhop);
        const double rhop0 = pd.Get_Rhop0();
        
        //pd.Get_Mass(npok,mass);
        for(unsigned p=0;p<npok;p++) {
          vol[p]=(data2d?dp*dp:dp*dp*dp)*rhop0/rhop[p];
          mass[p]=(data2d?dp*dp:dp*dp*dp)*rhop0;
        }
        if(possimple){
          pd.Get_Pos(npok,pos);
        }
        else{ 
          pd.Get_Posd(npok,posd);
          for(unsigned p=0;p<npok;p++)pos[p]=ToTFloat3(posd[p]);
        }
        for(unsigned i=0;i<npok;i++){
          maxx.x=std::max(maxx.x,pos[i].x);
          maxx.y=std::max(maxx.y,pos[i].y);
          maxx.z=std::max(maxx.z,pos[i].z);
          minx.x=std::min(minx.x,pos[i].x);
          minx.y=std::min(minx.y,pos[i].y);
          minx.z=std::min(minx.z,pos[i].z);
        }
        if(std::isinf(maxx.x) || std::isinf(maxx.y) || std::isinf(maxx.z) || std::isinf(minx.x) || std::isinf(minx.y) || std::isinf(minx.z)){
          maxx = posmax;
          minx = posmin;
        }
        //printf("domain dimensions:\nX:%f - %f\nY:%f - %f\nZ:%f - %f\n",minx.x,maxx.x,minx.y,maxx.y,minx.z,maxx.z);
        if(outdp){
          pd.Get_Dropr(npok,dropr);
          pd.Get_Grav(npok,grav);
        }
      }
      np=npok;
    }

    //-Loads other vars.
    if(pres)for(unsigned c=0;c<np;c++)pres[c]=(float)(cteb*(pow(rhop[c]/rhop0,gamma)-1.));
    if(type){
      for(unsigned p=0;p<np;p++){
        const unsigned id=idp[p];
        type[p]=(id<casenfixed? 0: (id<casenmoving? 1: (id<casenfloat? 2: 3)));
      }
    }
    if(ridp)for(unsigned c=0;c<np;c++)ridp[idp[c]]=c;
    if(mk){
      for(unsigned c=0;c<xmlparts->CountBlocks();c++){
        const JSpacePartBlock &block=xmlparts->GetBlock(c);
        const byte mkblock=(byte)block.GetMk();
        const unsigned ipend=block.GetBegin()+block.GetCount();
        for(unsigned ip=block.GetBegin();ip<ipend;ip++)mk[ridp[ip]]=mkblock;
      }
    }

    //printf("Allocate Memory for grid building\n");
    // generate grid and interpolate
    const double geomx = maxx.x - minx.x;
    const double geomy = maxx.y - minx.y;
    const double geomz = maxx.z - minx.z;

    const unsigned nx = geomx/dp + 1;
    const unsigned ny = geomy/dp + 1;
    const unsigned nz = geomz/dp + 1;
    gnp = nx*ny*nz;

    gpos = new tfloat3[gnp];
    gidp = new unsigned[gnp];
    gvel = new tfloat3[gnp];
    gvol = new float[gnp];
    grhop= new float[gnp];
    gpres= new float[gnp];
    gtype= new byte[gnp];
    gmass= new float[gnp];
    gmk  = new byte[gnp];
    gdropr=new float[gnp];
    ggrav= new float[gnp];

    neiblist=new unsigned*[gnp];
    neibNum=new unsigned short[gnp];

    //printf("Initializing grid\n");
    unsigned idx=0;
    tfloat3 mingx={FLT_MAX,FLT_MAX,FLT_MAX};
    tfloat3 maxgx={FLT_MIN,FLT_MIN,FLT_MIN};
    for (auto i = 0; i < nx; i++)
      for (auto j = 0; j < ny; j++)
        for (auto k = 0; k < nz; k++, idx++)
        {
          tdouble3 _posd =
          {minx.x + i*dp,minx.y + j*dp,minx.z + k*dp};
          tfloat3 _pos = ToTFloat3(_posd);
          gpos[idx] = _pos;
          tfloat3 _vel = {0.0f,0.0f,0.0f};
          gidp[idx] = idx + 1;
          gvel[idx] = _vel;
          grhop[idx] = 0;
          gpres[idx] = 0;
          gmk[idx] = 0;
          gmass[idx] = 0;
          gvol[idx] = 0;
          gtype[idx] = 0;
          neibNum[idx] = 0;
          neiblist[idx] = new unsigned[maxlist];
          //for(unsigned i=0;i<maxlist;i++) neiblist[idx][i]=0;
          memset(neiblist[idx],0,maxlist*sizeof(unsigned));
          maxgx.x = std::max(maxgx.x,_pos.x);
          maxgx.y = std::max(maxgx.y,_pos.y);
          maxgx.z = std::max(maxgx.z,_pos.z);
          mingx.x = std::min(mingx.x,_pos.x);
          mingx.y = std::min(mingx.y,_pos.y);
          mingx.z = std::min(mingx.z,_pos.z);
        }
        
    gconn = (nx-1)*(ny-1)*(nz-1);
    connectivity = new tuint8[gconn];
    idx = 0;
    
    for (auto ii = 0; ii < nx-1; ii++)
      for (auto jj = 0; jj < ny-1; jj++)
        for (auto kk = 0; kk < nz-1; kk++,idx++)
        {
          const unsigned i=ii*ny*nz+jj*nz+kk;
          const unsigned j=(ii+1)*ny*nz+jj*nz+kk;
          const unsigned k=(ii+1)*ny*nz+(jj+1)*nz+kk;
          const unsigned l=ii*ny*nz+(jj+1)*nz+kk;
          const unsigned m=ii*ny*nz+jj*nz+kk+1;
          const unsigned n=(ii+1)*ny*nz+jj*nz+kk+1;
          const unsigned o=(ii+1)*ny*nz+(jj+1)*nz+kk+1;
          const unsigned p=ii*ny*nz+(jj+1)*nz+kk+1;
          connectivity[idx] = {i,j,k,l,m,n,o,p};
        }
    
    //printf("Grid node: %d x %d x %d , %d\n",nx,ny,nz,gnp);
    //printf("Grid element: %d x %d x %d, %d\n", nx-1,ny-1,nz-1,gconn);
    
    //printf("building background grid nodes neighbor list linked\n");
    //printf("running in %d threads\n",omp_get_max_threads());

    const unsigned ngnp = np + gnp;
    const double dgeomx = maxx.x - minx.x;
    const double dgeomy = maxx.y - minx.y;
    const double dgeomz = maxx.z - minx.z;
    const double ggeomx = maxgx.x - mingx.x;
    const double ggeomy = maxgx.y - mingx.y;
    const double ggeomz = maxgx.z - mingx.z;
    const double minxx = std::min(minx.x,mingx.x);
    const double minyy = std::min(minx.y,mingx.y);
    const double minzz = std::min(minx.z,mingx.z);
    const double dggeomx = std::max(dgeomx,ggeomx);
    const double dggeomy = std::max(dgeomy,ggeomy);
    const double dggeomz = std::max(dgeomz,ggeomz);

    const double hh = 3.0*dp;
    int ngridx = std::min(int(dggeomx / hh) + 1, 1000);
    int ngridy = std::min(int(dggeomy / hh) + 1, 1000);
    int ngridz = std::min(int(dggeomz / hh) + 1, 1000);

    int ***grid = new int**[ngridx];
    for (auto i = 0; i < ngridx; i++)
    {
        grid[i] = new int*[ngridy];
        for (auto j = 0; j < ngridy; j++)
        {
            grid[i][j] = new int[ngridz];
            for (auto k = 0; k < ngridz; k++)
            {
                grid[i][j][k] = 0;
            }
            
        }
        
    }

    //----------------
    int* xgcell = new int[ngnp];
    int* ygcell = new int[ngnp];
    int* zgcell = new int[ngnp];
    int*celldata= new int[ngnp];
    bool*cvalid = new bool[ngnp];

    for (auto i = 0; i < ngnp; i++)
    {
        xgcell[i] = ygcell[i] = zgcell[i] = celldata[i] = 0;
        cvalid[i] = true;
    }

    for (int i = 1; i <= ngnp; i++)
    {
        const double xg = i<=np?pos[i-1].x:gpos[i-np-1].x;
        const double yg = i<=np?pos[i-1].y:gpos[i-np-1].y;
        const double zg = i<=np?pos[i-1].z:gpos[i-np-1].z;
        cvalid[i-1] = !(std::isinf(xg) || std::isinf(yg) || std::isinf(zg) ||
                std::isnan(xg) || std::isnan(yg) || std::isnan(zg));
        
        if (xg>posmax.x || xg<posmin.x || yg>posmax.y || yg<posmin.y || zg>posmax.z || zg<posmin.z) cvalid[i-1]=false;

        if (!cvalid[i-1]) continue;

        const int xxcell = int((xg - minxx) / hh + 1.0);
        const int yycell = int((yg - minyy) / hh + 1.0);
        const int zzcell = int((zg - minzz) / hh + 1.0);

        xgcell[i-1] = xxcell;
        ygcell[i-1] = yycell;
        zgcell[i-1] = zzcell;
        celldata[i-1] = grid[xxcell-1][yycell-1][zzcell-1];
        grid[xxcell-1][yycell-1][zzcell-1]=i;
    }

    clock_t starttime = clock();
    clock_t testtime = starttime;

    unsigned irun = 0;
    unsigned int** list = neiblist;
    unsigned short* nnum = neibNum;

    #pragma omp parallel for reduction(+:irun) shared(starttime,np,gnp) schedule(dynamic) firstprivate(list,nnum,cvalid)
    for (int i = 1; i <= np; i++)
    {
        #pragma omp atomic
        ++irun;
        testtime = clock();
        clock_t timepassed = testtime - starttime;
        double secondsPassed = timepassed / (double) CLOCKS_PER_SEC;
        
        if (omp_get_thread_num()==0 && secondsPassed>=1) {
            starttime = clock();
            int perc = (int)(((double)irun) / np * 100.0 );
            std::cout << "\rworking at " << irun << "/" << np << " - " << perc << "%" << std::flush;
        }

        if (!cvalid[i-1]) continue;

        int minxcell = 1;
        int maxxcell = 1;
        int minycell = 1;
        int maxycell = 1;
        int minzcell = 1;
        int maxzcell = 1;

        const int dnxgcell = xgcell[i-1] - 1;
        const int dnygcell = ygcell[i-1] - 1;
        const int dnzgcell = zgcell[i-1] - 1;

        const int dpxgcell = xgcell[i-1] + 1;
        const int dpygcell = ygcell[i-1] + 1;
        const int dpzgcell = zgcell[i-1] + 1;

        minxcell = std::max(dnxgcell,1);
        minycell = std::max(dnygcell,1);
        minzcell = std::max(dnzgcell,1);
        maxxcell = std::min(dpxgcell,ngridx);
        maxycell = std::min(dpygcell,ngridy);
        maxzcell = std::min(dpzgcell,ngridz);

        for (auto zcell = minzcell; zcell <= maxzcell; zcell++)
        {
            for (auto ycell = minycell; ycell <= maxycell; ycell++)
            {
                for (auto xcell = minxcell; xcell <= maxxcell; xcell++)
                {
                    int j = grid[xcell-1][ycell-1][zcell-1];

                    tag: if (j>i) 
                    {
                        const double xi = i<=np?pos[i-1].x:gpos[i-np-1].x;
                        const double xj = j<=np?pos[j-1].x:gpos[j-np-1].x;
                        const double yi = i<=np?pos[i-1].y:gpos[i-np-1].y;
                        const double yj = j<=np?pos[j-1].y:gpos[j-np-1].y;
                        const double zi = i<=np?pos[i-1].z:gpos[i-np-1].z;
                        const double zj = j<=np?pos[j-1].z:gpos[j-np-1].z;
                        const double dx = xi - xj;
                        const double dy = yi - yj;
                        const double dz = zi - zj;
                        const double r = sqrt(dx*dx+dy*dy+dz*dz);
                        const double horizon = 3.0*dp;
                        const double horizon2 = 1.5*dp;

                        if (r<horizon) {
                            if (i<=np && j>np) {
                                nnum[j - np - 1]++;
                                if (nnum[j-np-1]>maxlist) 
                                {
                                    std::cerr << "exceed max list length" << nnum[j-np-1];
                                    system("pause");
                                }
                                list[j-np-1][nnum[j-np-1]-1]=i-1;
                            }
                            if (j<=np && i>np) {
                                nnum[i - np - 1]++;
                                if (nnum[i-np-1]>maxlist)
                                {
                                    std::cerr << "exceed max list length" << nnum[i-np-1];
                                    system("pause");
                                }
                                list[i-np-1][nnum[i-np-1]-1]=j-1;
                            }
                        }
                        j = celldata[j-1];
                        goto tag;
                    }
                    
                }
                
            }
            
        }
        
    }

    std::cout << "\rworking at " << np << "/" << np << " - " << 100 << "%\n" << std::flush;

    for (auto i = 0; i < ngridx; i++)
    {
        for (auto j = 0; j < ngridy; j++)
        {
            delete[] grid[i][j];
        }
        delete[] grid[i];
    }
    delete[] grid;
    delete[] xgcell,ygcell,zgcell,celldata,cvalid;

    unsigned total = 0;
    unsigned short mlen = 0;
    for (unsigned i = 0; i < gnp; i++)
    {
        total += neibNum[i];
        mlen = mlen<neibNum[i]?neibNum[i]:mlen;
    }
    //std::cout << "\naveraged background grid node list length: " << total/gnp << std::endl;
    //std::cout << "max neighbor list length: " << mlen << std::endl;

    // neighbor list build done
    // interpolation from particles to grid nodes

    const double horizon = 3.0*dp;
    //unsigned int** list = neiblist;
    //unsigned short* nnum = neibNum;

    #pragma omp parallel for shared(np,gnp,_mk) schedule(dynamic) firstprivate(gmass,gpres,grhop,gvol,gtype,gvel)
    for (int i = 0; i < gnp; i++)
    {
        if (!neibNum[i]) continue;
        double cmass=0,cpress=0,cvolume=0,ctype=0,cmk=0;
        double crhop=0;
        double vx=0,vy=0,vz=0;
        double vcc=0;
        double sumw=0;
        #pragma omp parallel for reduction(+:cmass,cpress,cvolume,ctype,cmk,crhop,vx,vy,vz,vcc,sumw) schedule(dynamic)
        for (auto j = 0; j < neibNum[i]; j++)
        {
          const unsigned k = neiblist[i][j];
          const double dx = gpos[i].x - pos[k].x;
          const double dy = gpos[i].y - pos[k].y;
          const double dz = gpos[i].z - pos[k].z;
          const double r = sqrt(dx*dx + dy*dy + dz*dz);
          const double q = r/dp;
          if (r>horizon) continue;
          const double factor = 1.0/(M_PI * pow(dp,3.0)) * (1.0 - exp(-9.0)) / (1.0 - 10.0 * exp(-9.0));
          const double dfac = factor * exp(-q*q) * (-2.0/dp/dp);
          const double bweight = factor*( exp(-q*q)-exp(-9.0) );
          //if (i==10000) std::cout << j+1 << " " << idp[k] << " bweight: " << bweight << std::endl;
          const double vol_j = rhop[k]?mass[k]/rhop[k]:0;
          vcc += bweight*vol_j;
          sumw+= bweight;
          cmass+= bweight*mass[k]*vol_j;
          cpress+=bweight*pres[k]*vol_j;
          crhop+=bweight*rhop[k]*vol_j;
          cvolume+=bweight*vol[k]*vol_j;
          ctype+=bweight*type[k]*vol_j;
          cmk+=bweight*mk[k]*vol_j;
          vx+=bweight*vel[k].x*vol_j;
          vy+=bweight*vel[k].y*vol_j;
          vz+=bweight*vel[k].z*vol_j;
          //------------------------
        }
        gmass[i] = vcc?cmass/vcc:cmass;
        gpres[i] = vcc?cpress/vcc:cpress;
        grhop[i] = vcc?crhop/vcc:crhop;
        gvol[i] = vcc?cvolume/vcc:cvolume;
        gtype[i] = static_cast<short>(vcc?ctype/vcc:ctype);
        gmk[i] = static_cast<unsigned short>(vcc?cmk/vcc:cmk);
        gvel[i].x = vx;
        gvel[i].y = vy;
        gvel[i].z = vz;
    }

    //interpolation done

    //-Defines variables to save in VTk or CSV.
    JDataArrays arrays,arrays2;
    arrays.AddArray("Pos"  ,gnp,gpos ,false);
    arrays.AddArray("Idp"  ,gnp,gidp ,false);
    arrays.AddArray("Vel"  ,gnp,gvel ,false);
    arrays.AddArray("Mass", gnp,gmass,false);
    arrays.AddArray("Rhop" ,gnp,grhop,false);
    arrays.AddArray("Type" ,gnp,gtype,false);
    arrays.AddArray("Press",gnp,gpres,false);
    arrays.AddArray("Mk"   ,gnp,gmk  ,false);
    arrays.AddArray("DropletR",gnp,gdropr,false);
    arrays.AddArray("GravityCoeff",gnp,ggrav,false);
    arrays2.AddArray("Connectivity",gconn,connectivity,false);

    if(!cfg->SaveGrid.empty()){
      string fileout=(onefile? cfg->SaveGrid: fun::FileNameSec(cfg->SaveGrid,part));
      if(fun::GetExtension(fileout).empty())fileout=fun::AddExtension(fileout,"vtk");
      printf("SaveGridVTK> %s\n",fun::ShortFileName(fileout,68).c_str());
      JVtkLib::SaveVtkGridData(fileout,arrays,"Pos",arrays2);
    }

    // free grid memory
    delete[] gpos; gpos=NULL;
    delete[] gidp; gidp=NULL;
    delete[] gvel; gvel=NULL;
    delete[] gvol; gvol=NULL;
    delete[] grhop; grhop=NULL;
    delete[] gpres; gpres=NULL;
    delete[] gtype; gtype=NULL;
    delete[] gmass; gmass=NULL;
    delete[] gmk; gmk=NULL;
    delete[] gdropr; gdropr=NULL;
    delete[] ggrav; ggrav=NULL;
    for(uint i=0;i<gnp;i++){delete[] neiblist[i];neiblist[i]=NULL;}
    delete[] neiblist; neiblist=NULL;
    delete[] neibNum; neibNum=NULL;
    delete[] connectivity; connectivity=NULL;

    firstdata=false;
    if(!onefile){
      part++;
      file=dirin+JPartDataBi4::GetFileNamePart(part,0,npiece);
    }
    else break;
  }

  // free memory
  delete[] pos;  pos=NULL;
  delete[] posd; posd=NULL;
  delete[] vel;  vel=NULL;
  delete[] rhop; rhop=NULL;
  delete[] idp; idp=NULL;
  delete[] ridp; ridp=NULL;
  delete[] mk;   mk=NULL;
  delete[] type; type=NULL;
  delete[] pres; pres=NULL;
  delete[] dropr; dropr=NULL;
  delete[] grav; grav=NULL;
  delete xmlctes;   xmlctes=NULL;
  delete xmleparms; xmleparms=NULL;
  delete xmlparts;  xmlparts=NULL;
}

//==============================================================================
/// Processes files.
//==============================================================================
void RunFiles(const JCfgRun *cfg){
  const bool outtype=true,outmk=true,outpres=true;

  //-Reads XML file.
  JSpaceCtes* xmlctes=NULL;
  JSpaceEParms* xmleparms=NULL;
  JSpaceParts* xmlparts=NULL;
  if(outmk){
    if(!fun::FileExists(cfg->FileXml))ExceptionFile("XML file of case was not found.",(cfg->FileXml.empty()? string("???"): cfg->FileXml));
    JXml xml; xml.LoadFile(cfg->FileXml);
    xmlctes=new JSpaceCtes();     xmlctes->LoadXmlRun(&xml,"case.execution.constants");
    xmleparms=new JSpaceEParms(); xmleparms->LoadXml(&xml,"case.execution.parameters");
    xmlparts=new JSpaceParts();   xmlparts->LoadXml(&xml,"case.execution.particles");
  }

  //-Prepares input files.
  const bool onefile=!cfg->FileIn.empty();
  const string casein=cfg->FileIn;
  const string dirin=fun::GetDirWithSlash(cfg->DirIn);
  const int last=cfg->Last;
  int part=(onefile||cfg->First<0? 0: cfg->First);

  unsigned npiece=0;
  unsigned casenp=0,casenbound=0,casenfixed=0,casenmoving=0,casenfloat=0,casenfluid=0;
  double cteb=0,rhop0=0,gamma=0;

  double timestep=0;
  unsigned *idp=NULL;
  tfloat3 *pos=NULL,*vel=NULL;
  tdouble3 *posd=NULL;
  float *rhop=NULL;
  unsigned *ridp=NULL;
  byte *type=NULL,*mk=NULL;
  float *pres=NULL;

  //-Processes input files.
  byte npie=0;
  string file=JPartDataBi4::GetFileData(casein,dirin,part,npie);
  if(file.empty())ExceptionText("Error: Data files not found.");
  bool firstdata=true;
  while((last<0||part<=last) && fun::FileExists(file)){
    if(!npiece){//-Reads initial data.
      JPartDataBi4 pd;
      if(onefile)pd.LoadFileCase("",casein,0,npie);
      else pd.LoadFilePart(dirin,part,0,npie);
      npiece=pd.GetNpiece();
      if(npiece>1)ExceptionText("Error: The number of pieces is higher than 1.");
      casenp=std::max((unsigned)pd.Get_CaseNp(),pd.Get_Npok());
      if(std::max((unsigned)pd.Get_CaseNp(),pd.Get_Npok())!=casenp)ExceptionText("Error: The number of particles is too big.");
      if(!pd.Get_IdpSimple())ExceptionText("Error: Only Idp (32 bits) is valid at the moment.");
      casenfluid=(unsigned)pd.Get_CaseNfluid();
      casenfixed=(unsigned)pd.Get_CaseNfixed();
      casenmoving=(unsigned)pd.Get_CaseNmoving();
      casenfloat=(unsigned)pd.Get_CaseNfloat();
      casenbound=casenp-casenfluid;
      cteb=pd.Get_B();
      rhop0=pd.Get_Rhop0();
      gamma=pd.Get_Gamma();
      pos=new tfloat3[casenp];
      posd=new tdouble3[casenp];
      vel=new tfloat3[casenp];
      rhop=new float[casenp];
      idp=new unsigned[casenp];
      //-Allocates memory to other variables.
      if(outmk)mk=new byte[casenp];
      if(outtype)type=new byte[casenp];
      if(outpres)pres=new float[casenp];
      if(outmk)ridp=new unsigned[casenp];
    }//end if npiece

    //-Reads particle data.
    unsigned np=0;
    for(unsigned cp=0;cp<npiece && cp<1;cp++){
      file=(onefile? JPartDataBi4::GetFileNameCase(casein,cp,npiece): dirin+JPartDataBi4::GetFileNamePart(part,cp,npiece));
      if(!cp)printf("load> %s\n",file.c_str());
      JPartDataBi4 pd;
      if(onefile)pd.LoadFileCase("",casein,0,npiece);
      else pd.LoadFilePart(dirin,part,0,npiece);
      if(!cp)timestep=pd.Get_TimeStep();
      const bool possimple=pd.Get_PosSimple();
      const unsigned npok=pd.Get_Npok();
      if(npok){
        //-Loads data from PART.
        pd.Get_Idp(npok,idp);
        pd.Get_Vel(npok,vel);
        pd.Get_Rhop(npok,rhop);
        if(possimple){
          pd.Get_Pos(npok,pos);
        }
        else{ 
          pd.Get_Posd(npok,posd);
          for(unsigned p=0;p<npok;p++)pos[p]=ToTFloat3(posd[p]);
        }
      }
      np=npok;
    }

    //-Loads other vars.
    if(pres)for(unsigned c=0;c<np;c++)pres[c]=(float)(cteb*(pow(rhop[c]/rhop0,gamma)-1.));
    if(type){
      for(unsigned p=0;p<np;p++){
        const unsigned id=idp[p];
        type[p]=(id<casenfixed? 0: (id<casenmoving? 1: (id<casenfloat? 2: 3)));
      }
    }
    if(ridp)for(unsigned c=0;c<np;c++)ridp[idp[c]]=c;
    if(mk){
      for(unsigned c=0;c<xmlparts->CountBlocks();c++){
        const JSpacePartBlock &block=xmlparts->GetBlock(c);
        const byte mkblock=(byte)block.GetMk();
        const unsigned ipend=block.GetBegin()+block.GetCount();
        for(unsigned ip=block.GetBegin();ip<ipend;ip++)mk[ridp[ip]]=mkblock;
      }
    }//end if mk

    //-Defines variables to save in VTk or CSV.
    JDataArrays arrays;
    arrays.AddArray("Pos"  ,np,pos ,false);
    arrays.AddArray("Idp"  ,np,idp ,false);
    arrays.AddArray("Vel"  ,np,vel ,false);
    arrays.AddArray("Rhop" ,np,rhop,false);
    arrays.AddArray("Type" ,np,type,false);
    arrays.AddArray("Press",np,pres,false);
    arrays.AddArray("Mk"   ,np,mk  ,false);

    //-Saves VTK files.
    if(!cfg->SaveVtk.empty()){
      string fileout=(onefile? cfg->SaveVtk: fun::FileNameSec(cfg->SaveVtk,part));
      if(fun::GetExtension(fileout).empty())fileout=fun::AddExtension(fileout,"vtk");
      printf("SaveVTK> %s\n",fun::ShortFileName(fileout,68).c_str());
      JVtkLib::SaveVtkData(fileout,arrays,"Pos");
    }
    //-Saves CSV files.
    if(!cfg->SaveCsv.empty()){
      string fileout=(onefile? cfg->SaveCsv: fun::FileNameSec(cfg->SaveCsv,part));
      if(fun::GetExtension(fileout).empty())fileout=fun::AddExtension(fileout,"csv");
      printf("SaveCSV> %s\n",fun::ShortFileName(fileout,68).c_str());
      JOutputCsv ocsv(false,true);
      ocsv.SaveCsv(fileout,arrays);
    }

    firstdata=false;
    if(!onefile){
      part++;
      file=dirin+JPartDataBi4::GetFileNamePart(part,0,npiece);
    }
    else break;
  }

  //-Free memory.
  delete[] pos;  pos=NULL;
  delete[] posd; posd=NULL;
  delete[] vel;  vel=NULL;
  delete[] rhop; rhop=NULL;
  delete[] idp; idp=NULL;
  delete[] ridp; ridp=NULL;
  delete[] mk;   mk=NULL;
  delete[] type; type=NULL;
  delete[] pres; pres=NULL;
  delete xmlctes;   xmlctes=NULL;
  delete xmleparms; xmleparms=NULL;
  delete xmlparts;  xmlparts=NULL;
}

//==============================================================================
/// GPL License.
//==============================================================================
std::string getlicense_lgpl(const std::string &name){
  std::string tx="";
  tx=tx+"\n\n <"+fun::StrUpper(name)+">  Copyright (c) 2020 by Dr Jose M. Dominguez";
  tx=tx+"\n (see http://dual.sphysics.org/index.php/developers/)\n";
  tx=tx+"\n EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo";
  tx=tx+"\n School of Mechanical, Aerospace and Civil Engineering, University of Manchester\n";
  tx=tx+"\n DualSPHysics is free software: you can redistribute it and/or"; 
  tx=tx+"\n modify it under the terms of the GNU Lesser General Public License";
  tx=tx+"\n as published by the Free Software Foundation, either version 2.1 of"; 
  tx=tx+"\n the License, or (at your option) any later version.\n"; 
  tx=tx+"\n DualSPHysics is distributed in the hope that it will be useful,"; 
  tx=tx+"\n but WITHOUT ANY WARRANTY; without even the implied warranty of"; 
  tx=tx+"\n MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"; 
  tx=tx+"\n GNU Lesser General Public License for more details.\n";
  tx=tx+"\n You should have received a copy of the GNU Lesser General Public License"; 
  tx=tx+"\n along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.\n\n";
  return(tx);
}

//==============================================================================
//==============================================================================
int main(int argc, char** argv){
  int errcode=1;
  printf("%s",getlicense_lgpl("TOVTK4").c_str());
  printf("\n%s\n",APP_NAME);
  for(unsigned c=0;c<=strlen(APP_NAME);c++)printf("="); printf("\n");
  JCfgRun cfg;
  try{
    cfg.LoadArgv(argc,argv);
    if(!cfg.PrintInfo){
      cfg.ValidaCfg();
      RunFilesGrid(&cfg);
    }
    errcode=0;
  }
  catch(const char *cad){
    printf("\n*** Exception: %s\n",cad);
  }
  catch(const string &e){
    printf("\n*** Exception: %s\n",e.c_str());
  }
  catch (const exception &e){
    printf("\n*** %s\n",e.what());
  }
  catch(...){
    printf("\n*** Attention: Unknown exception...\n");
  }
  return(errcode);
}


