

//  Dendrite Growth.cpp
//
//  Created by zhu on 6/20/17.
//  Copyright © 2017 zhu. All rights reserved.


#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define NDP 351
#define PI 3.14159
#define RR 8.3145

int ndx=NDP-1,ndy=NDP-1;
int ndmx=NDP-2, ndmy=NDP-2;
double s1h[NDP][NDP], Th[NDP][NDP];


double Tm, Tini;
int numsave;

void ini000();
void datasave();

int main(void)
{
    
//************ VARIABLES ************//
    
    double s1h2[NDP][NDP], Th2[NDP][NDP];//補助配列
    //*** Length Scale ***//
    double dx, dy;//差分格子サイズ
    double al;//計算領域の１辺の長さ
    
    //*** Time Scale ***//
    double delt;
    double dtp;//時間刻み
    double dtt;	//時間刻み
    double time1;
    double time1max;//計算カウント数の最大値
    
    //*** Phase Field ***//
    double s1;//フェーズフィールド
    int i,j,ip,im,jp,jm;// 格子位置
    double s1ip,s1im,s1jp,s1jm;//周囲のフェーズフィールド
    double s1ipjp, s1ipjm, s1imjp, s1imjm;
    double ep, ep1p, ep2p;//異方性項微分
    double dx_s1, dy_s1;//フェーズフィールド微分
    double dxx_s1, dyy_s1, dxy_s1;
    double s1kai, s1kais;//ポテンシャル
    double s1ddtt;
    
    //*** Temperature Field ***//
    double TT;
    double Tip, Tim, Tjp, Tjm;
    double Tddtt;
    
    //*** Material Parameters ***//
    
    double cndct;//熱伝導率
    double speht;//比熱
    double rlate;//潜熱
    
    //*** Interface Parameters ***//
    double astre;//異方性強度
    double gamma;//界面エネルギー密度
    double delta;//界面幅
    double skine;//界面カイネティック係数
    double aaa;//勾配エネルギー係数
    double www;	//ペナルティー項のエネルギー障壁
    double ram;	//λ
    double bbb;	//界面幅に関係する係数
    double anois;//ノイズの振幅
    double dF;//駆動力
    double pmobi;//フェーズフィールドのモビリティ
    
    //*** Geometry Parameters ***//
    double th0;	//優先成長方向の角度
    double th;//界面の法線方向の角度
    double j_fold;//異方性モード
    
    //*** Flow Control ***//
    double dami1, dami2;//場の計算スキップ変数


    //*** Value Assignment ***//
    dx=dy=30.0e-9;
    delta=3.0*dx;
    //dx=dy=20.0e-9;
    //delta=4.0*dx;
    al = dx*(double)ndx;
    gamma = 0.37;
    ram = 0.1;
    bbb=2.0*log((1.0+(1.0-2.0*ram))/(1.0-(1.0-2.0*ram)))/2.0;
    j_fold=4.0;
    anois=0.1;
    
    //astre=0.005;
    astre=0.01;
    //astre=0.03;
    //astre=0.05;
    
    //pure Ni
    cndct=84.01;
    speht=5.42e+06;
    rlate=2.350e+09;
    Tm=1728.0;
    Tini=1511.2;
    skine=2.0;
    
    th0=0.7;
    
    aaa=sqrt(3.0*delta*gamma/bbb);
    www=6.0*gamma*bbb/delta;
    pmobi=bbb*Tm*skine/(3.0*delta*rlate);
    
    dtp=dx*dx/(5.0*pmobi*aaa*aaa);
    dtt=dx*dx/(5.0*cndct/speht);
    if(dtp>dtt){delt=dtt;} else{delt=dtp;}
    
    time1max = 10e4;
    
    numsave=0;
    ini000();
    

    
//************ CALULATIONS ************//
start: ;
    
    if ((int)time1 % 100==0){numsave++; datasave();}
    
    
    for(i=0;i<=ndx;i++){
        for(j=0;j<=ndy;j++){
            ip=i+1; im=i-1; jp=j+1; jm=j-1;
            
            //*** 境界条件 ***//
            if(i==ndx){ip=ndmx;}  if(i==0){im=1;}
            if(j==ndy){jp=ndmy;}  if(j==0){jm=1;}
            
            s1=s1h[i][j];
            s1ip=s1h[ip][j];
            s1im=s1h[im][j];
            s1jp=s1h[i][jp];
            s1jm=s1h[i][jm];
            s1ipjp=s1h[ip][jp];
            s1ipjm=s1h[ip][jm];
            s1imjp=s1h[im][jp];
            s1imjm=s1h[im][jm];
            
            TT=Th[i][j];
            Tip=Th[ip][j];
            Tim=Th[im][j];
            Tjp=Th[i][jp];
            Tjm=Th[i][jm];
            
            //*** 計算をスキャップ判断 ***//
            dami1=fabs(s1+s1ip+s1im+s1jp+s1im); dami2=fabs(TT+Tip+Tim+Tjp+Tjm-5.0*Tini);
            if ((dami1<=1.0e-20)&&(dami2<=1.0e-20))
            {
                s1h2[i][j]=s1h[i][j]; Th2[i][j]=Th[i][j];
                goto dami;
            }
            
            //*** フェーズフィールドの微分 ***//
            dx_s1=(s1ip-s1im)/2.0/dx;
            dy_s1=(s1jp-s1jm)/2.0/dy;
            dxx_s1=(s1ip+s1im-2.0*s1)/dx/dx;
            dyy_s1=(s1jp+s1jm-2.0*s1)/dy/dy;
            dxy_s1=(s1ipjp+s1imjm-s1imjp-s1ipjm)/4.0/dx/dy;
            
            //*** Derivative of Anisotropic Strength ***//
            th=atan(dy_s1/(dx_s1+1.0e-20));
            ep=aaa*(1.0+astre*cos(j_fold*(th-th0)));
            ep1p=-aaa*astre*j_fold*sin(j_fold*(th-th0));
            ep2p=-aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));
            
            
            //*** ポテンシャルの計算　***//
            s1kais=-ep*ep*(dxx_s1+dyy_s1)
            -ep*ep1p*((dyy_s1-dxx_s1)*sin(2.0*th)+2.0*dxy_s1*cos(2.0*th))
            +0.5*(ep1p*ep1p+ep*ep2p)*(2.0*dxy_s1*sin(2.0*th)
                                      -dxx_s1-dyy_s1-(dyy_s1-dxx_s1)*cos(2.0*th));
            
            dF=15.0/(2.0*www)*rlate*(TT-Tm)/Tm*s1*(1.0-s1);
            s1kai=4.0*www*s1*(1.0-s1)*(0.5-s1+dF+anois*(DRND(1)-0.5));
            
            //*** 発展方程式 ***//
            s1ddtt=-pmobi*(s1kai+s1kais);
            Tddtt=( cndct*( (Tip+Tim-2.0*TT)/dx/dx+(Tjp+Tjm-2.0*TT)/dy/dy )
                   +30.0*s1*(1.0-s1)*s1*(1.0-s1)*rlate*s1ddtt )/speht;
            
            
            s1h2[i][j]=s1+s1ddtt*delt;
            //s1h2[i][j]=s1+s1ddtt*delt+anois*(DRND(1)-0.5)*s1*(1.0-s1);
            if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}
            if(s1h2[i][j]<=0.0){s1h2[i][j]=0.0;}
            Th2[i][j]=Th[i][j]+Tddtt*delt;
            
        dami:;
        }
    }
    
    //場の補助配列を主配列に移動
    for(i=0;i<=ndx;i++){
        for(j=0;j<=ndy;j++){
            s1h[i][j]=s1h2[i][j];
            Th[i][j]=Th2[i][j];    }}
    
    time1 = time1+1;
    if (time1<time1max){ goto start;}
    
end:;
    
    return 0;
}

//************ FUNCTIONS ************//

//*** Initialization of Field ***//
void ini000()
{
    int i,j;
    srand(time(NULL));// 乱数初期化
    for(i=0;i<=ndx;i++){
        for(j=0;j<=ndy;j++){
            s1h[i][j]=0.0;
            if(((i-120)*(i-120)+(j-120)*(j-120))<4.){s1h[i][j]=0.9;}
            //if(((i-200)*(i-200)+(j)*(j))<100.){s1h[i][j]=0.9;}
            //if((i<10.)&&(j<10.)){s1h[i][j]=0.9;}
            
        }
    }
    for(i=0;i<=ndx;i++){
        for(j=0;j<=ndy;j++){
            Th[i][j]=Tini+s1h[i][j]*(Tm-Tini);
        }
    }
    
}


//*** Data Saving ***//
void datasave()
{
    FILE *s1stream, *tstream; // stream pointer
    char number[20];
    char filename_s1[256]="Phase_field/phi";
    char filename_t[256]="Temperature_field/t";
    
    
    sprintf(number,"%d",numsave);
    strcat(filename_s1, number);
    strcat(filename_t, number);
    strcat(filename_s1, ".dat");
    strcat(filename_t, ".dat");
    
    s1stream = fopen(filename_s1, "w+"); // open
    tstream = fopen(filename_t, "w+");
    for(int i=0; i<=ndx; i++){
        for(int j=0;j<=ndy;j++){
            fprintf(s1stream, "%e ",s1h[i][j]);
            fprintf(tstream, "%e ",Th[i][j]);
        }
        fprintf(s1stream, "\n");
        fprintf(tstream, "\n");
    }
    
    fprintf(s1stream, "\n");
    fclose(s1stream);
    fprintf(tstream, "\n");
    fclose(tstream);
    
}































