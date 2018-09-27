/***************************************************************************** 
 * Project: RooBernsteinEffi                                                 * 
 *                                                                           * 
 * P.Dini fecit, Anno Domini MMXVIII                                         *
 * "Est modus in rebus."                                                     * 
 *                                                                           * 
 * Class to describe 3D angular efficiency in B0->K*MuMU Analysis            * 
 *                                                                           * 
 *****************************************************************************/ 


#include "Riostream.h" 

#include "RooBernsteinEffi.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooListProxy.h"
#include "RooAbsReal.h"
#include <math.h> 
#include "TMath.h"  

ClassImp(RooBernsteinEffi); 
   
  
   RooBernsteinEffi::RooBernsteinEffi(const char *name, const char *title, 
                        RooAbsReal& x, RooAbsReal& y, RooAbsReal& z, 
			const RooArgList& coefList,
			int maxDegree1, int maxDegree2, int maxDegree3
			) :
   RooAbsReal(name,title), 
   _x("_x","_x",this,x),
   _y("_y","_y",this,y),
   _z("_z","_z",this,z), 
   _coefList("_coefList","_coefList",this),
   _maxDegree1(maxDegree1),
   _maxDegree2(maxDegree2),
   _maxDegree3(maxDegree3)
{ 
    TIterator* coefIter = coefList.createIterator() ;
     RooAbsArg* coef ;
     while((coef = (RooAbsArg*)coefIter->Next())) {
       if (!dynamic_cast<RooAbsReal*>(coef)) {
    	 std::cout << "RooBernsteinEffi::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
    	      << " is not of type RooAbsReal" << std::endl ;
    	 assert(0) ;
       }
       _coefList.add(*coef) ;
//   	 std::cout << "RooBernsteinEffi::ctor(" << GetName() << ") coefficient " << coef->GetName()<<std::endl ;
      }
     delete coefIter;
 } 

 RooBernsteinEffi::RooBernsteinEffi(const RooBernsteinEffi& other, const char* name) :  
   RooAbsReal(other,name), 
   _x("_x",this,other._x),
   _y("_y",this,other._y),
   _z("_z",this,other._z), 
   _coefList("_coefList",this,other._coefList),
   _maxDegree1(other._maxDegree1),
   _maxDegree2(other._maxDegree2),
   _maxDegree3(other._maxDegree3)
 { 
 } 
 Double_t RooBernsteinEffi::evaluate() const 
 { 
   fptype x    = _x; // x, y, z...
   fptype y    = _y; // x, y, z...
   fptype z    = _z; // x, y, z...
   int maxDegree1  = _maxDegree1;
   int maxDegree2  = _maxDegree2;
   int maxDegree3  = _maxDegree3;
   
    fptype xmin = _x.min();
    fptype xdif = _x.max()-xmin;
    x=(x-xmin)/xdif;
    fptype ymin = _y.min();
    fptype ydif = _y.max()-ymin;
    y=(y-ymin)/ydif;
    fptype zmin = _z.min();
    fptype zdif =_z.max()-zmin;
    z=(z-zmin)/zdif;
 
       int ipar =0;
       fptype func =0.0;
       for(int i = 0; i <= maxDegree1 ; ++i) {
         for(int j = 0; j <= maxDegree2 ; ++j) {
          for(int k = 0; k <= maxDegree3 ; ++k) {
           fptype bernknvalx =  device_bernsteinkn_func(x,maxDegree1,i);
	   fptype bernknvaly =  device_bernsteinkn_func(y,maxDegree2,j);
	   fptype bernknvalz =  device_bernsteinkn_func(z,maxDegree3,k);
	   func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) _coefList[ipar]).getVal()<<std::endl;
	   ipar++;
	  }
         }
       }
//        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//       exit(0);
      if(func<1.E-30)  return 1.E-30;
      return  func/(xdif*ydif*zdif);
  } 

fptype RooBernsteinEffi::device_coeffbinomial  (fptype enne, fptype kappa) const {
        fptype factor=1.;
        for(fptype i = 1; i <=kappa; ++i) {
          factor *= (enne+1-i)/i; 
        }	 
        if (factor<=0 ){
	 printf("Error in RooBernsteinEffi coeffbinomial=> factor = %f enne=%f kappa=%f",factor,enne,kappa);
         return 0;
	} 
       return factor;
}

fptype  RooBernsteinEffi::device_bernsteinkn_func  (fptype x, fptype enne, fptype kappa) const{
   return RooBernsteinEffi::device_coeffbinomial(enne,kappa)*pow(x,kappa)*pow(1.0-x,enne-kappa);
}

 fptype RooBernsteinEffi::evaluateInt(fptype xBinw,fptype yBinw,fptype zBinw) const {
    int maxDegree1      = _maxDegree1;
    int maxDegree2      = _maxDegree2;
    int maxDegree3      = _maxDegree3;
    fptype x    = _x; 
    fptype y    = _y; 
    fptype z    = _z; 

//     fptype xBinw = 2.0/25.;
//     fptype yBinw = 2.0/25.;
//     fptype zBinw = 2.0*TMath::Pi()/25.;


    fptype xmin = _x.min();
    fptype xdif = _x.max() - xmin ;

    fptype xLeft  = ((x-xBinw/2.)-xmin)/xdif;
    fptype xRight = ((x+xBinw/2.)-xmin)/xdif;
//    x=(x-xmin)/xdif;
    fptype ymin = _y.min();
    fptype ydif = _y.max()-ymin;
    fptype yLeft  = ((y-yBinw/2.)-ymin)/ydif;
    fptype yRight = ((y+yBinw/2.)-ymin)/ydif;
//    y=(y-ymin)/ydif;
    fptype zmin = _z.min();
    fptype zdif = _z.max()-zmin;
    fptype zLeft  = ((z-zBinw/2.)-zmin)/zdif;
    fptype zRight = ((z+zBinw/2.)-zmin)/zdif;
//    z=(z-zmin)/zdif;
    
       int ipar = 0 ;
       fptype ret   =0;
       for(int i = 0; i <= maxDegree1 ; ++i) {
         for(int j = 0; j <= maxDegree2 ; ++j) {
//	  std::cout<<"func = par["<<ipar<<"]*x^"<<kk<<"*y^"<<jj<<std::endl;
          for(int k = 0; k <= maxDegree3 ; ++k) {

           fptype bernknintgbinx = device_EffiBernsteinkn_intgBin(xLeft,xRight,maxDegree1,i);
           fptype bernknintgbiny = device_EffiBernsteinkn_intgBin(yLeft,yRight,maxDegree2,j);
           fptype bernknintgbinz = device_EffiBernsteinkn_intgBin(zLeft,zRight,maxDegree3,k);
           ret   +=((RooAbsReal&) _coefList[ipar]).getVal()*bernknintgbinx*bernknintgbiny*bernknintgbinz;
	   
	   ipar++;
	  }
	  
	 
         }
       }
    ret=ret/(xBinw*yBinw*zBinw);
    if(ret<1.E-30) ret = 1.E-30;

   return ret;

 }
fptype  RooBernsteinEffi::device_EffiBernsteinkn_intgBin( fptype xLeft, fptype xRight, fptype enne, fptype kappa) const{
 
      fptype integbernkn = 0.0;
      fptype ifactni = 0.0;
      fptype ifactik = 0.0;
      
       for(fptype i = kappa; i <=enne ; ++i) {
// n!/(i!(n-i)!)

        ifactni =  device_coeffbinomial(enne,i);
// i!/(k!(i-k)!)
 
        ifactik =  device_coeffbinomial(i,kappa);
        integbernkn += ifactni*ifactik*pow(-1.0,i-kappa)*(pow(xRight, i+1)-pow(xLeft,i+1))/(i+1);
       }

       if (integbernkn<=0.0 ){
          printf(" Error in EffiBernsteinkn_intgbin xLeft=%f xRight=%f kappa=%f enne=%f integral = %5.15f\n",xLeft,xRight,kappa,enne,integbernkn);
        integbernkn=1.E-30;
       }
       return integbernkn;
}
