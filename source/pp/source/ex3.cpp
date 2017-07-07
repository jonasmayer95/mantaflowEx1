




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/paul/Desktop/mantaflowEx1/source/ex3.cpp"

#include "kernel.h"
#include "grid.h"

using namespace std;
namespace Manta {

	
 struct ComputeDiv : public KernelBase { ComputeDiv(FlagGrid& flags, Grid<Real>& div, MACGrid& vel) :  KernelBase(&flags,1) ,flags(flags),div(div),vel(vel)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& div, MACGrid& vel )  {
	//compute divergence
		if (flags.isFluid(i, j,k)) 
			div(i,j,k) = vel(i+1,j,k).x - vel(i,j,k).x + vel(i,j+1,k).y - vel(i,j,k).y;
	}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return div; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return vel; } typedef MACGrid type2; void runMessage() { debMsg("Executing kernel ComputeDiv ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,div,vel);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,div,vel);  } }  } FlagGrid& flags; Grid<Real>& div; MACGrid& vel;   };
#line 9 "ex3.cpp"



	
 struct iterate_p : public KernelBase { iterate_p(const Grid<Real>& div, Grid<Real>& pressure, const Grid<Real> &A0, const Grid<Real> &A1) :  KernelBase(&div,1) ,div(div),pressure(pressure),A0(A0),A1(A1)   { runMessage(); run(); }  inline void op(int i, int j, int k, const Grid<Real>& div, Grid<Real>& pressure, const Grid<Real> &A0, const Grid<Real> &A1 )  {
		pressure(i,j,k) = 1.0/A1(i,j,k)
				* (div(i,j,k) - A1(i,j,k)*pressure(i,j,k) - A0(i-1,j,k)*pressure(i-1,j,k) - A0(i,j-1,k)*pressure(i,j-1,k) - A0(i+1,j,k)*pressure(i+1,j,k) - A0(i,j+1,k)*pressure(i,j+1,k))
				+ pressure(i,j,k);
	}   inline const Grid<Real>& getArg0() { return div; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return pressure; } typedef Grid<Real> type1;inline const Grid<Real> & getArg2() { return A0; } typedef Grid<Real>  type2;inline const Grid<Real> & getArg3() { return A1; } typedef Grid<Real>  type3; void runMessage() { debMsg("Executing kernel iterate_p ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,pressure,A0,A1);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,pressure,A0,A1);  } }  } const Grid<Real>& div; Grid<Real>& pressure; const Grid<Real> & A0; const Grid<Real> & A1;   };
#line 16 "ex3.cpp"



	
 struct iterate_r : public KernelBase { iterate_r(const Grid<Real>& div, Grid<Real>& pressure, Grid<Real>& residual, const Grid<Real> &A0, const Grid<Real> &A1) :  KernelBase(&div,1) ,div(div),pressure(pressure),residual(residual),A0(A0),A1(A1) ,res(0.0)  { runMessage(); run(); }  inline void op(int i, int j, int k, const Grid<Real>& div, Grid<Real>& pressure, Grid<Real>& residual, const Grid<Real> &A0, const Grid<Real> &A1 ,Real& res)  {
	//r = d - A*pnew
		res = div(i,j,k) - A1(i,j,k)*pressure(i,j,k) - A0(i-1,j,k)*pressure(i-1,j,k) - A0(i,j-1,k)*pressure(i,j-1,k) - A0(i+1,j,k)*pressure(i+1,j,k) - A0(i,j+1,k)*pressure(i,j+1,k);
		residual(i,j,k) = res;
	}   inline operator Real () { return res; } inline Real  & getRet() { return res; }  inline const Grid<Real>& getArg0() { return div; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return pressure; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return residual; } typedef Grid<Real> type2;inline const Grid<Real> & getArg3() { return A0; } typedef Grid<Real>  type3;inline const Grid<Real> & getArg4() { return A1; } typedef Grid<Real>  type4; void runMessage() { debMsg("Executing kernel iterate_r ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  Real res = 0.0; 
#pragma omp for nowait  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,pressure,residual,A0,A1,res); 
#pragma omp critical
{this->res += res; } } } else { const int k=0; 
#pragma omp parallel 
 {  Real res = 0.0; 
#pragma omp for nowait  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,pressure,residual,A0,A1,res); 
#pragma omp critical
{this->res += res; } } }  } const Grid<Real>& div; Grid<Real>& pressure; Grid<Real>& residual; const Grid<Real> & A0; const Grid<Real> & A1;  Real res;  };
#line 23 "ex3.cpp"



	
 struct iterate2 : public KernelBase { iterate2(const Grid<Real>& div, Grid<Real>& pressure, Grid<Real>& residual, const Grid<Real> &A0, const Grid<Real> &A1) :  KernelBase(&div,1) ,div(div),pressure(pressure),residual(residual),A0(A0),A1(A1) ,res(0.0)  { runMessage(); run(); }  inline void op(int i, int j, int k, const Grid<Real>& div, Grid<Real>& pressure, Grid<Real>& residual, const Grid<Real> &A0, const Grid<Real> &A1 ,Real& res)  {
		res = div(i,j,k) - A1(i,j,k)*pressure(i,j,k) - A0(i-1,j,k)*pressure(i-1,j,k) - A0(i,j-1,k)*pressure(i,j-1,k) - A0(i+1,j,k)*pressure(i+1,j,k) - A0(i,j+1,k)*pressure(i,j+1,k);
		residual(i,j,k) = res;
		pressure(i,j,k) = (residual(i,j,k))/A1(i,j,k) + pressure(i,j,k);
	}   inline operator Real () { return res; } inline Real  & getRet() { return res; }  inline const Grid<Real>& getArg0() { return div; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return pressure; } typedef Grid<Real> type1;inline Grid<Real>& getArg2() { return residual; } typedef Grid<Real> type2;inline const Grid<Real> & getArg3() { return A0; } typedef Grid<Real>  type3;inline const Grid<Real> & getArg4() { return A1; } typedef Grid<Real>  type4; void runMessage() { debMsg("Executing kernel iterate2 ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  Real res = 0.0; 
#pragma omp for nowait  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,pressure,residual,A0,A1,res); 
#pragma omp critical
{this->res += res; } } } else { const int k=0; 
#pragma omp parallel 
 {  Real res = 0.0; 
#pragma omp for nowait  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,pressure,residual,A0,A1,res); 
#pragma omp critical
{this->res += res; } } }  } const Grid<Real>& div; Grid<Real>& pressure; Grid<Real>& residual; const Grid<Real> & A0; const Grid<Real> & A1;  Real res;  };
#line 30 "ex3.cpp"



	
 struct UpdateVel : public KernelBase { UpdateVel(const FlagGrid& flags, MACGrid& vel, const Grid<Real>& pressure) :  KernelBase(&flags,1) ,flags(flags),vel(vel),pressure(pressure)   { runMessage(); run(); }  inline void op(int i, int j, int k, const FlagGrid& flags, MACGrid& vel, const Grid<Real>& pressure )  {
		if (flags.isFluid(i, j,k))
			vel(i,j,k) = vel(i,j,k) - Vec3(pressure(i-1,j,k)-pressure(i,j,k), pressure(i,j-1,k)-pressure(i,j,k), vel(i,j,k).z); 
	}   inline const FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline MACGrid& getArg1() { return vel; } typedef MACGrid type1;inline const Grid<Real>& getArg2() { return pressure; } typedef Grid<Real> type2; void runMessage() { debMsg("Executing kernel UpdateVel ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,vel,pressure);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,vel,pressure);  } }  } const FlagGrid& flags; MACGrid& vel; const Grid<Real>& pressure;   };
#line 37 "ex3.cpp"



	
 struct totalNumber : public KernelBase { totalNumber(Grid<Real>& div) :  KernelBase(&div,1) ,div(div) ,total(0.0)  { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real>& div ,Real& total)  {
		total +=1;
	}   inline operator Real () { return total; } inline Real  & getRet() { return total; }  inline Grid<Real>& getArg0() { return div; } typedef Grid<Real> type0; void runMessage() { debMsg("Executing kernel totalNumber ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  Real total = 0.0; 
#pragma omp for nowait  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,total); 
#pragma omp critical
{this->total += total; } } } else { const int k=0; 
#pragma omp parallel 
 {  Real total = 0.0; 
#pragma omp for nowait  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,div,total); 
#pragma omp critical
{this->total += total; } } }  } Grid<Real>& div;  Real total;  };
#line 43 "ex3.cpp"



	
 struct totalDivr : public KernelBase { totalDivr(FlagGrid& flags, Grid<Real>& div, MACGrid& vel) :  KernelBase(&flags,1) ,flags(flags),div(div),vel(vel) ,totalDiv(0.0)  { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& div, MACGrid& vel ,Real& totalDiv)  {
	//compute divergence
		if (flags.isFluid(i, j,k)) 
			totalDiv += vel(i+1,j,k).x - vel(i,j,k).x + vel(i,j+1,k).y - vel(i,j,k).y;
	}   inline operator Real () { return totalDiv; } inline Real  & getRet() { return totalDiv; }  inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return div; } typedef Grid<Real> type1;inline MACGrid& getArg2() { return vel; } typedef MACGrid type2; void runMessage() { debMsg("Executing kernel totalDivr ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  Real totalDiv = 0.0; 
#pragma omp for nowait  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,div,vel,totalDiv); 
#pragma omp critical
{this->totalDiv += totalDiv; } } } else { const int k=0; 
#pragma omp parallel 
 {  Real totalDiv = 0.0; 
#pragma omp for nowait  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,div,vel,totalDiv); 
#pragma omp critical
{this->totalDiv += totalDiv; } } }  } FlagGrid& flags; Grid<Real>& div; MACGrid& vel;  Real totalDiv;  };
#line 48 "ex3.cpp"



	
 struct ComputeA0 : public KernelBase { ComputeA0(const FlagGrid& flags, Grid<Real> &A0) :  KernelBase(&flags,0) ,flags(flags),A0(A0)   { runMessage(); run(); }  inline void op(int i, int j, int k, const FlagGrid& flags, Grid<Real> &A0 )  {
		if(flags.isFluid(i, j,k))
			A0(i,j,k) = -1;
		else 
			A0(i,j,k) = 0;
	}   inline const FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real> & getArg1() { return A0; } typedef Grid<Real>  type1; void runMessage() { debMsg("Executing kernel ComputeA0 ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,A0);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,A0);  } }  } const FlagGrid& flags; Grid<Real> & A0;   };
#line 55 "ex3.cpp"



	
 struct ComputeA1 : public KernelBase { ComputeA1(const Grid<Real> &A0, Grid<Real> &A1) :  KernelBase(&A0,1) ,A0(A0),A1(A1)   { runMessage(); run(); }  inline void op(int i, int j, int k, const Grid<Real> &A0, Grid<Real> &A1 )  {
		A1(i,j,k) = 0.0-(A0(i-1,j,k) + A0(i,j-1,k) + A0(i+1,j,k) + A0(i,j+1,k));
	}   inline const Grid<Real> & getArg0() { return A0; } typedef Grid<Real>  type0;inline Grid<Real> & getArg1() { return A1; } typedef Grid<Real>  type1; void runMessage() { debMsg("Executing kernel ComputeA1 ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,A0,A1);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,A0,A1);  } }  } const Grid<Real> & A0; Grid<Real> & A1;   };
#line 63 "ex3.cpp"



	void solvePressureGS(FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure, Real gsAccuracy = 1e-4) {
		FluidSolver* parent = flags.getParent();
		Grid<Real> div(parent);
		Real res = 1.0;
		Grid<Real> residual(parent);
		Grid<Real> A0(parent);
		Grid<Real> A1(parent);

		ComputeDiv(flags, div, vel);
		ComputeA0(flags, A0);
		ComputeA1(A0, A1);

		while(res > gsAccuracy)
		{
			iterate_p(div, pressure, A0, A1);
			res = iterate_r(div, pressure, residual, A0, A1);
			//res = iterate2(div, pressure, residual, A0, A1);
		}
		UpdateVel(flags, vel, pressure);

	} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "solvePressureGS" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); MACGrid& vel = *_args.getPtr<MACGrid >("vel",1,&_lock); Grid<Real>& pressure = *_args.getPtr<Grid<Real> >("pressure",2,&_lock); Real gsAccuracy = _args.getOpt<Real >("gsAccuracy",3,1e-4,&_lock);   _retval = getPyNone(); solvePressureGS(flags,vel,pressure,gsAccuracy);  _args.check(); } pbFinalizePlugin(parent,"solvePressureGS", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("solvePressureGS",e.what()); return 0; } } static const Pb::Register _RP_solvePressureGS ("","solvePressureGS",_W_0);  extern "C" { void PbRegister_solvePressureGS() { KEEP_UNUSED(_RP_solvePressureGS); } } 

	Real getMaxDivergence(MACGrid& vel, FlagGrid& flags){
		FluidSolver* parent = flags.getParent();
		Grid<Real> div(parent);
		return totalDivr(flags, div, vel)/totalNumber(div);
	} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "getMaxDivergence" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; MACGrid& vel = *_args.getPtr<MACGrid >("vel",0,&_lock); FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",1,&_lock);   _retval = toPy(getMaxDivergence(vel,flags));  _args.check(); } pbFinalizePlugin(parent,"getMaxDivergence", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("getMaxDivergence",e.what()); return 0; } } static const Pb::Register _RP_getMaxDivergence ("","getMaxDivergence",_W_1);  extern "C" { void PbRegister_getMaxDivergence() { KEEP_UNUSED(_RP_getMaxDivergence); } } 

} // end namespace


