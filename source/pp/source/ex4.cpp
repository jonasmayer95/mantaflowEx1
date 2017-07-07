




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/paul/Desktop/mantaflowEx1/source/ex4.cpp"
/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 ******************************************************************************/
#include <thread>
#include "kernel.h"
#include "grid.h"
#include "pcgsolver.h"

using namespace std;
namespace Manta {
	

 struct computeDiffuseTemperatureExplicit : public KernelBase { computeDiffuseTemperatureExplicit(FlagGrid& flags, Grid<Real>& T0, Real alpha) :  KernelBase(&flags,1) ,flags(flags),T0(T0),alpha(alpha)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, Grid<Real>& T0, Real alpha )  {
	Real dt = flags.getParent()->getDt();
	if(flags.isFluid(i,j,k)){
		T0(i,j,k) = T0(i,j,k) + dt * alpha * (T0(i+1,j+1,k) + T0(i-1,j-1,k) - 4 * T0(i,j,k) + T0(i-1,j+1,k) + T0(i+1,j-1,k));
	}
	else{
		T0(i,j,k) = 0;
	} 

}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline Grid<Real>& getArg1() { return T0; } typedef Grid<Real> type1;inline Real& getArg2() { return alpha; } typedef Real type2; void runMessage() { debMsg("Executing kernel computeDiffuseTemperatureExplicit ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,T0,alpha);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=1; j < _maxY; j++) for (int i=1; i < _maxX; i++) op(i,j,k,flags,T0,alpha);  } }  } FlagGrid& flags; Grid<Real>& T0; Real alpha;   };
#line 20 "ex4.cpp"



// 4.2 
void diffuseTemperatureExplicit(FlagGrid& flags, Grid<Real>& temperature, Real alpha){
	// don't overwrite values in T that will be read again
	// write a KERNEL and make sure that the temperature in boundary cells stays zero
	computeDiffuseTemperatureExplicit(flags, temperature, alpha);
} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "diffuseTemperatureExplicit" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& temperature = *_args.getPtr<Grid<Real> >("temperature",1,&_lock); Real alpha = _args.get<Real >("alpha",2,&_lock);   _retval = getPyNone(); diffuseTemperatureExplicit(flags,temperature,alpha);  _args.check(); } pbFinalizePlugin(parent,"diffuseTemperatureExplicit", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("diffuseTemperatureExplicit",e.what()); return 0; } } static const Pb::Register _RP_diffuseTemperatureExplicit ("","diffuseTemperatureExplicit",_W_0);  extern "C" { void PbRegister_diffuseTemperatureExplicit() { KEEP_UNUSED(_RP_diffuseTemperatureExplicit); } } 


// 4.3
 struct setupB : public KernelBase { setupB(Grid<Real>& T0,std::vector<Real>* b) :  KernelBase(&T0,0) ,T0(T0),b(b)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real>& T0,std::vector<Real>* b )  {
	(*b)[i*T0.getStrideX()+j*T0.getStrideY()] = T0(i,j,k);
}   inline Grid<Real>& getArg0() { return T0; } typedef Grid<Real> type0;inline std::vector<Real>* getArg1() { return b; } typedef std::vector<Real> type1; void runMessage() { debMsg("Executing kernel setupB ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,T0,b);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,T0,b);  } }  } Grid<Real>& T0; std::vector<Real>* b;   };
#line 40 "ex4.cpp"



 struct fillT : public KernelBase { fillT(FlagGrid& flags,std::vector<Real>& x,Grid<Real>& T) :  KernelBase(&flags,0) ,flags(flags),x(x),T(T)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags,std::vector<Real>& x,Grid<Real>& T )  {
	// make sure that temperature values in boundary cells are zero!
	if(flags.isFluid(i,j,k))
		T(i,j,k) = x[ i*flags.getStrideX()+j*flags.getStrideY()];
	else
		T(i,j,k) = 0;
		
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline std::vector<Real>& getArg1() { return x; } typedef std::vector<Real> type1;inline Grid<Real>& getArg2() { return T; } typedef Grid<Real> type2; void runMessage() { debMsg("Executing kernel fillT ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,x,T);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,flags,x,T);  } }  } FlagGrid& flags; std::vector<Real>& x; Grid<Real>& T;   };
#line 44 "ex4.cpp"



// use the single argument to prevent multithreading (multiple threads might interfere during the matrix setup)
 struct setupA : public KernelBase { setupA(FlagGrid& flags, SparseMatrix<Real>& A, int N, Real alpha, Real dt) :  KernelBase(&flags,0) ,flags(flags),A(A),N(N),alpha(alpha),dt(dt)   { runMessage(); run(); }  inline void op(int i, int j, int k, FlagGrid& flags, SparseMatrix<Real>& A, int N, Real alpha, Real dt )  { 
	// set with:  A.set_element( index1, index2 , value );
	// if needed, read with: A(index1, index2);

	// avoid zero rows in A -> set the diagonal value for boundary cells to 1.0
	if(flags.isFluid(i,j,k))
	{
		A.set_element( i*flags.getStrideX()+j*flags.getStrideY(), i*flags.getStrideX()+j*flags.getStrideY() , 1.+4*alpha*dt );

		A.set_element( i*flags.getStrideX()+j*flags.getStrideY(), i+1*flags.getStrideX()+j*flags.getStrideY() , -alpha*dt );
		A.set_element( i*flags.getStrideX()+j*flags.getStrideY(), i-1*flags.getStrideX()+j*flags.getStrideY() , -alpha*dt );
		A.set_element( i*flags.getStrideX()+j*flags.getStrideY(), i*flags.getStrideX()+j+1*flags.getStrideY() , -alpha*dt );
		A.set_element( i*flags.getStrideX()+j*flags.getStrideY(), i*flags.getStrideX()+j-1*flags.getStrideY() , -alpha*dt );
	}
	else
		A.set_element( i*flags.getStrideX()+j*flags.getStrideY(), i*flags.getStrideX()+j*flags.getStrideY() , 1. );
}   inline FlagGrid& getArg0() { return flags; } typedef FlagGrid type0;inline SparseMatrix<Real>& getArg1() { return A; } typedef SparseMatrix<Real> type1;inline int& getArg2() { return N; } typedef int type2;inline Real& getArg3() { return alpha; } typedef Real type3;inline Real& getArg4() { return dt; } typedef Real type4; void runMessage() { debMsg("Executing kernel setupA ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; for (int k=minZ; k< maxZ; k++) for (int j=0; j< _maxY; j++) for (int i=0; i< _maxX; i++) op(i,j,k, flags,A,N,alpha,dt);  } FlagGrid& flags; SparseMatrix<Real>& A; int N; Real alpha; Real dt;   };


void diffuseTemperatureImplicit(FlagGrid& flags, Grid<Real>& T0, Real alpha){
	// solve A T = b
	const int N = T0.getSizeX()*T0.getSizeY();
	Real dt = flags.getParent()->getDt();
	SparseMatrix<Real> A(N);
	std::vector<Real> b(N);

	setupA(flags, A, N, alpha, dt);
	setupB(T0, &b);

	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(N);
	for (int j = 0; j<N; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(A, b, x, ret_pcg_residual, ret_pcg_iterations, 0);
	// x contains the new temperature values
	fillT(flags, x, T0);
} static PyObject* _W_1 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "diffuseTemperatureImplicit" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; FlagGrid& flags = *_args.getPtr<FlagGrid >("flags",0,&_lock); Grid<Real>& T0 = *_args.getPtr<Grid<Real> >("T0",1,&_lock); Real alpha = _args.get<Real >("alpha",2,&_lock);   _retval = getPyNone(); diffuseTemperatureImplicit(flags,T0,alpha);  _args.check(); } pbFinalizePlugin(parent,"diffuseTemperatureImplicit", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("diffuseTemperatureImplicit",e.what()); return 0; } } static const Pb::Register _RP_diffuseTemperatureImplicit ("","diffuseTemperatureImplicit",_W_1);  extern "C" { void PbRegister_diffuseTemperatureImplicit() { KEEP_UNUSED(_RP_diffuseTemperatureImplicit); } } 

}; 


