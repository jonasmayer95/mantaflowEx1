




// DO NOT EDIT !
// This file is generated using the MantaFlow preprocessor (prep generate).




#line 1 "/home/paul/Desktop/mantaflowEx1/source/ex1.cpp"
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

#include "kernel.h"
#include "grid.h"

using namespace std;
namespace Manta {

	
 struct levelsetToDensityKN : public KernelBase { levelsetToDensityKN(Grid<Real>& phi, Grid<Real>& density) :  KernelBase(&phi,0) ,phi(phi),density(density)   { runMessage(); run(); }  inline void op(int i, int j, int k, Grid<Real>& phi, Grid<Real>& density )  {
		// map phi values [-5,5] to density values [1,0]
	}   inline Grid<Real>& getArg0() { return phi; } typedef Grid<Real> type0;inline Grid<Real>& getArg1() { return density; } typedef Grid<Real> type1; void runMessage() { debMsg("Executing kernel levelsetToDensityKN ", 3); debMsg("Kernel range" <<  " x "<<  maxX  << " y "<< maxY  << " z "<< minZ<<" - "<< maxZ  << " "   , 4); }; void run() {  const int _maxX = maxX; const int _maxY = maxY; if (maxZ > 1) { 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int k=minZ; k < maxZ; k++) for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,phi,density);  } } else { const int k=0; 
#pragma omp parallel 
 {  
#pragma omp for  
  for (int j=0; j < _maxY; j++) for (int i=0; i < _maxX; i++) op(i,j,k,phi,density);  } }  } Grid<Real>& phi; Grid<Real>& density;   };
#line 19 "ex1.cpp"



	void levelsetToDensity(Grid<Real>& phi, Grid<Real>& density){
		levelsetToDensityKN(phi, density);
	} static PyObject* _W_0 (PyObject* _self, PyObject* _linargs, PyObject* _kwds) { try { PbArgs _args(_linargs, _kwds); FluidSolver *parent = _args.obtainParent(); bool noTiming = _args.getOpt<bool>("notiming", -1, 0); pbPreparePlugin(parent, "levelsetToDensity" , !noTiming ); PyObject *_retval = 0; { ArgLocker _lock; Grid<Real>& phi = *_args.getPtr<Grid<Real> >("phi",0,&_lock); Grid<Real>& density = *_args.getPtr<Grid<Real> >("density",1,&_lock);   _retval = getPyNone(); levelsetToDensity(phi,density);  _args.check(); } pbFinalizePlugin(parent,"levelsetToDensity", !noTiming ); return _retval; } catch(std::exception& e) { pbSetError("levelsetToDensity",e.what()); return 0; } } static const Pb::Register _RP_levelsetToDensity ("","levelsetToDensity",_W_0);  extern "C" { void PbRegister_levelsetToDensity() { KEEP_UNUSED(_RP_levelsetToDensity); } } 
}; 

