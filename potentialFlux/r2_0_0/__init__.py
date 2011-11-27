#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#----------------------------------------------------------------------------
from Foam import ref, man

#----------------------------------------------------------------------------
def readControls( mesh ):

    potentialFlow = mesh.solutionDict().subDict( ref.word( "potentialFlow" ) )

    nNonOrthCorr = potentialFlow.lookupOrDefault( ref.word( "nNonOrthogonalCorrectors" ), 0 )
  
    return potentialFlow, nNonOrthCorr


#----------------------------------------------------------------------------
def _createFields( runTime, mesh, potentialFlow ):
       
    ref.ext_Info() << "Reading field p\n" << ref.nl
    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.NO_WRITE ),
                            mesh )
    
    p << ref.dimensionedScalar( ref.word( "zero" ), p.dimensions(), 0.0 )

    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    U <<  ref.dimensionedVector( ref.word( "0" ), U.dimensions(), ref.vector.zero )

    phi = man.surfaceScalarField( man.IOobject( ref.word( "phi" ),
                                                ref.fileName( runTime.timeName() ),
                                                mesh,
                                                ref.IOobject.NO_READ,
                                                ref.IOobject.AUTO_WRITE ),
                                  man.fvc.interpolate( U ) & man.surfaceVectorField( mesh.Sf(), man.Deps( mesh ) ) )


    pRefCell = 0
    pRefValue = 0.0

    pRefCell, pRefValue = ref.setRefCell( p, potentialFlow, pRefCell, pRefValue )

    return p, U, phi, pRefCell, pRefValue


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    ref.argList.validOptions.fget().insert( ref.word( "writep" ), "" )
    
    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )
    
    mesh = man.createMesh( runTime )
    
    potentialFlow, nNonOrthCorr = readControls( mesh )

    p, U, phi, pRefCell, pRefValue = _createFields( runTime, mesh, potentialFlow )

    ref.ext_Info() << ref.nl << "Calculating potential flow" << ref.nl
    
    # Since solver contains no time loop it would never execute
    # function objects so do it ourselves.
    runTime.functionObjects().start()
    
    ref.adjustPhi(phi, U, p)
    
    for nonOrth in range( nNonOrthCorr + 1):
        pEqn = ( ref.fvm.laplacian( ref.dimensionedScalar( ref.word( "1" ), ref.dimTime / p.dimensions() * ref.dimensionSet( 0.0, 2.0, -2.0, 0.0, 0.0 ), 1.0 ), p ) 
                == ref.fvc.div( phi ) )
        
        pEqn.setReference( pRefCell, pRefValue )
        pEqn.solve()

        if nonOrth == nNonOrthCorr:
           phi -= pEqn.flux()
           pass
        pass
    
    ref.ext_Info() << "continuity error = " << ref.fvc.div( phi ).mag().weightedAverage( mesh.V() ).value() << ref.nl

    U << ref.fvc.reconstruct( phi )
    U.correctBoundaryConditions()
    ref.ext_Info() << "Interpolated U error = " << ( ( ( ref.fvc.interpolate( U ) & mesh.Sf() ) - phi ).sqr().sum().sqrt()  /mesh.magSf().sum() ).value() << ref.nl

    # Force the write
    U.write()
    phi.write()
    
    if args.optionFound( ref.word( "writep" ) ):
       p.write()
       pass
       
    runTime.functionObjects().end()
    
    ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 or higher\n "
   pass


#--------------------------------------------------------------------------------------

