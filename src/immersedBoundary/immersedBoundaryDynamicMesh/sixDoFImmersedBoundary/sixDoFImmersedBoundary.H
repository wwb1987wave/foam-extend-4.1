/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sixDoFImmersedBoundary

Description
    Moving immesed boundary.  Motion is prescribed using a solid body motion
    function.

SourceFiles
    sixDoFImmersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFImmersedBoundary_H
#define sixDoFImmersedBoundary_H

#include "polyMesh.H"
#include "sixDoFRigidBodyMotion.H"
#include "triSurfaceMesh.H"
#include "OFstream.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sixDoFImmersedBoundary Declaration
\*---------------------------------------------------------------------------*/

class sixDoFImmersedBoundary
{
    // Private data

        //- Patch name
        const word name_;

        //- Mesh reference
        const polyMesh& mesh_;

        //- Overset region motion control function
        sixDoFRigidBodyMotion motion_;

        //- Reference tri surface mesh position
        triSurfaceMesh refIbSurface_;
        
        //- Orientation of initial state
        tensor initialQ_;
        
        scalar rhoInf_;
        
        septernion TR_;
        
        //- Forces/moment file ptr
        autoPtr<OFstream> stateFilePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDoFImmersedBoundary(const sixDoFImmersedBoundary&);

        //- Disallow default bitwise assignment
        void operator=(const sixDoFImmersedBoundary&);
        
        const tensor& initialQ() const
        {
        	return initialQ_;
		}
		
		//- If the forces file has not been created create it
        void makeFile();

        //- Output file header information
        void writeFileHeader();      

public:

    // Constructors

        //- Construct from dictionary
        sixDoFImmersedBoundary
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~sixDoFImmersedBoundary();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }
        
        //- Return name
        const septernion& TR() const
        {
            return TR_;
        }

        //- Return mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }
        
                
        //- Write the forces
        virtual void write();  
        
        void transformation();

        //- Move immersed boundary patch points
        void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
