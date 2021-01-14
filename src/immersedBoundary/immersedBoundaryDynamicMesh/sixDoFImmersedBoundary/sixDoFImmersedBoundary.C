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

\*---------------------------------------------------------------------------*/

#include "sixDoFImmersedBoundary.H"
#include "immersedBoundaryPolyPatch.H"
#include "mixedIbFvPatchFields.H"
#include "foamTime.H"
#include "transformField.H"
#include "forces.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFImmersedBoundary::sixDoFImmersedBoundary
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    motion_(dict),
    initialQ_
    (
        dict.lookupOrDefault("initialOrientation", I)
    ),
    rhoInf_(readScalar(dict.lookup("rhoInf"))),
    stateFilePtr_(nullptr),
    refIbSurface_
    (
        IOobject
        (
            name  + ".ftr",
            mesh.time().constant(),      // instance
            "triSurface",                // local
            mesh,                        // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFImmersedBoundary::~sixDoFImmersedBoundary()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::sixDoFImmersedBoundary::makeFile()
{
    // Create the forces file if not already created
    if (stateFilePtr_.empty())
    {

        // File update
        if (Pstream::master())
        {
            fileName stateDir;
            word startTimeName =
                mesh().time().timeName(mesh().time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                stateDir = mesh().time().path()/".."/name_/startTimeName;
            }
            else
            {
                stateDir = mesh().time().path()/name_/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(stateDir);

            // Open new file at start up
            stateFilePtr_.reset(new OFstream(stateDir/"motionState.dat"));

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::sixDoFImmersedBoundary::writeFileHeader()
{
    if (stateFilePtr_.valid())
    {
        stateFilePtr_()
            << "# Time" << tab
            << "CMx CMy CMz LVx LVy LVz AVx AVy AVz"
            << endl;
    }
}

void Foam::sixDoFImmersedBoundary::write()
{
	// Create the forces file if not already created
	makeFile();
	
	if (Pstream::master())
	{
		stateFilePtr_() << mesh().time().value() << tab << motion_.centreOfMass().x() << tab << motion_.centreOfMass().y() << tab << motion_.centreOfMass().z() << tab << motion_.currentVelocity(motion_.centreOfMass()).x() << tab << motion_.currentVelocity(motion_.centreOfMass()).y() << tab << motion_.currentVelocity(motion_.centreOfMass()).z()  << tab << motion_.omega().x()  << tab << motion_.omega().y()  << tab << motion_.omega().z() << endl;
    }
}

void Foam::sixDoFImmersedBoundary::transformation()
{
	const Time& t = mesh().time();
	
	motion_.updatePosition(t.deltaTValue());

    dictionary forcesDict;
    
    forcesDict.add("patches", wordList(1, name()));
    forcesDict.add("rhoName", "rhoInf");
    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("CofR", motion_.centreOfMass());

    forces f("forces", mesh_, forcesDict);

    forces::forcesMoments fm = f.calcForcesMoment();

    // Get the forces on the patch faces at the current positions

    vector gravity = vector::zero;

    if (mesh_.foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField g =
            mesh_.lookupObject<uniformDimensionedVectorField>("g");

        gravity = g.value();
    }

    motion_.updateForce
    (
        fm.first().first() + fm.first().second() + gravity*motion_.mass(),
        fm.second().first() + fm.second().second(),
        t.deltaTValue()
    );
    
    tensor rotationTensor = motion_.orientation() & initialQ().T();
    
    Info << "Q: " << motion_.orientation() << " " <<  initialQ().T() << endl;

    vector v(vector::zero);
    
    scalar w = 0.0;
    
    scalar trace =
      rotationTensor.xx()
    + rotationTensor.yy()
    + rotationTensor.zz();

    if (trace > 0)
    {
        scalar s = 0.5/Foam::sqrt(trace + 1.0);

        w = 0.25/s;
        v[0] = (rotationTensor.zy() - rotationTensor.yz())*s;
        v[1] = (rotationTensor.xz() - rotationTensor.zx())*s;
        v[2] = (rotationTensor.yx() - rotationTensor.xy())*s;
    }
    else
    {
        if
        (
            rotationTensor.xx() > rotationTensor.yy()
         && rotationTensor.xx() > rotationTensor.zz()
        )
        {
            const scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.xx()
              - rotationTensor.yy()
              - rotationTensor.zz()
            );

            w = (rotationTensor.zy() - rotationTensor.yz())/s;
            v[0] = 0.25*s;
            v[1] = (rotationTensor.xy() + rotationTensor.yx())/s;
            v[2] = (rotationTensor.xz() + rotationTensor.zx())/s;
        }
        else if
        (
            rotationTensor.yy() > rotationTensor.zz()
        )
        {
            const scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.yy()
              - rotationTensor.xx()
              - rotationTensor.zz()
            );

            w = (rotationTensor.xz() - rotationTensor.zx())/s;
            v[0] = (rotationTensor.xy() + rotationTensor.yx())/s;
            v[1] = 0.25*s;
            v[2] = (rotationTensor.yz() + rotationTensor.zy())/s;
        }
        else
        {
            const scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.zz()
              - rotationTensor.xx()
              - rotationTensor.yy()
            );

            w = (rotationTensor.yx() - rotationTensor.xy())/s;
            v[0] = (rotationTensor.xz() + rotationTensor.zx())/s;
            v[1] = (rotationTensor.yz() + rotationTensor.zy())/s;
            v[2] = 0.25*s;
        }
    }
    
    Foam::septernion TR(motion_.centreOfMass(), quaternion(w,v));
    
    TR_=TR;
    
    write();
}

void Foam::sixDoFImmersedBoundary::movePoints()
{
	transformation();
	
    // Get ibMesh from patch
    const label patchID = mesh().boundaryMesh().findPatchID(name());

    if (patchID < 0)
    {
        FatalErrorIn
        (
            "void sixDoFImmersedBoundary::movePoints() const"
        )   << "Patch " << name() << " not found.  Available patch names: "
            << mesh().boundaryMesh().names()
            << abort(FatalError);
    }

    const immersedBoundaryPolyPatch& cibPatch =
        refCast<const immersedBoundaryPolyPatch>
        (
            mesh().boundaryMesh()[patchID]
        );

    // Get non-const reference to patch
    immersedBoundaryPolyPatch& ibPatch =
        const_cast<immersedBoundaryPolyPatch&>(cibPatch);

    // Move points
    ibPatch.moveTriSurfacePoints
    (
        transform(TR_, refIbSurface_.points())
    );
}


// ************************************************************************* //
