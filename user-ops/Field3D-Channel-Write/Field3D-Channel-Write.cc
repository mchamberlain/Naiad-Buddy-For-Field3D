// ----------------------------------------------------------------------------
//
// Field3D-Channel-Write.cc
//
// Copyright (c) 2011 Exotic Matter AB.  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of Exotic Matter AB nor its contributors may be used to
//   endorse or promote products derived from this software without specific
//   prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
//    LIMITED TO,  THE IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS
//    FOR  A  PARTICULAR  PURPOSE  ARE DISCLAIMED.  IN NO EVENT SHALL THE
//    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//    BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE GOODS  OR  SERVICES;
//    LOSS OF USE,  DATA,  OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER
//    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT
//    LIABILITY,  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN
//    ANY  WAY OUT OF THE USE OF  THIS SOFTWARE,  EVEN IF ADVISED OF  THE
//    POSSIBILITY OF SUCH DAMAGE.
//
// ----------------------------------------------------------------------------

// Naiad Base API
#include <NbFilename.h>
#include <NbBlock.h>

// Naiad Graph API
#include <NgBodyOp.h>
#include <NgProjectPath.h>

// Naiad Interface
#include <Ni.h>

// Naiad to Field3D interface
#include <NbF3D.h>

class Field3D_Channel_Write : public Ng::BodyOp
{
public:
    Field3D_Channel_Write(const Nb::String& name)
        : Ng::BodyOp(name) {}
// ----------------------------------------------------------------------------
    virtual Nb::String
    typeName() const
    { return "Field3D-Channel-Write"; }
// ----------------------------------------------------------------------------
    virtual void
    stepAdmittedBody(Nb::Body*             body,
                     Ng::NelContext&       nelContext,
                     const Nb::TimeBundle& tb)
    {
        //Control that the body actually has a field shape
        if (!body->matches("Field"))
            NB_WARNING("No Field Shape found.");

        // Get Bodies String
        const Nb::String bodiesStr = param1s("Body Name(s)")->eval(tb);
        // skip bodies not listed in "Body Names" ...
        if(!body->name().listed_in(bodiesStr))
            return;

        //Initialize Field 3D I/O
        Field3D::initIO();

        //Write output
        Field3D::Field3DOutputFile out;

        //Filename like "path/body_name.#.f3d"
        const Nb::String bodyStr = body->name();
        const Nb::String path = param1s("Output File Path")->eval(tb);
        Nb::String filename = Nb::sequenceToFilename(
            Ng::projectPath(),
            path + Nb::String("/") + bodyStr + Nb::String(".#.f3d"),
            tb.frame,
            tb.timestep,
            param1i("Frame Padding")->eval(tb)
        );

        //Create the actual write output file
        out.create(filename.c_str());

        //Get type (Dense or Sparse)
        const Nb::String typeStr = param1e("Output Type")->eval(tb);
        int type;
        if (typeStr == Nb::String("Dense Field")){
            type = NbF3D::DENSE;
        }
        else if (typeStr == Nb::String("Sparse Field")) {
            type = NbF3D::SPARSE;

            if (!NbF3D::IsPower2(body->constLayout().tileSize())){
                NB_WARNING("Tile Size NOT a power of 2 in " << body->name()
                        << ". Blame Field3D. Skipping body.");
                return;
            }

            //Remove this once Side FX has fixed their plug-in
            NB_WARNING("Houdini's Field3D plug-in tend to crash when loading"
                    << " Sparse voxel fields");
        }

        // Get Channel String
        const Nb::String channelsStr = param1s("Channel(s)")->eval(tb);

        //Get Min and Max @ cell position
        Nb::Vec3i min,max;
        NbF3D::getMinMax(body, min ,max);

        //Get the size. Remember, size in Field3D is from (0,0) to (N,N)
        Nb::Vec3i size = NbF3D::getBBSize(min,max);

        //Get size of each cell
        const float cSize =
                Ng::Store::globalOp()->param1f("Master Cell Size")->eval(tb);

        //Get offset. Should map (-N,-N) to (0,0)
        Nb::Vec3i offset = NbF3D::getCellOffset(min);

        //Gather all properties in a NbF3D struct
        NbF3D::data data = {body, size, offset, cSize, out};

        //Access field shape
        const Nb::FieldShape & field = body->constFieldShape();

        //Loop over all channels and write to Field3D file.
        for (int chnIdx = 0; chnIdx < field.channelCount(); ++chnIdx){
            //Get the name of the current channel
            Nb::String tmpChan = field.channel(chnIdx)->name();

            //Default value in naiad channel. It is possible that future
            //of Naiad will return a def value equal to the type of the channel
            //and if that happens, this code needs some minor fixes.
            const Nb::Vec3f nDef = field.channel(chnIdx)->defaultValue();

            //If the channel is not present in the channel list passed by
            //the user, then skip it
            if (!tmpChan.listed_in(channelsStr))
                continue;

            //Call Naiad to Field3D interface to store data.
            switch (field.channel(chnIdx)->type()){
                case Nb::ValueBase::FloatType: {
                    NbF3D::chanToF3D<true>(type
                            , data
                            , field.constField1f(chnIdx)
                            , chnIdx
                            , nDef[0]);
                    break;
                }
                case Nb::ValueBase::Vec3fType: {
                    //Use Imath in instead of Nb to fit Field3D
                    Imath::V3f def(nDef[0],nDef[1],nDef[2]);

                    //This struct has an operator() that will store
                    //data in the center of each cell (even if field uses
                    // MaC structure.
                    NbF3D::Field3f fieldData = {field.constField3f(chnIdx,0)
                            , field.constField3f(chnIdx,1)
                            , field.constField3f(chnIdx,2)
                            , body->constLayout()};

                    NbF3D::chanToF3D<false>(type
                            , data
                            , fieldData
                            , chnIdx
                            , def);
                    break;
                }
                case Nb::ValueBase::IntType:{
                    NbF3D::chanToF3D<true>(type
                            , data
                            , field.constField1i(chnIdx)
                            , chnIdx
                            , int(nDef[0]));
                    break;
                }
                default: {
                    NB_WARNING("Can't write data from body: " << bodyStr <<
                               " and channel: " << tmpChan <<
                               ". Type not supported (see help)");
                }
            }
        }

        //Close Field3D write file
        out.close();
    }
};
// ----------------------------------------------------------------------------

// Register and upload this user op to the dynamics server

extern "C" Ng::Op*
NiUserOpAlloc(const Nb::String& name)
{
    return new Field3D_Channel_Write(name);
}

// ----------------------------------------------------------------------------

