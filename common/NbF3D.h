// ----------------------------------------------------------------------------
//
// NbF3D.h
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

//Field 3D API
#include <Field3D/DenseField.h>
#include <Field3D/SparseField.h>
#include <Field3D/Field3DFile.h>
#include <Field3D/InitIO.h>

namespace NbF3D{
// ----------------------------------------------------------------------------
    enum type{
        DENSE = 0,
        SPARSE = 1,
        MAC = 2 //Currently not supported
    };
// ----------------------------------------------------------------------------
    struct data
    {
        const Nb::Body * body;
        const Nb::Vec3i & size;
        const Nb::Vec3i & offset;
        const float cSize;
        Field3D::Field3DOutputFile & out;
    };
// ----------------------------------------------------------------------------
    //Not used at the moment
    class Vec3f : public Imath::V3f
    {
    public:
        Vec3f(): Imath::V3f() {};
        Vec3f(int i): Imath::V3f(i) {};
        Vec3f(const Nb::Vec3f & v) : Imath::V3f(v[0], v[1], v[2]){};
    };
// ----------------------------------------------------------------------------
    struct Field3f{
        const Nb::Field1f &  fU, & fV, & fW;
        const Nb::TileLayout & tileLayout;

        Imath::V3f
        operator()(const int i) const
        {
            //Get cell indices from global tile cell index i
            int ci, cj, ck;
            tileLayout.cijk(i, ci, cj, ck);

            //Get the global coord of the other side of the MAC
            //Will return -1 if outside fine tile layout
            int       idu = tileLayout.cell(ci - 1, cj, ck)
                    , idv = tileLayout.cell(ci, cj - 1, ck)
                    , idw = tileLayout.cell(ci, cj, ck - 1);

            //Float values (for example velocity)
            float u,v,w;

            //If negative index, clamp to the only valid one.
            //If positive, mean value of the two Marker values
            //(store in center ofCell)
            if (idu < 0)
                u = fU(i);
            else
                u = (fU(i) + fU(idu))/2;

            if (idv < 0)
                v = fV(i);
            else
                v = (fV(i) + fV(idv))/2;

            if (idw < 0)
                w = fW(i);
            else
                w = (fW(i) + fW(idw))/2;

            //Imath to be safe.
            return Imath::V3f(u,v,w);
        };
    };
// ----------------------------------------------------------------------------
    template<class F, class T> void
    initField(      F fieldPtr,
                const data & d,
              const int chnIdx,
                 const T & def)
    {
        //Field name after the body's name
        fieldPtr->name = d.body->name().c_str();

        //Access field shape
        const Nb::FieldShape & fieldShape = d.body->constFieldShape();

        //Attribute name after the channel name
        fieldPtr->attribute = fieldShape.channel(chnIdx)->name().c_str();

        //Bounding Box
        fieldPtr->setSize(Imath::V3i(d.size[0], d.size[1], d.size[2]));

        //Set default value everywhere to begin with
        fieldPtr->clear(def);

        //Move to the correct world coordinates
        Field3D::M44d localToWorld;

        //Scale the size according to Master Cell Size
        localToWorld.setScale(
                Imath::V3d(d.size[0], d.size[1], d.size[2]) * d.cSize);

        //We don't want the top bottom left corner to be at (0,0,0)
        localToWorld *= Field3D::M44d().setTranslation(
                Imath::V3d(d.offset[0], d.offset[1], d.offset[2]) * d.cSize);

        //Apply Matrix
        Field3D::MatrixFieldMapping::Ptr mapping(
                new Field3D::MatrixFieldMapping);
        mapping->setLocalToWorld(localToWorld);
        fieldPtr->setMapping(mapping);
    };
// ----------------------------------------------------------------------------
    template<class F, class nF> void
    chanToField(F            fieldPtr,
                const data &        d,
                const nF & naiadField)
    {
        //Access field shape
        const Nb::FieldShape & fieldShape = d.body->constFieldShape();

        //Get Tile Layout
        const Nb::TileLayout & tileLayout = d.body->constLayout();

        //Cell indices
        int ci, cj ,ck;

#pragma omp parallel for schedule(dynamic)
        for (int iTiles = 0; iTiles < tileLayout.fineTileCount(); ++iTiles){
            Nb::Tile tile(tileLayout.fineTile(iTiles));
            //We are only interested in the fine tiles.
            if(tile.type()==Nb::Tile::Fine) {
                for (int iTC = tile.cellsBegin(); iTC < tile.cellsEnd(); ++iTC){
                    tileLayout.cijk(iTC,ci,cj,ck);
                    ci -= d.offset[0];
                    cj -= d.offset[1];
                    ck -= d.offset[2];
                    fieldPtr->fastLValue(ci, cj, ck) = naiadField(iTC);
                }
            }
        }
    };
// ----------------------------------------------------------------------------
    template <bool scalar, class F, class T>
    class F3DFileWrite{
    public:
        static void
        write(Field3D::Field3DOutputFile &  out, T fieldPtr)
        {
            out.writeScalarLayer<F>(fieldPtr);
        };
    };
// ----------------------------------------------------------------------------
    template<class F, class T>
    class F3DFileWrite<false, F, T>
    {
    public:
        static void
        write(Field3D::Field3DOutputFile &  out, T fieldPtr)
        {
            out.writeVectorLayer<float>(fieldPtr);
        };
    };
// ----------------------------------------------------------------------------
    template <bool scalar, class T, class nF> void
    chanToDenseF3D(const data &        d,
                   const nF & naiadField,
                   const int      chnIdx,
                   const T &         def)
    {
        //Create a Dense Field3D Voxel Field
        typename Field3D::DenseField<T>::Ptr fieldPtr(
                new Field3D::DenseField<T>());

        //Init names and size
        initField(fieldPtr, d, chnIdx, def);

        //From Naiad field to Field3D field
        chanToField(fieldPtr, d, naiadField);

        //From Dense Field3D field to Field3D file
        F3DFileWrite<scalar, T, typename Field3D::DenseField<T>::Ptr>::write(
                d.out,fieldPtr);
    };
// ----------------------------------------------------------------------------
    template <bool scalar, class T, class nF> void
    chanToSparseF3D(const data &        d,
                    const nF & naiadField,
                    const int      chnIdx,
                    const T &         def)
    {
        //Create a Sparse Field3D Voxel Field
        typename Field3D::SparseField<T>::Ptr fieldPtr(
                new Field3D::SparseField<T>());

        //Init names and size
        initField(fieldPtr, d, chnIdx, def);

        //From Naiad field to Field3D field
        chanToField(fieldPtr, d, naiadField);

        //From Sparse Field3D field to Field3D file
        F3DFileWrite<scalar, T, typename Field3D::SparseField<T>::Ptr>::write(
                d.out,fieldPtr);
    };
// ----------------------------------------------------------------------------
    template <bool scalar, class T, class nF> void
    chanToF3D(const int        type,
              const data &        d,
              const nF & naiadField,
              const int      chnIdx,
              const T &         def)
    {
        if (type == DENSE){
            chanToDenseF3D<scalar, T>(d, naiadField, chnIdx, def);
        }
        else if (type == SPARSE){
            chanToSparseF3D<scalar, T>(d, naiadField, chnIdx, def);
        }
    };

    Nb::Vec3i
    getBBSize(const Nb::Vec3i & min, const Nb::Vec3i & max)
    {
        //Important to add a vector (1,1,1)
        return max - min + Nb::Vec3i(1,1,1);
    };
// ----------------------------------------------------------------------------
    Nb::Vec3i
    getCellOffset(const Nb::Vec3i & min)
    {
        return min;
    };
// ----------------------------------------------------------------------------
    void
    getMinMax(const Nb::Body * body,
              Nb::Vec3i & min,
              Nb::Vec3i & max)
    {
        //Cell indices
        int ci,cj,ck;

        //Get Tile Layout
        const Nb::TileLayout & tileLayout = body->constLayout();

        //MinMax values for int
        const int maxVal(std::numeric_limits<int>::max());
        const int minVal(std::numeric_limits<int>::min());

        //Set Vectors to opposite possible value.
        min = Nb::Vec3i(maxVal,maxVal,maxVal);
        max = Nb::Vec3i(minVal,minVal,minVal);

        //Loop all fine tiles
        for (int iTiles = 0; iTiles < tileLayout.fineTileCount(); ++iTiles){
            Nb::Tile tile(tileLayout.fineTile(iTiles));

            //Safety check
            if(tile.type()==Nb::Tile::Fine) {
                //Loop over all cells in the tile
                for (int iCells = tile.cellsBegin();
                        iCells < tile.cellsEnd();
                        ++iCells){

                    //Get world cell position (naiad can have cells at -1 etc)
                    tileLayout.cijk(iCells,ci,cj,ck);

                    //Perform min and max checks
                    if (ci < min[0])
                        min[0] = ci;
                    if (ci > max[0])
                        max[0] = ci;
                    if (cj < min[1])
                        min[1] = cj;
                    if (cj > max[1])
                        max[1] = cj;
                    if (ck < min[2])
                        min[2] = ck;
                    if (ck > max[2])
                        max[2] = ck;
                }
            }
        }
    };
// ----------------------------------------------------------------------------
    bool
    IsPower2(int x)
    {
        return ( (x > 0) && ((x & (x - 1)) == 0) );
    };
} // end namespace NbF3D
