/*************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body 
        problem for orbital mechanics applications. 

Authors: Wayne Schlei and Xavier Tricoche 

Copyright (c) 2013-2018, Purdue University 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/


#ifndef __VTK_IMAGE_HELPER_HPP__
#define __VTK_IMAGE_HELPER_HPP__

#include "vtk_macros.hpp"
#include "vtk_data_helper.hpp"

#ifdef __HAS_EXPRTK__
#include <exprtk.hpp>
#endif

#include <string>

#include <math/fixed_vector.hpp>

#include <vtkAxesActor.h>
#include <vtkBMPReader.h>
#include <vtkBMPWriter.h>
#include <vtkCaptionActor2D.h>
#include <vtkDataArray.h>
#include <vtkDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkJPEGReader.h>
#include <vtkJPEGWriter.h>
#if VTK_MAJOR_VERSION >= 6
#include <vtkNrrdReader.h>
#endif
#include <vtkOrientationMarkerWidget.h>
#include <vtkPNGReader.h>
#include <vtkPNGWriter.h>
#include <vtkPNMReader.h>
#include <vtkPNMWriter.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTexture.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkWindowToImageFilter.h>

namespace vtk_utils {

template<typename Reader>
inline vtkImageData* load_image_of_given_format(const std::string& filename)
{
    VTK_CREATE(Reader, reader);
    reader->SetFileName(filename.c_str());
    reader->Update();
    vtkImageData* _img = vtkImageData::SafeDownCast(reader->GetOutput());
    reader->GetOutput()->Register(_img);
    return _img;
    reader->Delete();
}

template<typename Writer>
inline void save_image_in_given_format(const vtkImageData* data,
                                       const std::string& filename,
                                       int quality=100)
{
    VTK_CREATE(Writer, writer);
    VTK_CONNECT(writer, const_cast<vtkImageData*>(data));
    writer->SetFileName(filename.c_str());
    writer->Write();
}

template<>
inline void save_image_in_given_format<vtkJPEGWriter>(const vtkImageData* data,
                                                      const std::string& filename, 
                                                      int quality) {
    VTK_CREATE(vtkJPEGWriter, writer);
    VTK_CONNECT(writer, const_cast<vtkImageData*>(data));
    writer->SetFileName(filename.c_str());
    writer->SetQuality(quality);
    writer->Write();
}

inline vtkImageData* load_image(const std::string& filename)
{
    std::string ext = xavier::filename::extension(filename);
    if (ext == "png") {
        return load_image_of_given_format<vtkPNGReader>(filename);
    } else if (ext == "tif" || ext == "tiff") {
        return load_image_of_given_format<vtkTIFFReader>(filename);
    } else if (ext == "bmp") {
        return load_image_of_given_format<vtkBMPReader>(filename);
    } else if (ext == "jpg" || ext == "jpeg" || ext == "jpg2" || ext == "jpeg2") {
        return load_image_of_given_format<vtkJPEGReader>(filename);
    }
#if VTK_MAJOR_VERSION >= 6
    else if (ext == "nrrd" || ext == "nhdr") {
        return load_image_of_given_format<vtkNrrdReader>(filename);
    }
#endif
    else if (ext == "pnm" || ext == "ppm") {
        return load_image_of_given_format<vtkPNMReader>(filename);
    } else {
        std::cerr << "image type of " << filename << " (" << ext << ") unrecognized\n";
        throw std::runtime_error("unrecognized image filename extension");
    }
}

inline void save_image(const vtkImageData* data, const std::string& filename, int quality=100)
{
    std::string ext = xavier::filename::extension(filename);
    if (ext == "png") {
        VTK_CREATE(vtkImageShiftScale, cast);
        VTK_CONNECT(cast, const_cast<vtkImageData*>(data));
        double* range = const_cast<vtkImageData*>(data)->GetPointData()->GetScalars()->GetRange();
        double _min = range[0];
        double _max = range[1];
        cast->SetShift(-_min);
        cast->SetScale(255./(_max-_min));  
        cast->SetOutputScalarTypeToUnsignedChar();
        cast->Update();
        save_image_in_given_format<vtkPNGWriter>(cast->GetOutput(), filename);
    } else if (ext == "tif" || ext == "tiff") {
        save_image_in_given_format<vtkTIFFWriter>(data, filename);
    } else if (ext == "bmp") {
        save_image_in_given_format<vtkBMPWriter>(data, filename);
    } else if (ext == "jpg" || ext == "jpeg" || ext == "jpg2" || ext == "jpeg2") {
        save_image_in_given_format<vtkJPEGWriter>(data, filename, quality);
    }
#if __VTK_NRRD_WRITER_WORKS__
    else if (ext == "nrrd" || ext == "nhdr") {
        save_image_in_given_format<vtkNrrdWriter>(data, filename);
    }
#endif
    else if (ext == "pnm" || ext == "ppm") {
        save_image_in_given_format<vtkPNMWriter>(data, filename);
    } else if (ext == "vtk") {
        save_image_in_given_format<vtkDataSetWriter>(data, filename);
    } else if (ext.empty()) {
        save_image_in_given_format<vtkPNGWriter>(data, filename + ".png");
    } else {
        std::cerr << "image type of " << filename << " (" << ext << ") unrecognized\n";
        throw std::runtime_error("unrecognized image filename extension");
    }
}

inline void save_frame(vtkRenderWindow* window, 
                       const std::string& filename,
                       int quality=100)
{
    // turn screen (frame buffer) content into an image
    VTK_CREATE(vtkWindowToImageFilter, capture);
    capture->SetInput(window);
    capture->Update();
    save_image(capture->GetOutput(), filename, quality);
}

inline vtkRenderer* fill_window(vtkRenderer* inout, const nvis::bbox2& bounds) {
    inout->GetActiveCamera()->SetParallelProjection(1);
    nvis::vec2 center=bounds.center();
    inout->GetActiveCamera()->SetPosition(center[0], center[1], 1);
    inout->GetActiveCamera()->SetFocalPoint(center[0], center[1], 0);
    inout->GetActiveCamera()->SetViewUp(0, 1, 0);
    inout->GetActiveCamera()->SetParallelScale(0.5*(bounds.max()[1]-bounds.min()[1]));
    return inout;
}

} // vtk_utils

#endif