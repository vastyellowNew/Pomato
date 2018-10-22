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


#ifndef __VTK_CAMERA_HELPER_HPP__
#define __VTK_CAMERA_HELPER_HPP__

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <vtkCommand.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>

namespace vtk_utils {

void print_camera_settings(std::ostream& os, vtkRenderer* ren) {
    double pos[3], foc[3], up[3];
    double clip[2];
    ren->GetActiveCamera()->GetPosition(pos);
    ren->GetActiveCamera()->GetFocalPoint(foc);
    ren->GetActiveCamera()->GetViewUp(up);
    ren->GetActiveCamera()->GetClippingRange(clip);
    double scale = ren->GetActiveCamera()->GetParallelScale();
    int* res = ren->GetRenderWindow()->GetSize();

    os << "# Current camera setting\n"
         << "position:       [" << pos[0] << ", " << pos[1] << ", " << pos[2] << "]\n"
         << "focal_point:    [" << foc[0] << ", " << foc[1] << ", " << foc[2] << "]\n"
         << "up_vector:      [" << up[0] << ", " << up[1] << ", " << up[2] << "]\n"
         << "clipping_range: [" << clip[0] << ", " << clip[1] << "]\n"
         << "parallel_scale: " << scale << '\n'
         << "window_resolution: [" << res[0] << ", " << res[1] << "]\n";
}

class camera_setting_callback : public vtkCommand {
public:
    static camera_setting_callback* New() {
        return new camera_setting_callback;
    }
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        vtkRenderer* ren = reinterpret_cast<vtkRenderer*>(caller);
        print_camera_settings(std::cout, ren);
    }
};

void track_camera_setting(vtkRenderer* renderer) {
    VTK_CREATE(camera_setting_callback, callback);
    renderer->AddObserver(vtkCommand::StartEvent, callback);
}

void export_camera_settings(const std::string& filename,
                            vtkRenderer* renderer) {
    std::ostringstream oss;
    print_camera_settings(oss, renderer);

    if (!filename.empty()) {
        std::ofstream fout(filename.c_str());
        if (fout.fail()) {
            throw std::runtime_error("unable to open " + filename);
        }
        fout << oss.str();
        fout.close();
    }
    else std::cout << oss.str() << '\n';
};

void import_camera_settings(const std::string& filename,
                            vtkRenderer* renderer) {
    std::ifstream fin(filename.c_str());
    if (fin.fail()) {
        throw std::runtime_error("unable to open " + filename);
    }

    std::string buffer, item;
    double _array[3];
    double _scale;
    char c;

    try {
        while (!fin.eof()) {
            std::streampos at = fin.tellg();
            std::getline(fin, buffer);
            if (buffer[0] == '#') continue;
            else fin.seekg(at);

            fin >> item;
            if (fin.fail()) break;
            else if (item == "position:") {
                fin >> c >> _array[0] >> c >> _array[1] >> c >> _array[2] >> c;
                renderer->GetActiveCamera()->SetPosition(_array);
            }
            else if (item == "focal_point:") {
                fin >> c >> _array[0] >> c >> _array[1] >> c >> _array[2] >> c;
                renderer->GetActiveCamera()->SetFocalPoint(_array);
            }
            else if (item == "up_vector:") {
                fin >> c >> _array[0] >> c >> _array[1] >> c >> _array[2] >> c;
                renderer->GetActiveCamera()->SetViewUp(_array);
            }
            else if (item == "clipping_range:") {
                fin >> c >> _array[0] >> c >> _array[1] >> c;
                renderer->GetActiveCamera()->SetClippingRange(_array);
            }
            else if (item == "parallel_scale:") {
                fin >> _scale;
                renderer->GetActiveCamera()->SetParallelScale(_scale);
            }
            else if (item == "window_resolution:") {
                int resx, resy;
                fin >> c >> resx >> c >> resy >> c;
                if (renderer->GetRenderWindow() != nullptr) {
                    renderer->GetRenderWindow()->SetSize(resx, resy);
                }
            }
            else {
                std::cerr << "WARNING: unrecognized token: " << item
                          << ". Skipping\n";
                std::getline(fin, buffer);
            }
        }
        fin.close();
    }
    catch(std::runtime_error& e) {
        fin.close();
        std::cerr << "Exception caught while importing camera file:\n"
                  << "\t" << e.what() << '\n';
    }
}

} // vtk_utils

#endif
