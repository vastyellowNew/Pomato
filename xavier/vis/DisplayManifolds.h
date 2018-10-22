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


///////////////////////////////////////////////////////////////////////
// CR3BP Manifold Test Header File
//    - A compute module for topological analysis of CR3BP Planar Maps
//    - Specifically tests the manifold propagation for a given saddle
// Author:  Wayne Schlei (Purdue University)
// Date:  3/27/2016
///////////////////////////////////////////////////////////////////////

#ifndef DISPLAY_MANIFOLD_H
#define DISPLAY_MANIFOLD_H

//Base Definition
#include <hxcore/HxModule.h>
//Graphics
#include <mclib/McHandle.h> //smart pointer template class
#include <Inventor/nodes/SoSeparator.h> //Scene object separator
//Ports
#include <hxcore/HxPortToggleList.h> // Allows for Check-boxes
#include <hxcore/HxPortButtonList.h> // Allows for button list
#include <hxtime/HxPortTime.h> // Allows for Time slider
#include <hxcore/HxPortFloatTextN.h> // Allows for the Float Text Boxes
#include <hxcore/HxPortFloatSlider.h> // Allows for a float slider
#include <hxcore/HxPortColorList.h> // Allows for a color list
#include <hxcolor/HxPortColormap.h> // Colormaps
#include <hxcore/HxPortMultiMenu.h> // Allows for a multi-menu
#include <hxcore/HxPortSeparator.h> // Separator port
#include <hxcore/HxPortInfo.h> // Info port
#include <hxcore/HxPortRadioBox.h> // Radio Boxes
#include <hxcore/HxPortIntTextN.h> // Allows for Integer Text Box
#include <hxcore/HxPortIntSlider.h> // Allows for Integer Slider
#include <hxcore/HxPortText.h>     // Text Box (QString)
#include <hxcore/HxPortTabBar.h>   // Tab bar
#include <hxcore/HxPortGeneric.h>  // Generic Port
#include <hxcore/HxConnection.h>   //Extra connections for input 
#include <hxcore/HxPortFontSelection.h> //Fonts

//Other Data Types
#include <hxcluster/internal/HxCluster.h> 
#include <hxlines/internal/HxLineSet.h>   
#include <mypackage/spaceflight/State.h>
#include <mypackage/CR3BPsystem.h>
#include <mypackage/HxTrajectory.h>
//#include <poincare/HxManifoldArc.h> //Maybe need special HxTrajectory version
#include <poincare/MapSpaceProjection.h>
#include <poincare/ManifoldData.h>
#include <poincare/WWDesign.h>
#include <poincare/poincareAPI.h>

//DisplayManifolds Class
class POINCARE_API DisplayManifolds : public HxModule
{
    HX_HEADER(DisplayManifolds);
    
    /// Comparison class for (lineIdx,SegIdx) pair in SoLineSet
    struct LineSegPair {
        int lineIdx, segIdx;
        ///Constructor
        LineSegPair(const int lid, const int sid) :
            lineIdx(lid), segIdx(sid) {}
        /// Exactly equal operator
        bool operator==(const LineSegPair& other) const {
          return (lineIdx == other.lineIdx && segIdx == other.segIdx);
        }
        /// Less than operator
        bool operator<(const LineSegPair& other) const {
          return (lineIdx < other.lineIdx ||
                 (lineIdx == other.lineIdx && segIdx < other.segIdx) );
        }
    };

  public:
    //Constructor
    DisplayManifolds();
    //Destructor
    ~DisplayManifolds();
    //Compute
    virtual void compute();

    //Avizo Ports
    HxPortInfo info;
    HxConnection initStateConnection;
    HxPortTabBar portTabs;
    HxPortMultiMenu system;
    HxPortText gravParam;
    HxPortText portJacobiConstant;
    HxPortFloatTextN portXbounds;
    HxPortFloatTextN portYbounds;
    HxPortSeparator sep0;
    HxPortToggleList displayOptions;
    HxPortMultiMenu  horzMenu;
    HxPortMultiMenu  vertMenu;
    HxPortFloatSlider horzScale;
    HxPortFloatSlider vertScale;

    //Object display options
    HxPortFloatSlider eigenLength;
    HxPortSeparator sep1;
    HxPortInfo dispInfo;
    HxPortColorList colors;
    HxPortColormap  stableColormap;
    HxPortColormap  unstableColormap;
    HxPortIntSlider portLineWidth;
    HxPortIntSlider portPointSize;
    HxPortFloatSlider portMarkerScale;
    HxPortSeparator sep2;
    HxPortInfo  filterInfo;
    HxPortToggleList filterOptions;
    HxPortRadioBox filterType;
    HxPortIntSlider filterOrbitID;
    HxPortIntSlider filterFixedPoint;
    
    //Selection
    HxPortInfo selectionInfo;
    HxPortToggleList simOptions;
    HxPortRadioBox simDirection;
    HxPortRadioBox simRevColor;
    HxPortButtonList simData;
    HxPortGeneric  hcOutput;
    //Debug elements
    HxPortIntTextN selectedParts;
    HxPortIntTextN selectedManSeg;
    HxPortIntTextN selectedPropInfo;
    HxPortText selectedParamText;
    HxPortIntTextN hcPartnerParts;
    HxPortIntTextN hcPartnerManSeg;
    HxPortSeparator sep3;
    HxPortInfo  returnInfo;
    HxPortToggleList returnDisplayOpts;
    HxPortFloatSlider returnMarkerScale;
    
    //Function listing
    HxPortInfo trajOutInfo;
    HxPortInfo trajOutInfo2;
    HxPortInfo returnsOutInfo;
    HxPortInfo hcOutInfo;
    HxPortSeparator sep4;
    HxPortInfo funcInfo;
    HxPortIntSlider  funcOrbitID;
    HxPortIntTextN   funcManifoldID;
    HxPortButtonList clearFuncs;
    HxPortText  portOrbitIDList;
    HxPortButtonList extractFuncs;
    
    //Topology design ports
    HxPortInfo          portDesignInfo;
    HxPortText          jcIS;
    HxPortMultiMenu     portDesignOpt;
    HxPortToggleList    portDesignDisplayOpts;
    HxPortRadioBox      portDeltaVUnits;
    HxPortColorList     designColors;
    HxPortFontSelection portDesignFont;
    HxPortButtonList    portDesignOutput;
    HxPortButtonList    portWtoWDesignOutput;

    
    /// The picked point
    SbVec3f pickedPoint;
    /// Selected Line Index
    int selectedLineID;
    /// Selected segment (part) index
    int selectedPartID;
    /// Linear parameter along a segment based on selection point 
    double selectedLinearParameter;
    /// Bool indicating if selection is unstable manifold
    bool selectedUnstable;
    /// Selected Line Index of a Heteroclinic/homoclinic partner
    int hcPartnerLineID;
    /// Selected segment (part) index of a Heteroclinic/homoclinic partner 
    int hcPartnerPartID;
    /// Function to construct the heteroclinic/homoclinic data pairing (returns true if valid HC/HO)
    bool makeHHConnection();
    /// Heteroclinic (or homoclinic) connection information [ALWAYS Wu to Ws]
    WWDesign hhConn;
    /// Linear parameter along a segment based on selection point 
    //double hcPartnerLinearParameter;
    /// Flag to show integrated arcs 
    bool show;
    /// Flag indicating whether or not to perform numerical propagation of a selection
    bool integrate;
    /// Progenitor State point (UNUSED!)
    SbVec3f progStatePoint; //In (x,y,z) space
    
    /// Pointer to appropriate SoLineSet Node in a target scene (for reference with pick)
    SoNode *targetNode,*wuTarget,*wsTarget;
    
    /// InitialState input is connected
    bool hasInput;
    /// Flag to indicate if the button is pressed
    bool buttonPressed;
    /// Target node for deltaV line
    SoNode *topoDesignTargetNode;
    /// Function for rendering the topology design scene 
    void renderTopoDesignScene();
    /// The picked velocity differential from topology design
    SbVec3d deltaV;
    /// Function to translate selected point to topology design maneuver (Only with hasInput=true)
    bool topologyDeltaV();
    /// Flag to indicate that a topology design solution was computed
    bool topoDesignSelected;
    /// Topology design flag (indicating that manifold was selected for design)
    bool wwTopoDesignON;
    /// Design object for a Manifold to Manifold Design creation
    WWDesign  augHConn;
    /// Render the DeltaV constraint line for Manifold-to-Manifold design following click actions (also removes)
    //void renderWWDeltaVConstraint(bool remove);
    /// Render the partially constructed design to indicate selection location (continuously altered)
    void renderWWTempScene(bool remove);
    //void renderWWDeltaVLine(bool remove);
    /// Render the Manifold-to-Manifold Design scene for a valid computed design
    void renderWWTopoDesignScene(bool remove);
    /// Initiating a Manifold-to-Manifold Design by storing information with starting selection
    void startAugmentedConnection();
    /// Compute the Manifold-to-Manifold Design based on a valid 'hover-over' a secondary manifold
    bool computeAugmentedConnection(); //Returns true if a valid selection (intersects DeltaVLine)
    
    
    /// Integrate A Manifold Arc from a selection
    void integrateManifoldArc(const bool useLinearParam=false);
    ///Render the manifold objects
    void renderManifolds();
    ///Perform the propagation, rendering and data output if necessary
    void arcSimulation(const bool useLinearParam=false);

  protected:
     typedef std::pair<int,int> IntPair;
     McHandle<SoSeparator> scene;
     McHandle<SoSeparator> trajScene;
     McHandle<SoSeparator> tdScene;
     McHandle<SoSeparator> wwdvScene;
     McHandle<SoEventCallback> eventCB;
     bool initCB;
     SoSeparator* root;
     ///Hash Table for (lineIdx,segIdx) -> (manifoldID,segID) on Unstable Manifold
     std::map<LineSegPair,IntPair> uLinesetToManifold;
     ///Hash Table for (lineIdx,segIdx) -> (manifoldID,segID) on Stable Manifold
     std::map<LineSegPair,IntPair> sLinesetToManifold;
  private:
     CR3BPsystem *sys;
     MapSpaceProjection mapSpaceProjector;
     double mup, C, lstar;
     std::vector<State> uStates, dStates;
     std::vector<nvis::vec2> uIterates, uSubIts, dIterates, dSubIts;
     std::vector<double> uTimes, dTimes, uItTimes, dItTimes, uSubItsT, dSubItsT;
     /// List of orbit IDs for extraction
     std::list<int> orbitIDList;
     /// Print the orbitID list to the list text port
     void printListToPort();
     ///Output selected manifold arc as HxLineSet
     void exportArcToLineSet(HxLineSet *line);
     ///Output selected manifold arc as HxTrajectory object
     void exportArcToTrajectory(HxTrajectory *traj);
     ///Output returns to a HxCluster object 
     void exportReturnsToCluster(HxCluster *cluster);
     ///Output selected manifold arc as HxManifoldArc object
     //void exportArcToManifoldArc(HxManifoldArc *traj) const;
     
};

#endif // DISPLAY_MANIFOLD_H
